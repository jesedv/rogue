use crate::field::{Field1D, Field2D, C};
use crate::random::seeded;
use rand::Rng;
use std::f64::consts::PI;

/// Wrap a distance onto the periodic interval `[−L/2, L/2)` (nearest image).
#[inline]
fn wrap(x: f64, l: f64) -> f64 {
    let w = x.rem_euclid(l);
    if w >= l * 0.5 {
        w - l
    } else {
        w
    }
}

/// Plane wave `ψ = amp · exp(i k0 x)`.
pub fn plane_wave(nx: usize, lx: f64, k0: f64, amp: f64) -> Field1D {
    Field1D {
        data: (0..nx)
            .map(|i| C::from_polar(amp, k0 * i as f64 * lx / nx as f64))
            .collect(),
        nx,
        lx,
    }
}

/// Gaussian hump in 1D.
pub fn gaussian_1d(nx: usize, lx: f64, amp: f64, sigma: f64, center: f64) -> Field1D {
    Field1D {
        data: (0..nx)
            .map(|i| {
                let x = i as f64 * lx / nx as f64;
                let dx = wrap(x - center, lx);
                C::new(amp * (-dx * dx / (2.0 * sigma * sigma)).exp(), 0.0)
            })
            .collect(),
        nx,
        lx,
    }
}

/// Gaussian hump in 2D.
#[allow(clippy::too_many_arguments)]
pub fn gaussian_2d(
    nx: usize,
    ny: usize,
    lx: f64,
    ly: f64,
    amp: f64,
    sigma_x: f64,
    sigma_y: f64,
    cx: f64,
    cy: f64,
) -> Field2D {
    let mut data = vec![C::new(0.0, 0.0); nx * ny];
    for j in 0..ny {
        let y = j as f64 * ly / ny as f64;
        let dy = wrap(y - cy, ly);
        for i in 0..nx {
            let x = i as f64 * lx / nx as f64;
            let dx = wrap(x - cx, lx);
            let g = (-dx * dx / (2.0 * sigma_x * sigma_x) - dy * dy / (2.0 * sigma_y * sigma_y)).exp();
            data[j * nx + i] = C::new(amp * g, 0.0);
        }
    }
    Field2D { data, nx, ny, lx, ly }
}

/// Bright soliton of the focusing cubic NLS, `ψ = A sech(A(x−x0)) e^{i v x}`.
pub fn soliton_1d(nx: usize, lx: f64, amp: f64, vel: f64, center: f64) -> Field1D {
    // Exact soliton of iψ_t + ψ_xx + |ψ|²ψ = 0: ψ = √2·A·sech(A·(x−x0))·e^{ivx}.
    Field1D {
        data: (0..nx)
            .map(|i| {
                let x = i as f64 * lx / nx as f64;
                let dx = wrap(x - center, lx);
                let s = std::f64::consts::SQRT_2 * amp / (amp * dx).cosh();
                C::from_polar(s, vel * dx)
            })
            .collect(),
        nx,
        lx,
    }
}

/// Stokes wave train plus Benjamin–Feir sidebands:
/// `ψ = A0 e^{i k0 x} + ε A0 (e^{i(k0−ν)x} + e^{i(k0+ν)x})`, seeded phases.
pub fn stokes_sidebands(
    nx: usize,
    lx: f64,
    amp: f64,
    k0: f64,
    nu: f64,
    eps: f64,
    seed: u64,
) -> Field1D {
    let mut rng = seeded(seed);
    let mut data = vec![C::new(0.0, 0.0); nx];
    for i in 0..nx {
        let x = i as f64 * lx / nx as f64;
        let ph0: f64 = rng.gen_range(0.0..2.0 * PI);
        let phs: f64 = rng.gen_range(0.0..2.0 * PI);
        data[i] = amp * C::from_polar(1.0, k0 * x + ph0)
            + amp * eps * (C::from_polar(1.0, (k0 - nu) * x + phs) + C::from_polar(1.0, (k0 + nu) * x + phs));
    }
    Field1D { data, nx, lx }
}

/// Zero-mass perturbation on a uniform background: used to seed breather
/// evolutions. `ψ = A + ρ(x)` where ρ is a compact bump.
pub fn background_plus_bump(
    nx: usize,
    lx: f64,
    background: f64,
    amp: f64,
    sigma: f64,
    center: f64,
) -> Field1D {
    Field1D {
        data: (0..nx)
            .map(|i| {
                let x = i as f64 * lx / nx as f64;
                let dx = wrap(x - center, lx);
                C::new(background + amp * (-dx * dx / (2.0 * sigma * sigma)).exp(), 0.0)
            })
            .collect(),
        nx,
        lx,
    }
}

/// 2D JONSWAP-style random sea envelope (seeded). `|ψ|` is the envelope
/// amplitude; scaled so its RMS matches the narrow-band expectation
/// `rms(|ψ|) ≈ √2 · Hs/4`.
pub fn random_sea_2d(nx: usize, ny: usize, lx: f64, hs: f64, seed: u64) -> Field2D {
    let mut rng = seeded(seed);
    let dk = 2.0 * PI / lx;
    let mut data = vec![C::new(0.0, 0.0); nx * ny];
    let nmodes = 600;
    for _ in 0..nmodes {
        let kx = rng.gen_range(-4.0..4.0) * dk;
        let ky = rng.gen_range(-4.0..4.0) * dk;
        let k = (kx * kx + ky * ky).sqrt().max(1e-6);
        let f = (9.81 * k).sqrt() / (2.0 * PI);
        let fp = 0.15;
        let f_r = f / fp;
        let sigma = if f <= fp { 0.07 } else { 0.09 };
        let peak = (-((f - fp).powi(2)) / (2.0 * sigma * sigma * fp * fp)).exp();
        let s_k = f.powi(-5) * (-1.25 * f_r.powi(-4)).exp() * 3.3f64.powf(peak);
        let a = (s_k * 2.0 * dk * dk).max(0.0).sqrt() * 4.0;
        let ph: f64 = rng.gen_range(0.0..2.0 * PI);
        for j in 0..ny {
            let y = j as f64 * lx / ny as f64;
            for i in 0..nx {
                let x = i as f64 * lx / nx as f64;
                data[j * nx + i] += C::from_polar(a, kx * x + ky * y + ph);
            }
        }
    }
    let rms: f64 = (data.iter().map(|c| c.norm_sqr()).sum::<f64>() / data.len() as f64).sqrt();
    let target = std::f64::consts::SQRT_2 * hs / 4.0;
    let scale = if rms > 0.0 { target / rms } else { 1.0 };
    for c in data.iter_mut() {
        *c = c.scale(scale);
    }
    Field2D { data, nx, ny, lx, ly: lx }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gaussian_mass_positive() {
        let f = gaussian_1d(512, 20.0, 1.0, 1.0, 10.0);
        assert!(f.mass() > 0.0);
    }

    #[test]
    fn soliton_shape_is_sech() {
        let lx = 40.0;
        let f = soliton_1d(4096, lx, 2.0, 0.0, lx / 2.0);
        // exact soliton of iψ_t + ψ_xx + |ψ|²ψ = 0: |ψ(0)| = √2·A.
        assert!((f.sup_norm() - std::f64::consts::SQRT_2 * 2.0).abs() < 1e-2);
        let xc = f.center_of_mass();
        assert!(((xc - lx / 2.0).rem_euclid(lx)).abs() < 0.1);
    }

    #[test]
    fn jonswap_envelope_has_expected_rms() {
        let f = crate::random::jonswap_envelope_1d(2048, 256.0, 4.0, 12.0, 42);
        let rms = (f.data.iter().map(|c| c.norm_sqr()).sum::<f64>() / f.data.len() as f64).sqrt();
        // For a narrow-band envelope, rms(|ψ|) ≈ √2 · Hs/4.
        let expected = std::f64::consts::SQRT_2 * 4.0 / 4.0;
        assert!(
            (rms - expected).abs() / expected < 0.2,
            "rms={rms} expected~{expected}"
        );
    }
}
