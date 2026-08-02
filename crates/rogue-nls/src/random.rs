use crate::field::{Field1D, C};
use rand::Rng;
use rand_pcg::Pcg64Mcg;
use std::f64::consts::PI;

/// Seeded deterministic RNG (reproducible runs).
pub type SeededRng = Pcg64Mcg;

pub fn seeded(seed: u64) -> SeededRng {
    Pcg64Mcg::new(seed as u128)
}

/// Complex Gaussian white noise envelope, seeded.
pub fn gaussian_noise_1d(rng: &mut SeededRng, nx: usize, amp: f64) -> Vec<C> {
    let two = std::f64::consts::PI;
    (0..nx)
        .map(|_| {
            let u: f64 = rng.gen();
            let v: f64 = rng.gen();
            let mag = (-2.0 * u.ln()).sqrt();
            let ph = two * 2.0 * v;
            C::from_polar(amp * mag, ph)
        })
        .collect()
}

/// Deep-water JONSWAP directional-free spectrum in `k` (variance (Hs/4)²).
pub fn jonswap_spec_1d(k: &[f64], hs: f64, tp: f64, gamma: f64) -> Vec<f64> {
    let g = 9.81;
    let fp = 1.0 / tp;
    let mut s: Vec<f64> = k
        .iter()
        .map(|&kj| {
            if kj.abs() < 1e-12 {
                return 0.0;
            }
            let w = (g * kj.abs()).sqrt();
            let f = w / (2.0 * PI);
            let f_r = f / fp;
            let sigma = if f <= fp { 0.07 } else { 0.09 };
            let peak = (-((f - fp) * (f - fp)) / (2.0 * sigma * sigma * fp * fp)).exp();
            let s_f = f.powi(-5) * (-1.25 * f_r.powi(-4)).exp() * gamma.powf(peak);
            // S(f) df = S(k) dk ; deep-water k = w²/g ⇒ df/dk = g/(8π² f)
            let dfdk = g / (8.0 * PI * PI * f);
            s_f * dfdk
        })
        .collect();
    let var0: f64 = s.iter().map(|x| x * (k[1] - k[0]).abs()).sum();
    let target = (hs / 4.0) * (hs / 4.0);
    if var0 > 0.0 {
        let scale = target / var0;
        for x in s.iter_mut() {
            *x *= scale;
        }
    }
    s
}

/// Random-phase realization of the analytic envelope from a JONSWAP spectrum.
/// The real part is the free-surface elevation η(x); |ψ| is the wave envelope.
pub fn jonswap_envelope_1d(nx: usize, lx: f64, hs: f64, tp: f64, seed: u64) -> Field1D {
    let mut rng = seeded(seed);
    let dk = 2.0 * PI / lx;
    let k: Vec<f64> = (0..nx)
        .map(|j| (j as i64 - (nx / 2) as i64) as f64 * dk)
        .collect();
    let s = jonswap_spec_1d(&k, hs, tp, 3.3);
    let mut data = vec![C::new(0.0, 0.0); nx];
    for (kj, s_j) in k.iter().zip(s.iter()) {
        let ph: f64 = rng.gen_range(0.0..2.0 * PI);
        let a = (2.0 * s_j * dk).max(0.0).sqrt();
        if a == 0.0 {
            continue;
        }
        for i in 0..nx {
            let x = i as f64 * lx / nx as f64;
            data[i] += C::from_polar(a, kj * x + ph);
        }
    }
    Field1D { data, nx, lx }
}
