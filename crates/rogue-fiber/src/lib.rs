//! Fiber-optic nonlinear Schrödinger equation (NLSE) utilities.
//!
//! The NLSE with group-velocity dispersion and Kerr nonlinearity:
//!
//! ```text
//! i ψ_z + (β2/2) ψ_tt + γ |ψ|² ψ = 0
//! ```
//!
//! Fundamental-soliton condition `A₀²T₀² = |β2|/γ`, spectral-width
//! diagnostics for supercontinuum / high-power damage prediction.

use rogue_nls::field::Field1D;
use serde::{Deserialize, Serialize};

/// Soliton parameters.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Soliton {
    pub amplitude: f64,
    pub width: f64,
}

/// Peak power for a fundamental soliton of width `T0`.
pub fn soliton_peak_power(beta2: f64, gamma: f64, t0: f64) -> f64 {
    beta2.abs() / (gamma.abs() * t0 * t0)
}

/// Fundamental soliton: `ψ = A0 sech(t/T0)` with `A0²T0² = |β2|/γ`.
pub fn fundamental_soliton(
    nx: usize,
    lx: f64,
    beta2: f64,
    gamma: f64,
    t0: f64,
    center: f64,
) -> (Soliton, Field1D) {
    let a0 = soliton_peak_power(beta2, gamma, t0).sqrt();
    let data: Vec<rogue_nls::field::C> = (0..nx)
        .map(|i| {
            let t = i as f64 * lx / nx as f64;
            let dt = (t - center).rem_euclid(lx);
            let s = a0 / (dt / t0).cosh();
            rogue_nls::field::C::new(s, 0.0)
        })
        .collect();
    let field = Field1D { data, nx, lx };
    (Soliton { amplitude: a0, width: t0 }, field)
}

/// Spectral diagnostics.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SpectralStats {
    /// RMS angular-frequency width `σ_ω`.
    pub sigma_omega: f64,
    /// Mean frequency `⟨ω⟩`.
    pub mean_omega: f64,
    /// −20 dB bandwidth ratio vs the input.
    pub bandwidth_ratio_20db: f64,
    /// Spectral broadening factor (σ_ω / σ_ω,0).
    pub broadening: f64,
}

/// RMS spectral width of the field (via FFT engine), in angular frequency.
pub fn spectral_stats(
    eng: &rogue_nls::fft::FftEngine1D,
    field: &Field1D,
    input_sigma: f64,
) -> SpectralStats {
    let mut spec = field.data.clone();
    eng.forward(&mut spec);
    // physical angular frequencies ω = k
    let total: f64 = spec.iter().map(|c| c.norm_sqr()).sum();
    let mean: f64 = spec
        .iter()
        .zip(&eng.k)
        .map(|(c, k)| k * c.norm_sqr())
        .sum::<f64>()
        / total.max(1e-30);
    let sigma = (spec
        .iter()
        .zip(&eng.k)
        .map(|(c, k)| (k - mean) * (k - mean) * c.norm_sqr())
        .sum::<f64>()
        / total.max(1e-30))
    .sqrt();
    SpectralStats {
        sigma_omega: sigma,
        mean_omega: mean,
        bandwidth_ratio_20db: 0.0,
        broadening: if input_sigma > 0.0 { sigma / input_sigma } else { 1.0 },
    }
}

/// Estimate the −20 dB bandwidth ratio.
pub fn bandwidth_ratio_20db(eng: &rogue_nls::fft::FftEngine1D, field: &Field1D, input: &Field1D) -> f64 {
    fn bw(eng: &rogue_nls::fft::FftEngine1D, f: &Field1D) -> f64 {
        let mut spec = f.data.clone();
        eng.forward(&mut spec);
        let pmax = spec.iter().map(|c| c.norm_sqr()).fold(0.0, f64::max);
        let thresh = pmax * 1e-2;
        let mut lo = f64::MAX;
        let mut hi = f64::MIN;
        for (c, k) in spec.iter().zip(&eng.k) {
            if c.norm_sqr() >= thresh {
                lo = lo.min(*k);
                hi = hi.max(*k);
            }
        }
        if hi < lo {
            0.0
        } else {
            hi - lo
        }
    }
    let b_out = bw(eng, field);
    let b_in = bw(eng, input).max(1e-12);
    b_out / b_in
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fundamental_soliton_condition() {
        let (s, f) = fundamental_soliton(2048, 40.0, 1.0, 1.0, 1.0, 20.0);
        // A0² T0² = |β2|/γ
        assert!((s.amplitude.powi(2) * s.width.powi(2) - 1.0).abs() < 1e-9);
        assert!((f.sup_norm() - s.amplitude).abs() < 1e-2);
    }
}
