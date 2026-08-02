//! Benjamin–Feir (modulational) instability of the focusing cubic NLS in
//! the solver convention `i ψ_t + ψ_xx + |ψ|²ψ = 0` (β2 = 2, focusing).
//!
//! A uniform background `A` is unstable to sidebands at perturbation
//! wavenumber `ν` when `0 < ν < √2·A`, with growth rate
//!
//! ```text
//! γ(ν) = ν · √(2A² − ν²)
//! ```
//!
//! maximal at `ν = A`: `γ_max = A²`.

use serde::{Deserialize, Serialize};

/// Sideband growth rate `γ(ν)` for background amplitude `A`.
pub fn bf_growth_rate(nu: f64, amp: f64) -> f64 {
    let d = 2.0 * amp * amp - nu * nu;
    if d <= 0.0 {
        0.0
    } else {
        nu * d.sqrt()
    }
}

/// Instability band: sideband is amplified iff `ν < √2·A`.
pub fn is_unstable(nu: f64, amp: f64) -> bool {
    nu > 0.0 && nu < std::f64::consts::SQRT_2 * amp
}

/// Wavenumber of maximum growth: `ν = A`.
pub fn bf_max_growth_wavenumber(amp: f64) -> f64 {
    amp
}

/// Maximum growth rate `A²`.
pub fn bf_max_growth(amp: f64) -> f64 {
    amp * amp
}

/// Growth timescale `1/γ_max`.
pub fn bf_timescale(amp: f64) -> f64 {
    1.0 / bf_max_growth(amp).max(1e-12)
}

/// Spectral analysis of a Stokes-wave field: identifies the carrier mode
/// and quantifies sideband growth potential.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BfAnalysis {
    pub carrier_wavenumber: f64,
    pub carrier_amplitude: f64,
    pub band_max_nu: f64,
    pub band_max_growth: f64,
    pub unstable_band_halfwidth: f64,
}

/// Analyze the spectrum `spec` (row-major, physical wavenumbers `k`) of a
/// quasi-monochromatic field: find the dominant carrier and report the
/// Benjamin–Feir instability band around it.
pub fn analyze_spectrum(k: &[f64], spec: &[num_complex::Complex64]) -> Option<BfAnalysis> {
    if spec.is_empty() {
        return None;
    }
    let mut best = (0usize, 0.0f64);
    for (i, c) in spec.iter().enumerate() {
        let a = c.norm();
        if a > best.1 {
            best = (i, a);
        }
    }
    if best.1 < 1e-12 {
        return None;
    }
    let (k0, a0) = (k[best.0], best.1);
    let nu_max = bf_max_growth_wavenumber(a0);
    let half = std::f64::consts::SQRT_2 * a0;
    Some(BfAnalysis {
        carrier_wavenumber: k0,
        carrier_amplitude: a0,
        band_max_nu: nu_max,
        band_max_growth: bf_max_growth(a0),
        unstable_band_halfwidth: half,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_growth_location() {
        let a = 1.0;
        let nu_star = bf_max_growth_wavenumber(a);
        assert!((nu_star - 1.0).abs() < 1e-9);
        let gmax = bf_max_growth(a);
        assert!((gmax - 1.0).abs() < 1e-9);
        assert!((bf_growth_rate(nu_star, a) - gmax).abs() < 1e-9);
    }

    #[test]
    fn band_edges_zero_growth() {
        assert!((bf_growth_rate(0.0, 1.0) - 0.0).abs() < 1e-12);
        let nu_max = std::f64::consts::SQRT_2;
        assert!((bf_growth_rate(nu_max, 1.0) - 0.0).abs() < 1e-12);
        assert!(bf_growth_rate(1.5, 1.0) == 0.0);
        assert!(!is_unstable(1.5, 1.0));
        assert!(is_unstable(1.0, 1.0));
    }
}
