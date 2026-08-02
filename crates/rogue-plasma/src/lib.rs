//! Langmuir-wave envelope dynamics via the NLS reduction.
//!
//! In the subsonic limit the Zakharov system reduces to a focusing cubic NLS
//! for the Langmuir envelope `ψ` (normalized electric field):
//!
//! ```text
//! i ψ_t + (3/2) ω_pe λ_D² ψ_xx + (ω_pe / (2 n₀)) |ψ|² ψ = 0
//! ```
//!
//! with `ω_pe` the plasma frequency, `λ_D` the Debye length, `n₀` the
//! background density. The modulational (Langmuir collapse) instability is
//! exactly the Benjamin–Feir instability of this NLS.

use rogue_nls::params::NlsParams;
use serde::{Deserialize, Serialize};

/// Plasma parameters.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PlasmaParams {
    pub omega_pe: f64,
    pub n0: f64,
    pub lambda_de: f64,
}

impl PlasmaParams {
    /// NLS parameters for the Langmuir envelope equation.
    pub fn to_nls(&self, lx: f64, dt: f64) -> NlsParams {
        let beta2 = 3.0 * self.omega_pe * self.lambda_de * self.lambda_de;
        let gamma = self.omega_pe / (2.0 * self.n0);
        NlsParams::custom(lx, dt, beta2, -gamma, 3.0)
    }

    /// Modulational-instability gain `γ(ν)` for a Langmuir wave of envelope
    /// amplitude `E0` at perturbation wavenumber `ν`.
    ///
    /// From the general NLS gain
    /// `γ(ν) = ν √(β2 g A² − (β2/2)² ν²)` (see `rogue-rogue::benjamin_feir`).
    pub fn mi_gain(&self, nu: f64, e0: f64) -> f64 {
        let beta2 = 3.0 * self.omega_pe * self.lambda_de * self.lambda_de;
        let g = self.omega_pe / (2.0 * self.n0);
        let inside = beta2 * g * e0 * e0 - 0.25 * beta2 * beta2 * nu * nu;
        if inside <= 0.0 {
            0.0
        } else {
            nu * inside.sqrt()
        }
    }

    /// Maximum modulational-instability growth rate and its wavenumber:
    /// `γ_max = g A²` at `ν = √(2g/β2) A`.
    pub fn max_mi_growth(&self, e0: f64) -> (f64, f64) {
        let beta2 = 3.0 * self.omega_pe * self.lambda_de * self.lambda_de;
        let g = self.omega_pe / (2.0 * self.n0);
        let nu_star = (2.0 * g / beta2).sqrt() * e0;
        let gamma_max = g * e0 * e0;
        (gamma_max, nu_star)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn params_reduce_to_focusing_cubic() {
        let p = PlasmaParams {
            omega_pe: 1.0,
            n0: 1e3,
            lambda_de: 0.05,
        };
        let nls = p.to_nls(20.0, 1e-3);
        assert_eq!(nls.power, 3.0);
        assert!(nls.gamma() < 0.0, "Langmuir envelope must be focusing");
    }

    #[test]
    fn max_growth_matches_gain_at_star() {
        let p = PlasmaParams {
            omega_pe: 2.0,
            n0: 500.0,
            lambda_de: 0.1,
        };
        let e0 = 1.0;
        let (gmax, nu_star) = p.max_mi_growth(e0);
        let g_at_star = p.mi_gain(nu_star, e0);
        assert!((gmax - g_at_star).abs() / gmax < 1e-9);
    }

    #[test]
    fn stable_beyond_band() {
        let p = PlasmaParams {
            omega_pe: 1.0,
            n0: 1e3,
            lambda_de: 0.05,
        };
        let e0 = 1.0;
        let (_gmax, nu_star) = p.max_mi_growth(e0);
        assert!(p.mi_gain(nu_star * 2.0, e0) == 0.0);
    }
}
