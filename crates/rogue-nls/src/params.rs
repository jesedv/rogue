use serde::{Deserialize, Serialize};

/// Nonlinearity sign convention: defocusing = +|ψ|^(p-1) repulsion,
/// focusing = −|ψ|^(p-1) attraction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Focus {
    Focusing,
    Defocusing,
}

/// Solver parameters.
///
/// The equation solved is `i ∂t ψ = -(β2/2) ∇²ψ + γ |ψ|^(p-1) ψ` with
/// `γ = -|γ|` for focusing and `γ = +|γ|` for defocusing. For the
/// Merle–Raphael–Rodnianski–Szeftel form `i ∂t ψ + Δψ − |ψ|^(p-1) ψ = 0`
/// (β2 = 2), `γ = ±1` (default).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct NlsParams {
    /// Physical box length (per dimension). Grid spacing `dx = L / n`.
    pub lx: f64,
    /// Temporal step.
    pub dt: f64,
    /// Dispersion coefficient `β2`.
    pub beta2: f64,
    /// Nonlinearity power `p` (e.g. cubic `p = 3`).
    pub power: f64,
    /// Focusing/defocusing sign.
    pub focus: Focus,
    /// Nonlinearity strength magnitude `|γ|` (default 1).
    pub nonlin: f64,
}

impl NlsParams {
    /// Defocusing NLS `i ∂t ψ + Δψ − |ψ|^(p-1) ψ = 0` (MRRS convention, β2 = 2).
    pub fn defocusing(lx: f64, dt: f64, power: f64) -> Self {
        Self {
            lx,
            dt,
            beta2: 2.0,
            power,
            focus: Focus::Defocusing,
            nonlin: 1.0,
        }
    }

    /// Focusing NLS `i ∂t ψ + Δψ + |ψ|^(p-1) ψ = 0`.
    pub fn focusing(lx: f64, dt: f64, power: f64) -> Self {
        Self {
            lx,
            dt,
            beta2: 2.0,
            power,
            focus: Focus::Focusing,
            nonlin: 1.0,
        }
    }

    /// Optical-fiber NLSE scaling: `i ∂t ψ + (β2/2) ∂x²ψ + γ |ψ|² ψ = 0`.
    pub fn fiber(beta2: f64, gamma: f64, lx: f64, dt: f64) -> Self {
        Self {
            lx,
            dt,
            beta2,
            power: 3.0,
            focus: if gamma < 0.0 { Focus::Focusing } else { Focus::Defocusing },
            nonlin: gamma.abs(),
        }
    }

    /// Explicit coefficients: `i ∂t ψ + (β2/2) ∇²ψ + γ |ψ|^(p-1) ψ = 0`.
    pub fn custom(lx: f64, dt: f64, beta2: f64, gamma: f64, power: f64) -> Self {
        Self {
            lx,
            dt,
            beta2,
            power,
            focus: if gamma < 0.0 { Focus::Focusing } else { Focus::Defocusing },
            nonlin: gamma.abs(),
        }
    }

    /// Nonlinearity coefficient `γ` with sign (=-|γ| focusing, =+|γ| defocusing).
    #[inline]
    pub fn gamma(&self) -> f64 {
        match self.focus {
            Focus::Focusing => -self.nonlin,
            Focus::Defocusing => self.nonlin,
        }
    }

    /// `L²`-critical power exponent for dimension `d`: `p_c = 1 + 4/d`.
    pub fn critical_power(&self, d: usize) -> f64 {
        1.0 + 4.0 / d as f64
    }

    /// True when the problem is `L²`-supercritical in dimension `d`,
    /// i.e. the regime in which finite-time blow-up is possible (focusing).
    pub fn is_supercritical(&self, d: usize) -> bool {
        self.power > self.critical_power(d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn critical_power_table() {
        let p = NlsParams::focusing(1.0, 0.01, 5.0);
        assert!((p.critical_power(1) - 5.0).abs() < 1e-12);
        assert!((p.critical_power(2) - 3.0).abs() < 1e-12);
        assert!((p.critical_power(3) - (1.0 + 4.0 / 3.0)).abs() < 1e-12);
    }

    #[test]
    fn supercriticality() {
        // 1D cubic is subcritical (p=3 < 5); 2D cubic is critical (p=3 == 3);
        // 2D quintic is supercritical (p=5 > 3).
        let cubic = NlsParams::focusing(1.0, 0.01, 3.0);
        assert!(!cubic.is_supercritical(1));
        assert!(!cubic.is_supercritical(2));
        let quintic = NlsParams::focusing(1.0, 0.01, 5.0);
        assert!(quintic.is_supercritical(2));
    }
}
