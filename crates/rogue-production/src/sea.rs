//! Physical bridge: real-world sea-state parameters -> NLS coefficients.
//!
//! The ocean NLS is a scaled equation on the slowly-varying envelope of the
//! carrier wave. Given the peak period `Tp`, deep-water dispersion gives the
//! dominant wavenumber `k0` via `omega^2 = g*k0`, from which we derive the
//! group-velocity dispersion `beta` and nonlinearity `gamma` that the solver
//! consumes.
//!
//! Reference: D. Zakharov, *Stability of periodic waves of finite amplitude
//! on deep water* (1968); the narrow-band Dysthe/Trulsen expansion.

use serde::{Deserialize, Serialize};

const G: f64 = 9.81; // standard gravity [m/s^2]

/// A dimensional, real-world sea state, as measured by a buoy or forecast.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct SeaState {
    /// Significant wave height Hs [m].
    pub hs: f64,
    /// Peak period Tp [s].
    pub tp: f64,
    /// JONSWAP peak-enhancement gamma (default 3.3).
    pub gamma: f64,
}

impl SeaState {
    pub fn new(hs: f64, tp: f64) -> Self {
        Self { hs, tp, gamma: 3.3 }
    }

    /// Peak angular frequency omega0 = 2*pi/Tp [rad/s].
    pub fn omega0(&self) -> f64 {
        2.0 * std::f64::consts::PI / self.tp
    }

    /// Dominant carrier wavenumber k0 from deep-water dispersion `omega^2 = g*k` [1/m].
    pub fn k0(&self) -> f64 {
        self.omega0() * self.omega0() / G
    }

    /// Peak wavelength lambda0 = 2*pi/k0 [m].
    pub fn wavelength(&self) -> f64 {
        2.0 * std::f64::consts::PI / self.k0()
    }

    /// Group-velocity dispersion coefficient beta = omega0/(8 k0^2) [m^2/s].
    pub fn beta(&self) -> f64 {
        self.omega0() / (8.0 * self.k0() * self.k0())
    }

    /// Nonlinearity coefficient gamma = (omega0 k0^2)/2 [1/(m s)].
    pub fn gamma_c(&self) -> f64 {
        self.omega0() * self.k0() * self.k0() / 2.0
    }

    /// Benjamin-Feir gain bandwidth scale: the dimensionless steepness
    /// `2*A*k0` sets the growth-rate peak (see rogue-rogue bf module).
    pub fn steepness(&self) -> f64 {
        2.0 * (self.hs / 4.0) * self.k0()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deep_water_dispersion_consistent() {
        let s = SeaState::new(3.0, 10.0);
        let err = (s.omega0() * s.omega0() / G - s.k0()).abs();
        assert!(err < 1e-9);
    }

    #[test]
    fn finite_and_positive() {
        let s = SeaState::new(3.0, 10.0);
        assert!(s.beta().is_finite());
        assert!(s.gamma_c().is_finite());
        assert!(s.steepness() > 0.0);
        assert!(s.wavelength() > 0.0);
    }
}