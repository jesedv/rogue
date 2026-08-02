//! Hydrodynamic (Madelung) reformulation of NLS — the compressible-Euler
//! bridge.
//!
//! Write `ψ = √ρ e^{iφ}` with velocity `u = ∇φ`. Then the defocusing
//! cubic NLS `i ψ_t + ψ_xx − |ψ|²ψ = 0` (β2 = 2, γ = +1, 1D) is
//! equivalent to the isothermal compressible Euler system plus a
//! "quantum pressure" term:
//!
//! ```text
//! ρ_t + (ρu)_x = 0
//! u_t + u u_x + ρ_x = (1/2) ∂x( (√ρ)_xx / √ρ )
//! ```
//!
//! When the quantum-pressure term is negligible (slowly varying envelope),
//! this is exactly the isothermal Euler equations with sound speed `c_s = 1`.
//! This is the mathematical bridge used to couple NLS blow-up / rogue-wave
//! dynamics to fluid (compressible Euler) behavior.

use rogue_nls::field::Field1D;
use serde::{Deserialize, Serialize};

/// Hydrodynamic state reconstructed from the NLS wave function.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MadelungState {
    /// Density `ρ = |ψ|²`.
    pub rho: Vec<f64>,
    /// Velocity `u = ∂_x φ = Im(conj ψ · ψ_x)/ρ`.
    pub u: Vec<f64>,
    /// Quantum pressure `q = (1/2)(√ρ)_xx / √ρ`.
    pub quantum_pressure: Vec<f64>,
    pub nx: usize,
}

/// Central second difference on a periodic grid.
fn ddx(f: &[f64], dx: f64) -> Vec<f64> {
    let n = f.len();
    let mut out = vec![0.0; n];
    for i in 0..n {
        let ip = (i + 1) % n;
        let im = (i + n - 1) % n;
        out[i] = (f[ip] - f[im]) / (2.0 * dx);
    }
    out
}

fn d2dx(f: &[f64], dx: f64) -> Vec<f64> {
    let n = f.len();
    let mut out = vec![0.0; n];
    for i in 0..n {
        let ip = (i + 1) % n;
        let im = (i + n - 1) % n;
        out[i] = (f[ip] - 2.0 * f[i] + f[im]) / (dx * dx);
    }
    out
}

/// Reconstruct `(ρ, u)` and the quantum-pressure term from the 1D field.
/// Assumes the solver convention `β2 = 2` (unit dispersion), matching the
/// derivation above. For other `β2`, pass `beta2` to rescale the velocity
/// (`u = (β2/2) φ_x` scaling).
pub fn madelung(field: &Field1D, beta2: f64) -> MadelungState {
    let dx = field.dx();
    let n = field.nx;
    let mut rho = vec![0.0; n];
    for i in 0..n {
        rho[i] = field.data[i].norm_sqr();
    }
    // ψ_x via central difference
    let mut psix = vec![rogue_nls::field::C::new(0.0, 0.0); n];
    for i in 0..n {
        let ip = (i + 1) % n;
        let im = (i + n - 1) % n;
        psix[i] = (field.data[ip] - field.data[im]) / (2.0 * dx);
    }
    let mut u = vec![0.0; n];
    for i in 0..n {
        let j = field.data[i];
        let jx = psix[i];
        // Im( conj ψ · ψ_x ) / ρ
        let im = j.re * jx.im - j.im * jx.re;
        u[i] = if rho[i] > 1e-14 {
            im / rho[i] * beta2 / 2.0
        } else {
            0.0
        };
    }
    // quantum pressure: (β2/4)·(√ρ)_xx / √ρ  (β2=2 ⇒ (1/2)(√ρ)_xx/√ρ)
    let sqrt_rho: Vec<f64> = rho.iter().map(|r| r.sqrt()).collect();
    let d2 = d2dx(&sqrt_rho, dx);
    let quantum_pressure: Vec<f64> = (0..n)
        .map(|i| {
            if rho[i] > 1e-14 {
                beta2 * 0.25 * d2[i] / sqrt_rho[i]
            } else {
                0.0
            }
        })
        .collect();
    MadelungState {
        rho,
        u,
        quantum_pressure,
        nx: n,
    }
}

/// `∂x(ρu)` for the continuity equation (periodic, central).
pub fn div_rho_u(rho: &[f64], u: &[f64], dx: f64) -> Vec<f64> {
    let n = rho.len();
    let mut flux = vec![0.0; n];
    for i in 0..n {
        flux[i] = rho[i] * u[i];
    }
    ddx(&flux, dx)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_field_gives_rest_state() {
        let mut f = Field1D {
            data: vec![rogue_nls::field::C::new(1.0, 0.0); 64],
            nx: 64,
            lx: 20.0,
        };
        f.data[0] = rogue_nls::field::C::new(1.0, 0.0);
        let s = madelung(&f, 2.0);
        for r in &s.rho {
            assert!((r - 1.0).abs() < 1e-12);
        }
        for v in &s.u {
            assert!(v.abs() < 1e-12);
        }
    }

    #[test]
    fn plane_wave_gives_constant_velocity() {
        // ψ = exp(i k0 x) ⇒ u = (β2/2) k0, with k0 on the discrete grid.
        // The central-difference ψ_x picks up sinc(k0·dx), which we fold
        // into the expectation.
        let (nx, lx, beta2) = (128, 20.0, 2.0);
        let k0 = 2.0 * std::f64::consts::PI * 4.0 / lx;
        let f = rogue_nls::scenario::plane_wave(nx, lx, k0, 1.0);
        let s = madelung(&f, beta2);
        let dx = lx / nx as f64;
        let expected = beta2 * 0.5 * k0 * (k0 * dx).sin() / (k0 * dx);
        for v in &s.u {
            assert!((v - expected).abs() < 1e-9);
        }
    }

    #[test]
    fn quantum_pressure_vanishes_for_constant_field() {
        let mut f = Field1D {
            data: vec![rogue_nls::field::C::new(2.0, 0.0); 32],
            nx: 32,
            lx: 8.0,
        };
        f.data[0] = rogue_nls::field::C::new(2.0, 0.0);
        let s = madelung(&f, 2.0);
        for q in &s.quantum_pressure {
            assert!(q.abs() < 1e-12);
        }
    }
}
