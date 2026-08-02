//! Isothermal compressible Euler solver (Rusanov finite volume, periodic).
//!
//! Used to cross-check the NLS hydrodynamic limit: with `c_s = 1` the
//! pressure is `p(ρ) = ρ` and the equations match the NLS Madelung system
//! when the quantum-pressure term is negligible.

use serde::{Deserialize, Serialize};

/// Conservative state `w = (ρ, ρu)`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EulerState {
    pub rho: Vec<f64>,
    pub u: Vec<f64>,
    pub n: usize,
}

impl EulerState {
    pub fn from_parts(rho: Vec<f64>, u: Vec<f64>) -> Self {
        assert_eq!(rho.len(), u.len());
        let n = rho.len();
        Self { rho, u, n }
    }

    fn cons(&self) -> (Vec<f64>, Vec<f64>) {
        let m: Vec<f64> = self.rho.iter().zip(&self.u).map(|(r, u)| r * u).collect();
        (self.rho.clone(), m)
    }

    fn flux(r: f64, m: f64) -> (f64, f64) {
        let u = if r > 1e-14 { m / r } else { 0.0 };
        let p = r;
        (m, m * u + p)
    }

    fn wave_speed(r: f64, m: f64) -> f64 {
        let u = if r > 1e-14 { m / r } else { 0.0 };
        u.abs() + 1.0 // c_s = 1
    }
}

/// Rusanov numerical flux between states `wL = (rL, mL)` and `wR`.
fn rusanov(wl: (f64, f64), wr: (f64, f64)) -> (f64, f64) {
    let (fl, f2l) = EulerState::flux(wl.0, wl.1);
    let (fr, f2r) = EulerState::flux(wr.0, wr.1);
    let s = EulerState::wave_speed(wl.0, wl.1).max(EulerState::wave_speed(wr.0, wr.1));
    (
        0.5 * (fl + fr) - 0.5 * s * (wr.0 - wl.0),
        0.5 * (f2l + f2r) - 0.5 * s * (wr.1 - wl.1),
    )
}

impl EulerState {
    /// Heun (RK2) step of the isothermal Euler equations on a periodic grid.
    pub fn step(&mut self, dx: f64, dt: f64) {
        let stage = |w: &(Vec<f64>, Vec<f64>), out: &mut (Vec<f64>, Vec<f64>)| {
            let n = self.n;
            for i in 0..n {
                let wl = (w.0[i], w.1[i]);
                let wr = (w.0[(i + 1) % n], w.1[(i + 1) % n]);
                let f = rusanov(wl, wr);
                let f_prev = rusanov(
                    (w.0[(i + n - 1) % n], w.1[(i + n - 1) % n]),
                    wl,
                );
                out.0[i] = -dt / dx * (f.0 - f_prev.0);
                out.1[i] = -dt / dx * (f.1 - f_prev.1);
            }
        };
        let (r0, m0) = self.cons();
        let (mut k1, mut k2) = (self.cons(), (vec![0.0; self.n], vec![0.0; self.n]));
        // k1
        {
            let w0 = (r0.clone(), m0.clone());
            stage(&w0, &mut k1);
        }
        // w0 + dt*k1
        for i in 0..self.n {
            k2.0[i] = r0[i] + dt * k1.0[i];
            k2.1[i] = m0[i] + dt * k1.1[i];
        }
        // k2 = f(w0 + dt k1)
        {
            let w1 = (k2.0.clone(), k2.1.clone());
            stage(&w1, &mut k2);
        }
        // w = w0 + dt/2 (k1 + k2)
        for i in 0..self.n {
            let r = r0[i] + 0.5 * dt * (k1.0[i] + k2.0[i]);
            let m = m0[i] + 0.5 * dt * (k1.1[i] + k2.1[i]);
            self.rho[i] = r;
            self.u[i] = if r > 1e-14 { m / r } else { 0.0 };
        }
    }

    /// CFL-limited stable step.
    pub fn stable_dt(&self, dx: f64, cfl: f64) -> f64 {
        let s = self
            .u
            .iter()
            .map(|u| u.abs() + 1.0)
            .fold(0.0, f64::max);
        cfl * dx / s.max(1e-12)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_state_stays_uniform() {
        let mut s = EulerState::from_parts(vec![2.0; 64], vec![0.3; 64]);
        let dx = 0.2;
        for _ in 0..100 {
            let dt = s.stable_dt(dx, 0.4);
            s.step(dx, dt);
        }
        for r in &s.rho {
            assert!((r - 2.0).abs() < 1e-9);
        }
        for u in &s.u {
            assert!((u - 0.3).abs() < 1e-9);
        }
    }

    #[test]
    fn mass_conserved() {
        let mut s = EulerState::from_parts(vec![1.0; 128], vec![0.0; 128]);
        s.rho[60] = 2.5;
        s.rho[68] = 2.5;
        let m0: f64 = s.rho.iter().sum();
        let dx = 0.1;
        for _ in 0..200 {
            let dt = s.stable_dt(dx, 0.4);
            s.step(dx, dt);
        }
        let m1: f64 = s.rho.iter().sum();
        assert!((m0 - m1).abs() / m0 < 1e-9);
    }
}
