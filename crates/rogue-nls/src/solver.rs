use crate::field::{Field1D, Field2D, C};
use crate::fft::{FftEngine1D, FftEngine2D};
use crate::params::NlsParams;

/// 1D NLS solver (Strang split-step, symplectic, second-order).
pub struct NlsSolver1D {
    pub params: NlsParams,
    eng: FftEngine1D,
    spec: Vec<C>,
    lin_mult: Option<(f64, Vec<C>)>,
    step_count: u64,
    t: f64,
}

impl NlsSolver1D {
    pub fn new(params: NlsParams, nx: usize) -> Self {
        let eng = FftEngine1D::new(nx, params.lx);
        let spec = vec![C::default(); nx];
        Self {
            params,
            eng,
            spec,
            lin_mult: None,
            step_count: 0,
            t: 0.0,
        }
    }

    pub fn step_count(&self) -> u64 {
        self.step_count
    }

    pub fn time(&self) -> f64 {
        self.t
    }

    pub fn engine(&self) -> &FftEngine1D {
        &self.eng
    }

    fn linear_multiplier(&mut self) -> &[C] {
        let dt = self.params.dt;
        let beta2 = self.params.beta2;
        if self.lin_mult.as_ref().map(|(d, _)| *d) != Some(dt) {
            let phase_scale = -0.5 * beta2 * dt * 0.5;
            let mult: Vec<C> = self
                .eng
                .ksq
                .iter()
                .map(|k2| {
                    let (s, co) = (phase_scale * k2).sin_cos();
                    C::new(co, s)
                })
                .collect();
            self.lin_mult = Some((dt, mult));
        }
        &self.lin_mult.as_ref().unwrap().1
    }

    /// Advance the field by one `dt` (Strang: linear·½, nonlinear·1, linear·½).
    pub fn step(&mut self, field: &mut Field1D) {
        let mult = self.linear_multiplier().to_vec();
        self.spec.copy_from_slice(&field.data);

        self.eng.forward(&mut self.spec);
        for (c, m) in self.spec.iter_mut().zip(&mult) {
            *c *= m;
        }
        self.eng.inverse(&mut self.spec);

        field.data.copy_from_slice(&self.spec);
        field.apply_nonlinear_phase(self.params.gamma(), self.params.power, self.params.dt);

        self.spec.copy_from_slice(&field.data);
        self.eng.forward(&mut self.spec);
        for (c, m) in self.spec.iter_mut().zip(&mult) {
            *c *= m;
        }
        self.eng.inverse(&mut self.spec);

        field.data.copy_from_slice(&self.spec);
        self.step_count += 1;
        self.t += self.params.dt;
    }

    pub fn step_n(&mut self, field: &mut Field1D, n: usize) {
        for _ in 0..n {
            self.step(field);
        }
    }
}

/// 2D NLS solver.
pub struct NlsSolver2D {
    pub params: NlsParams,
    eng: FftEngine2D,
    spec: Vec<C>,
    lin_mult: Option<(f64, Vec<C>)>,
    step_count: u64,
    t: f64,
}

impl NlsSolver2D {
    pub fn new(params: NlsParams, nx: usize, ny: usize) -> Self {
        let eng = FftEngine2D::new(nx, ny, params.lx, params.lx);
        let spec = vec![C::default(); nx * ny];
        Self {
            params,
            eng,
            spec,
            lin_mult: None,
            step_count: 0,
            t: 0.0,
        }
    }

    pub fn step_count(&self) -> u64 {
        self.step_count
    }

    pub fn time(&self) -> f64 {
        self.t
    }

    pub fn engine(&mut self) -> &mut FftEngine2D {
        &mut self.eng
    }

    fn linear_multiplier(&mut self) -> &[C] {
        let dt = self.params.dt;
        let beta2 = self.params.beta2;
        if self.lin_mult.as_ref().map(|(d, _)| *d) != Some(dt) {
            let phase_scale = -0.5 * beta2 * dt * 0.5;
            let mult: Vec<C> = self
                .eng
                .ksq
                .iter()
                .map(|k2| {
                    let (s, co) = (phase_scale * k2).sin_cos();
                    C::new(co, s)
                })
                .collect();
            self.lin_mult = Some((dt, mult));
        }
        &self.lin_mult.as_ref().unwrap().1
    }

    pub fn step(&mut self, field: &mut Field2D) {
        let mult = self.linear_multiplier().to_vec();
        self.spec.copy_from_slice(&field.data);

        self.eng.forward(&mut self.spec);
        for (c, m) in self.spec.iter_mut().zip(&mult) {
            *c *= m;
        }
        self.eng.inverse(&mut self.spec);

        field.data.copy_from_slice(&self.spec);
        field.apply_nonlinear_phase(self.params.gamma(), self.params.power, self.params.dt);

        self.spec.copy_from_slice(&field.data);
        self.eng.forward(&mut self.spec);
        for (c, m) in self.spec.iter_mut().zip(&mult) {
            *c *= m;
        }
        self.eng.inverse(&mut self.spec);

        field.data.copy_from_slice(&self.spec);
        self.step_count += 1;
        self.t += self.params.dt;
    }

    pub fn step_n(&mut self, field: &mut Field2D, n: usize) {
        for _ in 0..n {
            self.step(field);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scenario;

    #[test]
    fn linear_solve_is_exact_dispersion() {
        // With nonlinearity disabled, a single mode ψ = exp(i k0 x) must
        // reproduce the analytic phase shift exp(−i (β2/2) k0² dt).
        let lx = 16.0;
        let nx = 128;
        let k0 = 2.0 * std::f64::consts::PI * 3.0 / lx;
        let mut params = NlsParams::defocusing(lx, 0.1, 3.0);
        params.nonlin = 0.0;
        let mut solver = NlsSolver1D::new(params, nx);
        let mut field = scenario::plane_wave(nx, lx, k0, 1.0);
        let before = field.data.clone();
        solver.step(&mut field);
        for (i, (c, b)) in field.data.iter().zip(&before).enumerate() {
            let x = i as f64 * lx / nx as f64;
            let expected = *b * C::from_polar(1.0, -0.5 * params.beta2 * k0 * k0 * params.dt);
            assert!((c - expected).norm() < 1e-9, "i={i}: {} vs {}", c, expected);
        }
    }

    #[test]
    fn soliton_walks_with_constant_velocity() {
        // Focusing cubic NLS soliton propagates rigidly; speed = 2 v0.
        let (lx, nx) = (80.0, 2048);
        let params = NlsParams::focusing(lx, 1e-3, 3.0);
        let mut solver = NlsSolver1D::new(params, nx);
        let mut field = scenario::soliton_1d(nx, lx, 1.0, 0.0, lx / 2.0);
        let m0 = field.mass();
        let sup0 = field.sup_norm();
        for _ in 0..2000 {
            solver.step(&mut field);
        }
        assert!((field.mass() - m0).abs() / m0 < 1e-6);
        assert!((field.sup_norm() - sup0).abs() / sup0 < 1e-4);
        let t = solver.time();
        let xc = field.center_of_mass();
        let expected = lx / 2.0 + 0.0 * t;
        // periodicity: nearest-image distance
        let mut drift = (xc - expected).rem_euclid(lx);
        if drift > lx * 0.5 {
            drift = lx - drift;
        }
        assert!(drift < lx * 0.01, "soliton drifted too far: {drift}");
    }

    #[test]
    fn energy_conserved_in_linear_regime() {
        let (lx, nx) = (10.0, 256);
        let params = NlsParams::defocusing(lx, 1e-3, 3.0);
        let mut solver = NlsSolver1D::new(params, nx);
        let mut field = scenario::gaussian_1d(nx, lx, 1.0, 1.0, lx / 2.0);
        let e0 = crate::spectrum::Diagnostics::measure_1d(solver.engine(), &field, &params).energy;
        for _ in 0..500 {
            solver.step(&mut field);
        }
        let e1 = crate::spectrum::Diagnostics::measure_1d(solver.engine(), &field, &params).energy;
        assert!((e0 - e1).abs() / e0.abs().max(1e-12) < 1e-6);
    }
}
