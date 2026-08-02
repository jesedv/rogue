//! Cross-validation: when the quantum-pressure term is negligible, the
//! defocusing NLS density follows the isothermal compressible Euler
//! equations (`c_s = 1`).

use rogue_fluid::euler::EulerState;
use rogue_fluid::madelung::madelung;
use rogue_nls::params::NlsParams;
use rogue_nls::scenario::gaussian_1d;
use rogue_nls::solver::NlsSolver1D;

#[test]
fn nls_density_matches_isothermal_euler() {
    let lx = 100.0;
    let nx = 2048;
    let sigma = 4.0;
    let amp = 0.5;
    let dx = lx / nx as f64;
    let params = NlsParams::defocusing(lx, 2e-3, 3.0); // β2=2, γ=+1 (MRRS)
    let mut solver = NlsSolver1D::new(params, nx);
    let mut field = gaussian_1d(nx, lx, amp, sigma, lx / 2.0);

    let m0 = madelung(&field, params.beta2);
    let mut euler = EulerState::from_parts(m0.rho.clone(), m0.u.clone());

    let t_end = 0.5;
    let n_steps = (t_end / params.dt) as usize;
    for _ in 0..n_steps {
        solver.step(&mut field);
        let dt_e = euler.stable_dt(dx, 0.4);
        euler.step(dx, dt_e.min(params.dt));
    }

    let m1 = madelung(&field, params.beta2);
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..nx {
        let diff = m1.rho[i] - euler.rho[i];
        num += diff * diff;
        den += m1.rho[i] * m1.rho[i];
    }
    let rel = (num / den).sqrt();
    assert!(
        rel < 0.08,
        "NLS vs Euler density relative L2 error too large: {rel}"
    );
}
