//! End-to-end validation of the split-step solver against exact solutions.
//!
//! 1. The focusing cubic NLS solver reproduces the **Akhmediev breather**
//!    evolution exactly (via the `ψ(x,t) = u(x, 2t)` scaling between the
//!    solver's ½-convention and the reference `u_xx + 2|u|²u` convention).
//! 2. The measured sideband growth of a Benjamin–Feir-unstable Stokes wave
//!    matches the analytic growth law `γ(ν) = ν√(4A² − ν²)`.

use rogue_nls::params::NlsParams;
use rogue_nls::solver::NlsSolver1D;
use rogue_rogue::benjamin_feir::bf_growth_rate;
use rogue_rogue::breathers::akhmediev;
use std::f64::consts::PI;

#[test]
fn solver_reproduces_akhmediev_breather() {
    let phi: f64 = 1.0;
    let period = PI / phi.sin();
    let lx = period;
    let nx = 1024;
    let params = NlsParams::focusing(lx, 2e-3, 3.0);
    let mut solver = NlsSolver1D::new(params, nx);
    let mut field = rogue_rogue::breathers::akhmediev_field(nx, lx, phi, 0.0);

    let n_steps = 500; // T = 1.0
    let mut err_max = 0.0f64;
    let s2 = std::f64::consts::SQRT_2;
    for step in 0..=n_steps {
        let t = step as f64 * params.dt;
        for i in 0..nx {
            let x = i as f64 * lx / nx as f64;
            let expected = s2 * akhmediev(x, t, phi, 0.0);
            let err = (field.data[i] - expected).norm();
            if err > err_max {
                err_max = err;
            }
        }
        if step < n_steps {
            solver.step(&mut field);
        }
    }
    assert!(
        err_max < 5e-3,
        "breather error too large: {err_max}"
    );
}

#[test]
fn benjamin_feir_growth_rate_measured() {
    let a0 = 1.0;
    let nu = 1.0;
    let lx = 16.0 * PI; // makes nu·lx/(2π) = 8 integer
    let nx = 2048;
    let eps = 1e-3;
    let params = NlsParams::focusing(lx, 1e-3, 3.0);
    let mut solver = NlsSolver1D::new(params, nx);

    // Stokes wave on a uniform background + sidebands at ±ν.
    let mut field = rogue_nls::scenario::plane_wave(nx, lx, 0.0, a0);
    for i in 0..nx {
        let x = i as f64 * lx / nx as f64;
        field.data[i] +=
            rogue_nls::field::C::new(eps * a0 * (nu * x).cos(), 0.0);
    }

    // Initial sideband amplitude in Fourier space (standard DIT ordering:
    // positive wavenumbers in the first half, negative in the second).
    let mut spec = field.data.clone();
    solver.engine().forward(&mut spec);
    let ks = &solver.engine().k;
    let j = ks
        .iter()
        .position(|&k| (k - nu).abs() < 1e-9)
        .expect("seed wavenumber on-grid");
    let a_init = spec[j].norm();
    assert!(a_init > 0.0, "sideband not seeded");

    // Sample ln(sideband) through the linear-growth stage; fit the slope
    // (the growth rate) over the late window where the cosh transient of
    // the non-exact-AB seed has settled onto the exponential e^{γt}.
    let fit_from = 4.0;
    let t_end = 6.0;
    let mut ts: Vec<f64> = Vec::new();
    let mut ls: Vec<f64> = Vec::new();
    let mut t = 0.0;
    while t < t_end {
        solver.step(&mut field);
        t += params.dt;
        if t >= fit_from && (t * 100.0).round() as i64 % 100 == 0 {
            let mut spec = field.data.clone();
            solver.engine().forward(&mut spec);
            ts.push(t);
            ls.push((spec[j].norm() / a_init).ln());
        }
    }
    assert!(ts.len() >= 4, "no fit samples");
    let n = ts.len() as f64;
    let (mt, ml) = (ts.iter().sum::<f64>() / n, ls.iter().sum::<f64>() / n);
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for (ti, li) in ts.iter().zip(&ls) {
        sxx += (ti - mt) * (ti - mt);
        sxy += (ti - mt) * (li - ml);
    }
    let gamma_measured = sxy / sxx;
    let gamma_theory = bf_growth_rate(nu, a0);
    assert!(
        (gamma_measured - gamma_theory).abs() / gamma_theory < 0.15,
        "measured {gamma_measured} vs theory {gamma_theory}"
    );
}
