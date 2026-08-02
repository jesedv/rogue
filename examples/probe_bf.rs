//! Regression probe: measure Benjamin–Feir sideband growth in the actual
//! solver for several on-grid wavenumbers and compare against the analytic
//! law `γ(ν) = ν√(2A² − ν²)` (solver convention `iψ_t + ψ_xx + |ψ|²ψ = 0`).

use rogue_nls::params::NlsParams;
use rogue_nls::solver::NlsSolver1D;
use rogue_rogue::benjamin_feir::bf_growth_rate;

fn main() {
    let a0 = 1.0;
    let lx = 16.0 * std::f64::consts::PI;
    let nx = 4096;
    let params = NlsParams::focusing(lx, 1e-3, 3.0);
    let ks = rogue_nls::fft::FftEngine1D::new(nx, lx).k;

    for m in [5usize, 6, 7, 8, 9, 10, 11] {
        let nu = m as f64 * 2.0 * std::f64::consts::PI / lx;
        let j = ks.iter().position(|&k| (k - nu).abs() < 1e-9).unwrap();
        let mut solver = NlsSolver1D::new(params, nx);
        let mut field = rogue_nls::scenario::plane_wave(nx, lx, 0.0, a0);
        for i in 0..nx {
            let x = i as f64 * lx / nx as f64;
            field.data[i] += rogue_nls::field::C::from_polar(1e-3 * a0, nu * x);
        }
        let mut spec = field.data.clone();
        solver.engine().forward(&mut spec);
        let a0m = spec[j].norm();

        let mut history: Vec<(f64, f64)> = Vec::new();
        let mut t = 0.0;
        while t < 6.0 {
            solver.step(&mut field);
            t += params.dt;
            if t >= 4.0 && (t * 100.0).round() as i64 % 100 == 0 {
                let mut spec = field.data.clone();
                solver.engine().forward(&mut spec);
                history.push((t, (spec[j].norm() / a0m).ln()));
            }
        }
        let n = history.len();
        let (mt, my) = (
            history.iter().map(|x| x.0).sum::<f64>() / n as f64,
            history.iter().map(|x| x.1).sum::<f64>() / n as f64,
        );
        let (mut sxx, mut sxy) = (0.0, 0.0);
        for (t, y) in &history {
            sxx += (t - mt) * (t - mt);
            sxy += (t - mt) * (y - my);
        }
        let gamma = sxy / sxx;
        let theory = bf_growth_rate(nu, a0);
        println!(
            "nu={nu:.4} (m={m})  measured_gamma={gamma:.4}  theory={theory:.4}  rel_err={:.3}",
            (gamma - theory).abs() / theory
        );
    }
}
