//! rogue — finite-time blow-up & rogue-wave predictor CLI.
//!
//! Usage:
//!   rogue akhmediev [phi] [T]
//!   rogue peregrine [T]
//!   rogue ocean [hs] [T]
//!   rogue blowup [power] [T]
//!   rogue soliton [T]
//!   rogue fiber [T]
//!
//! Every run is seeded and reproducible; mass/momentum/energy conservation
//! is reported, plus blow-up state and rogue-wave events.

use rogue_blow_up::{BlowUpDetector, BlowUpSignals};
use rogue_nls::params::NlsParams;
use rogue_nls::scenario;
use rogue_nls::solver::NlsSolver1D;
use rogue_nls::spectrum::Diagnostics;
use rogue_rogue::breathers::{akhmediev_field, peregrine_field};
use rogue_rogue::detect::RogueDetector;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).map(|s| s.as_str()).unwrap_or("help");
    match cmd {
        "akhmediev" => {
            let phi = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1.0);
            let t_end = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1.0);
            run_akhmediev(phi, t_end);
        }
        "peregrine" => {
            let t_end = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1.0);
            run_peregrine(t_end);
        }
        "ocean" => {
            let _hs = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4.0);
            let t_end = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(60.0);
            run_ocean(t_end);
        }
        "blowup" => {
            let power = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(7.0);
            let t_end = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.6);
            run_blowup(power, t_end);
        }
        "soliton" => {
            let t_end = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4.0);
            run_soliton(t_end);
        }
        "fiber" => {
            let t_end = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(5.0);
            run_fiber(t_end);
        }
        "predict" => {
            let path = args.get(2).map(|s| s.to_string()).unwrap_or_default();
            let json = args.iter().any(|a| a == "--json");
            run_predict(&path, json);
        }
        _ => {
            println!(
                "usage: rogue <akhmediev|peregrine|ocean|blowup|soliton|fiber|predict> [params]"
            );
            println!("  predict <csv> [--json]   production sea-state forecast");
        }
    }
}

fn report(label: &str, solver: &NlsSolver1D, _field: &rogue_nls::field::Field1D, d: &Diagnostics) {
    println!(
        "{label:>12} t={:8.4}  mass={:9.6}  momentum={:10.6}  energy={:11.6}  H1={:10.4}  sup={:9.4}",
        solver.time(),
        d.mass,
        d.momentum,
        d.energy,
        d.h1,
        d.sup
    );
}

fn run_akhmediev(phi: f64, t_end: f64) {
    let period = std::f64::consts::PI / phi.sin();
    let nx = 2048;
    let dt = 2e-3;
    let params = NlsParams::focusing(period, dt, 3.0);
    let mut solver = NlsSolver1D::new(params, nx);
    let mut field = akhmediev_field(nx, period, phi, 0.0);
    let d0 = Diagnostics::measure_1d(solver.engine(), &field, &params);
    println!("== Akhmediev breather (phi={phi}, period={period:.4}) ==");
    report("t=0", &solver, &field, &d0);
    let mut max_amp = field.sup_norm();
    let steps = (t_end / dt) as usize;
    for s in 0..steps {
        solver.step(&mut field);
        if (s % 100) == 0 {
            let d = Diagnostics::measure_1d(solver.engine(), &field, &params);
            report("", &solver, &field, &d);
        }
        max_amp = max_amp.max(field.sup_norm());
    }
    println!("peak |psi| reached: {max_amp:.4}");
}

fn run_peregrine(t_end: f64) {
    let lx = 60.0;
    let nx = 4096;
    let dt = 2e-3;
    let params = NlsParams::focusing(lx, dt, 3.0);
    let mut solver = NlsSolver1D::new(params, nx);
    let mut field = peregrine_field(nx, lx, 0.0);
    let d0 = Diagnostics::measure_1d(solver.engine(), &field, &params);
    println!("== Peregrine soliton ==");
    report("t=0", &solver, &field, &d0);
    let steps = (t_end / dt) as usize;
    for s in 0..steps {
        solver.step(&mut field);
        if (s % 100) == 0 {
            let d = Diagnostics::measure_1d(solver.engine(), &field, &params);
            report("", &solver, &field, &d);
        }
    }
}

fn run_ocean(t_end: f64) {
    let lx = 512.0;
    let nx = 8192;
    let dt = 5e-3;
    let params = NlsParams::focusing(lx, dt, 3.0);
    let mut solver = NlsSolver1D::new(params, nx);
    let mut field = scenario::stokes_sidebands(nx, lx, 1.0, 0.0, 0.6, 0.05, 42);
    let mut rogue = RogueDetector::new(nx * 2, 2.2, 40);
    println!("== Ocean envelope (focusing cubic NLS, seeded) ==");
    let steps = (t_end / dt) as usize;
    for s in 0..steps {
        solver.step(&mut field);
        if s % 40 == 0 {
            let mut amax = 0.0f64;
            let mut loc = 0.0f64;
            for i in 0..nx {
                let a = field.data[i].norm();
                if a > amax {
                    amax = a;
                    loc = i as f64 * lx / nx as f64;
                }
            }
            let rms = (field.data.iter().map(|c| c.norm_sqr()).sum::<f64>() / nx as f64).sqrt();
            rogue.observe(solver.time(), amax, loc, Some(rms));
        }
    }
    println!("rogue events: {}", rogue.event_count());
    for e in rogue.events() {
        println!(
            "  t={:8.3}  x={:7.3}  amp={:6.3}  crest_factor={:5.2}  Hs={:5.2}",
            e.t, e.location, e.amplitude, e.crest_factor, e.significant_height
        );
    }
}

fn run_blowup(power: f64, t_end: f64) {
    let lx = 40.0;
    let nx = 8192;
    let dt = 2e-4;
    let params = NlsParams::focusing(lx, dt, power);
    let mut solver = NlsSolver1D::new(params, nx);
    // Mass ~1.5× critical so the collapse stays in the clean self-similar
    // stage (1/H¹² ≈ T*−t) long enough for the ETA fit to lock on.
    let mut field = scenario::gaussian_1d(nx, lx, 1.2, 3.0, lx / 2.0);
    let mut detector = BlowUpDetector::new();
    let d0 = Diagnostics::measure_1d(solver.engine(), &field, &params);
    println!("== Supercritical blow-up (p={power}, supercritical={}) ==", params.is_supercritical(1));
    report("t=0", &solver, &field, &d0);
    let steps = (t_end / dt) as usize;
    let start = Instant::now();
    let mut ever_active = false;
    let mut best_eta: Option<(f64, f64)> = None; // (eta, r2)
    for s in 0..steps {
        solver.step(&mut field);
        if s % 100 == 0 {
            let d = Diagnostics::measure_1d(solver.engine(), &field, &params);
            let sig = BlowUpSignals::from_diag_1d(
                solver.time(),
                &d,
                field.virial_about_com(),
                field.center_of_mass(),
                power,
                1,
            );
            detector.push(sig);
            let st = detector.state();
            ever_active |= st.active;
            if let Some(e) = st.eta {
                if best_eta.map_or(true, |(_, r)| st.eta_r2 > r) {
                    best_eta = Some((e, st.eta_r2));
                }
            }
            let eta = st.eta.map(|e| format!("{e:.4}")).unwrap_or_else(|| "-".into());
            report("", &solver, &field, &d);
            if s % 500 == 0 {
                println!(
                    "  blowup: active={} growth={:+.4} eta={} r2={:.3} conc={:.3}",
                    st.active, st.growth_rate, eta, st.eta_r2, st.concentration
                );
            }
        }
    }
    let st = detector.state();
    let elapsed = start.elapsed();
    println!(
        "final: active={} guaranteed={} growth={:.4} eta={:?} r2={:.3} ({} steps in {:.2}s)",
        st.active,
        st.guaranteed_by_energy,
        st.growth_rate,
        st.eta,
        st.eta_r2,
        steps,
        elapsed.as_secs_f64()
    );
    let alert = if best_eta.is_some() {
        format!(
            "best ETA {:.4} (r²={:.2})",
            best_eta.map(|(e, _)| e).unwrap_or(0.0),
            best_eta.map(|(_, r)| r).unwrap_or(0.0)
        )
    } else {
        "no ETA (resolve further to catch pre-collapse stage)".to_string()
    };
    println!("alert: ever_active={ever_active} · {alert}");
}

fn run_soliton(t_end: f64) {
    let lx = 60.0;
    let nx = 8192;
    let dt = 5e-4;
    let params = NlsParams::focusing(lx, dt, 3.0);
    let mut solver = NlsSolver1D::new(params, nx);
    let mut field = scenario::soliton_1d(nx, lx, 1.0, 0.5, lx / 2.0);
    let d0 = Diagnostics::measure_1d(solver.engine(), &field, &params);
    println!("== Focusing soliton ==");
    report("t=0", &solver, &field, &d0);
    let steps = (t_end / dt) as usize;
    for s in 0..steps {
        solver.step(&mut field);
        if s % 200 == 0 {
            let d = Diagnostics::measure_1d(solver.engine(), &field, &params);
            report("", &solver, &field, &d);
        }
    }
}

fn run_predict(path: &str, json: bool) {
    use rogue_production::ForecastConfig;
    let cfg = ForecastConfig::default();
    if path.is_empty() {
        println!("usage: rogue predict <sea-state.csv> [--json]");
        return;
    }
    match std::fs::read_to_string(path) {
        Ok(text) => match rogue_production::parse_observations(&text) {
            Ok(obs) => {
                let report = rogue_production::forecast_observations(&obs, &cfg, 42);
                if json {
                    match serde_json::to_string_pretty(&report) {
                        Ok(s) => println!("{s}"),
                        Err(e) => println!("json serialization error: {e}"),
                    }
                } else {
                    print_report(&report);
                }
            }
            Err(e) => println!("parse error: {e}"),
        },
        Err(e) => println!("cannot read '{path}': {e}"),
    }
}

fn print_report(report: &rogue_production::ForecastReport) {
    println!("rogue production forecast — {}", report.description);
    println!(
        "summary: {} observation(s), {} elevated, max crest factor {:.2} ({:?})",
        report.summary.observations,
        report.summary.elevated_events,
        report.summary.max_crest_factor,
        report.summary.risk_level
    );
    for f in &report.forecasts {
        println!(
            "  t={:8.1}s  Hs={:4.1}m  Tp={:4.1}s  λ={:6.1}m  maxCF={:5.2}  {note}",
            f.obs_t,
            f.hs,
            f.tp,
            f.wavelength_m,
            f.max_crest_factor,
            note = f.note
        );
        for e in f.events.iter().take(5) {
            println!(
                "      event t={:8.1}s  x={:7.1}m  crest={:5.2}σ  Hs={:4.1}m",
                e.t, e.x_m, e.crest_factor, e.hs_m
            );
        }
    }
}

fn run_fiber(t_end: f64) {
    let lx = 40.0;
    let nx = 8192;
    let dt = 5e-4;
    let beta2 = -1.0;
    let gamma = 1.0;
    let t0 = 1.0;
    let params = NlsParams::fiber(beta2, gamma, lx, dt);
    let (_sol, mut field) = rogue_fiber::fundamental_soliton(nx, lx, beta2, gamma, t0, lx / 2.0);
    let input_field = field.clone();
    let mut solver = NlsSolver1D::new(params, nx);
    let d0 = Diagnostics::measure_1d(solver.engine(), &field, &params);
    let stats0 = rogue_fiber::spectral_stats(solver.engine(), &field, 0.0);
    println!("== Fiber fundamental soliton (β2={beta2}, γ={gamma}) ==");
    report("t=0", &solver, &field, &d0);
    println!("  input σ_ω = {:.4}", stats0.sigma_omega);
    let steps = (t_end / dt) as usize;
    for s in 0..steps {
        solver.step(&mut field);
        if s % 200 == 0 {
            let d = Diagnostics::measure_1d(solver.engine(), &field, &params);
            report("", &solver, &field, &d);
        }
    }
    let stats = rogue_fiber::spectral_stats(solver.engine(), &field, stats0.sigma_omega);
    let b20 = rogue_fiber::bandwidth_ratio_20db(solver.engine(), &field, &input_field);
    println!(
        "  output σ_ω = {:.4}  broadening = {:.3}  −20 dB bandwidth ratio = {:.2}",
        stats.sigma_omega, stats.broadening, b20
    );
}
