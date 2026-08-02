//! wasm-bindgen bindings for browser deployment of the rogue predictor.

use rogue_blow_up::{BlowUpDetector, BlowUpSignals};
use rogue_nls::params::NlsParams;
use rogue_nls::scenario;
use rogue_nls::solver::NlsSolver1D;
use rogue_rogue::detect::RogueDetector;
use serde::Serialize;
use wasm_bindgen::prelude::*;

#[wasm_bindgen(start)]
pub fn init() {
    console_error_panic_hook::set_once();
}

#[derive(Serialize)]
struct JsDiag {
    mass: f64,
    momentum: f64,
    energy: f64,
    h1: f64,
    sup: f64,
}

#[derive(Serialize)]
struct JsBlowUp {
    active: bool,
    guaranteed_by_energy: bool,
    growth_rate: f64,
    eta: Option<f64>,
    eta_r2: f64,
    concentration: f64,
    virial_sign: i8,
    supercritical: bool,
    sample_count: usize,
}

#[derive(Serialize)]
struct JsRogue {
    sigma: f64,
    crest_factor: f64,
    event_count: usize,
    events: Vec<JsRogueEvent>,
    kurtosis: f64,
}

#[derive(Serialize)]
struct JsRogueEvent {
    t: f64,
    location: f64,
    amplitude: f64,
    crest_factor: f64,
    significant_height: f64,
}

/// 1D NLS simulation wrapped for the browser dashboard.
#[wasm_bindgen]
pub struct NlsSim {
    solver: NlsSolver1D,
    field: rogue_nls::field::Field1D,
    blow: BlowUpDetector,
    rogue: RogueDetector,
    counter: usize,
    report_every: usize,
    nx: usize,
}

#[wasm_bindgen]
impl NlsSim {
    /// `scenario`: `akhmediev | peregrine | ocean | blowup | soliton | stokes`.
    #[wasm_bindgen(constructor)]
    pub fn new(nx: usize, lx: f64, dt: f64, scenario: &str, seed: u32) -> Result<NlsSim, JsValue> {
        let params = match scenario {
            "blowup" => NlsParams::focusing(lx, dt, 5.0),
            _ => NlsParams::focusing(lx, dt, 3.0),
        };
        let field = match scenario {
            "akhmediev" => {
                let phi: f64 = 1.0;
                let period = std::f64::consts::PI / phi.sin();
                rogue_rogue::breathers::akhmediev_field(nx, period, phi, 0.0)
            }
            "peregrine" => {
                rogue_rogue::breathers::peregrine_field(nx, lx, 0.0)
            }
            "ocean" => scenario::stokes_sidebands(nx, lx, 1.0, 0.0, 0.8, 0.05, seed as u64),
            "blowup" => scenario::gaussian_1d(nx, lx, 2.5, 1.2, lx / 2.0),
            "soliton" => scenario::soliton_1d(nx, lx, 1.0, 0.0, lx / 2.0),
            "stokes" => scenario::stokes_sidebands(nx, lx, 1.0, 0.0, 0.8, 0.05, seed as u64),
            other => return Err(JsValue::from_str(&format!("unknown scenario {other}"))),
        };
        let solver = NlsSolver1D::new(params, nx);
        Ok(NlsSim {
            solver,
            field,
            blow: BlowUpDetector::new(),
            rogue: RogueDetector::new(nx * 8, 2.2, nx / 4),
            counter: 0,
            report_every: 10,
            nx,
        })
    }

    pub fn step(&mut self) {
        self.solver.step(&mut self.field);
        self.counter += 1;
        if self.counter % self.report_every == 0 {
            let diag = rogue_nls::spectrum::Diagnostics::measure_1d(
                self.solver.engine(),
                &self.field,
                &self.solver.params,
            );
            let center = self.field.center_of_mass();
            let virial = self.field.virial_about_com();
            let s = BlowUpSignals::from_diag_1d(
                self.solver.time(),
                &diag,
                virial,
                center,
                self.solver.params.power,
                1,
            );
            self.blow.push(s);
            // rogue observation: global amplitude peak
            let mut amax = 0.0f64;
            let mut loc = 0.0f64;
            for i in 0..self.nx {
                let a = self.field.data[i].norm();
                if a > amax {
                    amax = a;
                    loc = i as f64 * self.field.lx / self.nx as f64;
                }
            }
            let rms = (self.field.data.iter().map(|c| c.norm_sqr()).sum::<f64>()
                / self.nx as f64)
                .sqrt();
            self.rogue.observe(self.solver.time(), amax, loc, Some(rms));
        }
    }

    pub fn step_n(&mut self, n: usize) {
        for _ in 0..n {
            self.step();
        }
    }

    pub fn time(&self) -> f64 {
        self.solver.time()
    }

    pub fn step_count(&self) -> u64 {
        self.solver.step_count()
    }

    pub fn amplitudes(&self) -> Vec<f32> {
        self.field.data.iter().map(|c| c.norm() as f32).collect()
    }

    pub fn surface(&self) -> Vec<f32> {
        self.field.data.iter().map(|c| c.re as f32).collect()
    }

    pub fn spectrum(&self) -> Vec<f32> {
        let mut spec = self.field.data.clone();
        self.solver.engine().forward(&mut spec);
        spec.iter().map(|c| c.norm() as f32).collect()
    }

    pub fn diagnostics(&self) -> JsValue {
        let d = rogue_nls::spectrum::Diagnostics::measure_1d(
            self.solver.engine(),
            &self.field,
            &self.solver.params,
        );
        serde_wasm_bindgen::to_value(&JsDiag {
            mass: d.mass,
            momentum: d.momentum,
            energy: d.energy,
            h1: d.h1,
            sup: d.sup,
        })
        .unwrap_or(JsValue::NULL)
    }

    pub fn blow_up_state(&self) -> JsValue {
        let s = self.blow.state();
        serde_wasm_bindgen::to_value(&JsBlowUp {
            active: s.active,
            guaranteed_by_energy: s.guaranteed_by_energy,
            growth_rate: s.growth_rate,
            eta: s.eta,
            eta_r2: s.eta_r2,
            concentration: s.concentration,
            virial_sign: s.virial_sign,
            supercritical: self.solver.params.power > 5.0,
            sample_count: s.sample_count,
        })
        .unwrap_or(JsValue::NULL)
    }

    pub fn rogue_stats(&self) -> JsValue {
        let events: Vec<JsRogueEvent> = self
            .rogue
            .events()
            .iter()
            .map(|e| JsRogueEvent {
                t: e.t,
                location: e.location,
                amplitude: e.amplitude,
                crest_factor: e.crest_factor,
                significant_height: e.significant_height,
            })
            .collect();
        let amps: Vec<f64> = self.field.data.iter().map(|c| c.norm()).collect();
        let kurt = rogue_turbulence_kurtosis(&amps);
        serde_wasm_bindgen::to_value(&JsRogue {
            sigma: self.rogue_est_sigma(),
            crest_factor: self.rogue_est_sigma().max(1e-12).recip() * self.field.sup_norm(),
            event_count: events.len(),
            events,
            kurtosis: kurt,
        })
        .unwrap_or(JsValue::NULL)
    }

    fn rogue_est_sigma(&self) -> f64 {
        (self.field.data.iter().map(|c| c.norm_sqr()).sum::<f64>() / self.nx as f64).sqrt()
    }
}

fn rogue_turbulence_kurtosis(amp: &[f64]) -> f64 {
    let n = amp.len() as f64;
    let mean = amp.iter().sum::<f64>() / n;
    let mut m2 = 0.0;
    let mut m4 = 0.0;
    for a in amp {
        let d = a - mean;
        m2 += d * d;
        m4 += d * d * d * d;
    }
    let var = m2 / n;
    if var > 1e-30 {
        m4 / n / (var * var) - 3.0
    } else {
        0.0
    }
}

#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}
