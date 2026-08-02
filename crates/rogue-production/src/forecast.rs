//! The production forecast pipeline.
//!
//! Takes real sea-state observations (Hs, Tp), derives the carrier-scaling
//! NLS coefficients, seeds a JONSWAP envelope on a physical box, runs the
//! split-step solver, and reports rogue events + blow-up risk in measurable,
//! reproducible sea state (metres, seconds).

use crate::io::Observation;
use crate::sea::SeaState;
use rogue_blow_up::{BlowUpDetector, BlowUpSignals};
use rogue_nls::params::NlsParams;
use rogue_nls::solver::NlsSolver1D;
use rogue_nls::spectrum::Diagnostics;
use rogue_rogue::detect::RogueDetector;
use serde::{Deserialize, Serialize};

/// Box length, in carrier wavelengths.
const BOX_LAMBDAS: f64 = 12.0;
/// Grid points per peak wavelength.
const PPW: usize = 64;
/// Integration horizon, in carrier periods.
const FORECAST_PERIODS: f64 = 60.0;

/// One reported rogue event, in dimensional physical units.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ForecastEvent {
    /// Forecast time [s] from the start of the record.
    pub t: f64,
    /// Along-track location [m].
    pub x_m: f64,
    /// Crest elevation amplitude [m].
    pub amplitude_m: f64,
    /// Crest factor in units of sigma (rogue >= 2.2).
    pub crest_factor: f64,
    /// Significant wave height at the event [m].
    pub hs_m: f64,
}

/// Forecast output for a single observation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObservationForecast {
    pub obs_t: f64,
    pub hs: f64,
    pub tp: f64,
    pub wavelength_m: f64,
    pub max_crest_factor: f64,
    pub max_amplitude_m: f64,
    pub events: Vec<ForecastEvent>,
    pub blowup_risk: f64,
    pub blowup_eta_s: Option<f64>,
    pub note: String,
}

/// Entire report across the intake record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ForecastReport {
    pub version: &'static str,
    pub description: &'static str,
    pub forecasts: Vec<ObservationForecast>,
    pub summary: ForecastSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ForecastSummary {
    pub observations: usize,
    pub elevated_events: usize,
    pub max_crest_factor: f64,
    pub risk_level: RiskLevel,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RiskLevel {
    Low,
    Elevated,
    High,
    Severe,
}

impl RiskLevel {
    pub fn from_crest(cf: f64) -> Self {
        if cf >= 3.0 {
            RiskLevel::Severe
        } else if cf >= 2.6 {
            RiskLevel::High
        } else if cf >= 2.2 {
            RiskLevel::Elevated
        } else {
            RiskLevel::Low
        }
    }
}

/// Tuning for a forecast run.
#[derive(Debug, Clone, Copy)]
pub struct ForecastConfig {
    pub ppw: usize,
    pub box_lambdas: f64,
    pub forecast_periods: f64,
    pub rogue_threshold: f64,
}

impl Default for ForecastConfig {
    fn default() -> Self {
        Self {
            ppw: PPW,
            box_lambdas: BOX_LAMBDAS,
            forecast_periods: FORECAST_PERIODS,
            rogue_threshold: 2.2,
        }
    }
}

/// Run a physical forecast for a single observed sea state.
pub fn forecast_observation(
    obs: &Observation,
    cfg: &ForecastConfig,
    seed: u64,
) -> ObservationForecast {
    let sea = obs.sea_state();
    forecast_sea_state(obs.t, &sea, cfg, seed)
}

/// Core of the pipeline: run one physical sea state.
pub fn forecast_sea_state(
    obs_t: f64,
    sea: &SeaState,
    cfg: &ForecastConfig,
    seed: u64,
) -> ObservationForecast {
    let lam = sea.wavelength();
    let lx = cfg.box_lambdas * lam;
    let nx = (cfg.box_lambdas * cfg.ppw as f64).round() as usize;
    let per = sea.tp;
    let dt = per / 200.0;
    let t_end = cfg.forecast_periods * per;

    // Physical focusing cubic NLS scaling (beta from carrier, gamma negative).
    let params = NlsParams::custom(lx, dt, sea.beta(), -sea.gamma_c(), 3.0);
    let mut solver = NlsSolver1D::new(params, nx);
    // Seed a JONSWAP random sea at the observed Hs on the physical box.
    let mut field = rogue_nls::random::jonswap_envelope_1d(nx, lx, sea.hs, sea.tp, seed);

    let mut rogue = RogueDetector::new(nx * 8, cfg.rogue_threshold, nx / 4);
    let mut blow = BlowUpDetector::new();
    let steps = (t_end / dt) as usize;

    let mut max_cf = 0.0f64;
    for s in 0..steps {
        solver.step(&mut field);
        if s % 40 == 0 {
            let d = Diagnostics::measure_1d(solver.engine(), &field, &params);
            let sig = BlowUpSignals::from_diag_1d(
                solver.time(),
                &d,
                field.virial_about_com(),
                field.center_of_mass(),
                3.0,
                1,
            );
            blow.push(sig);

            let (mut amax, mut loc) = (0.0f64, 0.0f64);
            for i in 0..nx {
                let a = field.data[i].norm();
                if a > amax {
                    amax = a;
                    loc = i as f64 * lx / nx as f64;
                }
            }
            let rms = (field.data.iter().map(|c| c.norm_sqr()).sum::<f64>() / nx as f64).sqrt();
            let cf = if rms > 0.0 { amax / rms } else { 0.0 };
            if cf > max_cf {
                max_cf = cf;
            }
            rogue.observe(solver.time(), amax, loc, Some(rms));
        }
    }

    let events: Vec<ForecastEvent> = rogue
        .events()
        .iter()
        .map(|e| ForecastEvent {
            t: obs_t + e.t,
            x_m: e.location,
            amplitude_m: e.amplitude,
            crest_factor: e.crest_factor,
            hs_m: 4.0 * e.sigma,
        })
        .collect();

    let b = blow.state();
    ObservationForecast {
        obs_t,
        hs: sea.hs,
        tp: sea.tp,
        wavelength_m: lam,
        max_crest_factor: max_cf,
        max_amplitude_m: field.sup_norm(),
        events,
        blowup_risk: if b.active { b.growth_rate.max(0.0) } else { 0.0 },
        blowup_eta_s: b.eta,
        note: note(max_cf),
    }
}

fn note(cf: f64) -> String {
    match RiskLevel::from_crest(cf) {
        RiskLevel::Severe => "SEVERE rogue regime - crest >= 3 sigma".to_string(),
        RiskLevel::High => "high rogue risk (>= 2.6 sigma)".to_string(),
        RiskLevel::Elevated => "elevated rogue risk (>= 2.2 sigma)".to_string(),
        RiskLevel::Low => "benign sea".to_string(),
    }
}

/// Run forecasts over every observation in a record and summarize.
pub fn forecast_observations(
    obs: &[Observation],
    cfg: &ForecastConfig,
    seed: u64,
) -> ForecastReport {
    let forecasts: Vec<ObservationForecast> = obs
        .iter()
        .map(|o| forecast_observation(o, cfg, seed))
        .collect();
    let elevated = forecasts
        .iter()
        .filter(|f| f.max_crest_factor >= cfg.rogue_threshold)
        .count();
    let max_cf = forecasts
        .iter()
        .map(|f| f.max_crest_factor)
        .fold(0.0f64, f64::max);
    ForecastReport {
        version: env!("CARGO_PKG_VERSION"),
        description: "rogue finite-time blow-up & rogue-wave predictor (production)",
        summary: ForecastSummary {
            observations: forecasts.len(),
            elevated_events: elevated,
            max_crest_factor: max_cf,
            risk_level: RiskLevel::from_crest(max_cf),
        },
        forecasts,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runs_one_record() {
        use crate::io::Observation;
        let obs = Observation {
            t: 0.0,
            hs: 1.0,
            tp: 9.0,
            gamma: 3.3,
        };
        let cfg = ForecastConfig {
            box_lambdas: 6.0,
            ppw: 32,
            forecast_periods: 20.0,
            rogue_threshold: 2.2,
        };
        let f = forecast_observation(&obs, &cfg, 7);
        assert!(f.wavelength_m > 0.0);
        assert!(f.max_crest_factor >= 0.0);
    }

    #[test]
    fn summary_reflects_input() {
        use crate::io::Observation;
        let list = vec![
            Observation { t: 0.0, hs: 1.0, tp: 9.0, gamma: 3.3 },
            Observation { t: 3600.0, hs: 3.5, tp: 11.0, gamma: 3.3 },
        ];
        let rep = forecast_observations(&list, &ForecastConfig::default(), 5);
        assert_eq!(rep.summary.observations, 2);
    }
}