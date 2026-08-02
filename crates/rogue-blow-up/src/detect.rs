use crate::signals::{least_squares, BlowUpSignals, SignalBuffer};
use serde::{Deserialize, Serialize};

/// Result of the blow-up detector at the current time.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BlowUpState {
    /// A blow-up precursor is active (energy-negative or accelerating H¹).
    pub active: bool,
    /// Rigorous guarantee: focusing with `E < 0`.
    pub guaranteed_by_energy: bool,
    /// Log-growth rate of `H¹`: `d ln‖∇ψ‖/dt` over the fitting window.
    pub growth_rate: f64,
    /// Estimated time to blow-up `T_* − t` (from `1/‖∇ψ‖²` linear fit).
    pub eta: Option<f64>,
    /// R² of the ETA fit (0..1).
    pub eta_r2: f64,
    /// Concentration signal `‖ψ‖_∞² / mass` (mass-conserved ⇒ grows on
    /// blow-up).
    pub concentration: f64,
    /// Virial curvature sign (−1 = collapsing).
    pub virial_sign: i8,
    /// Number of samples collected.
    pub sample_count: usize,
}

/// Streaming blow-up detector.
///
/// Feed it one [`BlowUpSignals`] snapshot per observation (every `report_dt`
/// of simulation time). It keeps a ring buffer for the ETA fit.
pub struct BlowUpDetector {
    buffer: SignalBuffer,
    window: usize,
    active_threshold: f64,
}

impl Default for BlowUpDetector {
    fn default() -> Self {
        Self::new()
    }
}

impl BlowUpDetector {
    pub fn new() -> Self {
        Self {
            buffer: SignalBuffer::new(128),
            window: 48,
            active_threshold: 0.15,
        }
    }

    /// Minimum log-growth rate of `H¹` (per unit time) that counts as
    /// blow-up acceleration.
    pub fn with_active_threshold(mut self, t: f64) -> Self {
        self.active_threshold = t;
        self
    }

    pub fn push(&mut self, s: BlowUpSignals) {
        self.buffer.push(s);
    }

    pub fn sample_count(&self) -> usize {
        self.buffer.len()
    }

    pub fn state(&self) -> BlowUpState {
        let samples = self.buffer.ordered();
        let n = samples.len();
        let last = samples.last();
        let Some(last) = last else {
            return BlowUpState {
                active: false,
                guaranteed_by_energy: false,
                growth_rate: 0.0,
                eta: None,
                eta_r2: 0.0,
                concentration: 0.0,
                virial_sign: 0,
                sample_count: 0,
            };
        };

        let concentration = if last.mass > 0.0 {
            last.sup * last.sup / last.mass
        } else {
            0.0
        };

        let win = samples.len().min(self.window);
        let slice = &samples[n - win..];

        // Growth rate of ln(H¹).
        let ts: Vec<f64> = slice.iter().map(|s| s.t).collect();
        let ln_h1: Vec<f64> = slice.iter().map(|s| s.h1.max(1e-12).ln()).collect();
        let (_, growth_rate, _) = least_squares(&ts, &ln_h1);

        // ETA fit: y = 1/‖∇ψ‖² ≈ a + b·t, blow-up at t* = −a/b.
        let mut eta = None;
        let mut eta_r2 = 0.0;
        if slice.len() >= 4 {
            let ys: Vec<f64> = slice.iter().map(|s| 1.0 / s.h1_sq.max(1e-30)).collect();
            let (a, b, r2) = least_squares(&ts, &ys);
            // y = a + b·t → 0 at t* = −a/b: valid only when a, b have
            // opposite signs (y decreasing to zero as t → t*).
            if a * b < 0.0 && r2 > 0.5 {
                let t_star = -a / b;
                let et = t_star - last.t;
                if et > 0.0 && et < 1e6 {
                    eta = Some(et);
                    eta_r2 = r2;
                }
            }
        }

        let virial_sign = crate::virial::virial_second_derivative(
            &ts,
            &samples[n - win..].iter().map(|s| s.virial).collect::<Vec<_>>(),
        )
        .map(|c| c.sign)
        .unwrap_or(0);

        let active = last.energy_negative || (growth_rate > self.active_threshold && eta.is_some());

        BlowUpState {
            active,
            guaranteed_by_energy: last.energy_negative,
            growth_rate,
            eta,
            eta_r2,
            concentration,
            virial_sign,
            sample_count: n,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synth(sup: f64, h1: f64, t: f64) -> BlowUpSignals {
        BlowUpSignals {
            t,
            mass: 1.0,
            energy: -0.1,
            h1,
            h1_sq: h1 * h1,
            sup,
            virial: 10.0,
            center: 0.0,
            supercritical: true,
            energy_negative: true,
        }
    }

    #[test]
    fn energy_negative_immediately_active() {
        let mut d = BlowUpDetector::new();
        d.push(synth(1.0, 1.0, 0.0));
        let s = d.state();
        assert!(s.guaranteed_by_energy);
        assert!(s.active);
    }

    #[test]
    fn eta_recovers_exact_blowup_time() {
        // h1² = 1/(T* − t) with T* = 1 ⇒ 1/h1² = 1 − t.
        let t_star = 1.0;
        let mut d = BlowUpDetector::new();
        for i in 0..32 {
            let t = 0.02 * i as f64;
            let h1 = (1.0 / (t_star - t)).sqrt();
            d.push(synth(2.0, h1, t));
        }
        let s = d.state();
        let eta = s.eta.expect("expected an ETA");
        assert!((eta - (t_star - 0.02 * 31.0)).abs() < 0.05, "eta={eta}");
        assert!(s.eta_r2 > 0.9);
    }
}
