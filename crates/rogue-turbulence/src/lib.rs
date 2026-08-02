//! Clear-air turbulence (CAT) detection on NLS-modelled atmospheric
//! amplitude fields.
//!
//! The NLS serves as a universal amplitude equation for slowly-varying
//! wave trains; intermittent extreme "bursts" in the envelope (measured
//! through the crest factor and the excess kurtosis of `|ψ|`) are the
//! same mathematical signature as oceanic rogue waves and are linked to
//! CAT severity for aviation safety.

use rogue_rogue::detect::{RogueDetector, RogueEvent, ROGUE_THRESHOLD};
use serde::{Deserialize, Serialize};

/// Kurtosis-based intermittency signals.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct Intermittency {
    /// Excess kurtosis of `|ψ|` (0 = Gaussian, > 0 = intermittent).
    pub excess_kurtosis: f64,
    /// Crest factor `max|ψ| / rms|ψ|` over the field.
    pub crest_factor: f64,
    /// Severity classification based on kurtosis.
    pub severity: Severity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum Severity {
    #[default]
    Smooth,
    Light,
    Moderate,
    Severe,
}

impl Severity {
    pub fn from_kurtosis(k: f64) -> Self {
        if k < 0.5 {
            Severity::Smooth
        } else if k < 2.0 {
            Severity::Light
        } else if k < 6.0 {
            Severity::Moderate
        } else {
            Severity::Severe
        }
    }
}

/// Compute the excess kurtosis and crest factor of the envelope amplitudes.
pub fn intermittency(amp: &[f64]) -> Intermittency {
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
    let kurtosis = if var > 1e-30 { m4 / n / (var * var) } else { 3.0 };
    let rms = var.sqrt();
    let max = amp.iter().cloned().fold(0.0, f64::max);
    Intermittency {
        excess_kurtosis: kurtosis - 3.0,
        crest_factor: if rms > 0.0 { max / rms } else { 0.0 },
        severity: Severity::from_kurtosis(kurtosis - 3.0),
    }
}

/// Burst detector for CAT: reuses the rogue crest-factor detector, fed with
/// per-step field statistics from an NLS run.
pub struct BurstDetector {
    inner: RogueDetector,
}

impl Default for BurstDetector {
    fn default() -> Self {
        Self {
            inner: RogueDetector::new(4096, ROGUE_THRESHOLD, 128),
        }
    }
}

impl BurstDetector {
    pub fn new(window: usize, threshold: f64, refractory: usize) -> Self {
        Self {
            inner: RogueDetector::new(window, threshold, refractory),
        }
    }

    pub fn observe(&mut self, t: f64, amp_max: f64, loc: f64, sigma: Option<f64>) {
        self.inner.observe(t, amp_max, loc, sigma);
    }

    pub fn events(&self) -> &[RogueEvent] {
        self.inner.events()
    }

    pub fn event_count(&self) -> usize {
        self.inner.event_count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn gaussian_noise_has_zero_excess_kurtosis() {
        let mut rng = rogue_nls::random::seeded(7);
        let amps: Vec<f64> = (0..20000)
            .map(|_| {
                let u: f64 = rand::Rng::gen(&mut rng);
                let v: f64 = rand::Rng::gen(&mut rng);
                let mag = (-2.0 * u.ln()).sqrt();
                let ph = std::f64::consts::TAU * v;
                mag * ph.cos()
            })
            .collect();
        let i = intermittency(&amps);
        assert!(i.excess_kurtosis.abs() < 0.15, "{}", i.excess_kurtosis);
        assert_eq!(i.severity, Severity::Smooth);
    }

    #[test]
    fn sparse_spikes_give_severe_kurtosis() {
        let mut amps = vec![0.5; 1000];
        for a in amps.iter_mut().take(5) {
            *a = 5.0;
        }
        let i = intermittency(&amps);
        assert_eq!(i.severity, Severity::Severe);
        assert!(i.excess_kurtosis > 6.0);
    }
}
