//! Crest-factor rogue-wave detection.
//!
//! A rogue wave is conventionally a crest (or trough-to-crest height) that
//! exceeds **2.2× the significant wave height**. For a narrow-band sea the
//! significant wave height `Hs = 4σ`, so the threshold on the NLS envelope
//! `|ψ|` is `|ψ|/σ > 2.2` (σ = long-window RMS of the envelope).

use serde::{Deserialize, Serialize};

/// Default rogue crest-factor threshold (2.2).
pub const ROGUE_THRESHOLD: f64 = 2.2;

/// A detected rogue-wave event.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct RogueEvent {
    pub t: f64,
    pub location: f64,
    pub amplitude: f64,
    pub crest_factor: f64,
    pub sigma: f64,
    pub significant_height: f64,
}

/// Streaming crest-factor detector with hysteresis (refractory period so a
/// single extreme peak yields a single event).
pub struct RogueDetector {
    window: usize,
    threshold: f64,
    sigma: f64,
    refractory: usize,
    samples_since_peak: usize,
    peak_t: f64,
    peak_loc: f64,
    peak_amp: f64,
    peak_cf: f64,
    last_amp: f64,
    events: Vec<RogueEvent>,
}

impl Default for RogueDetector {
    fn default() -> Self {
        Self::new(4096, ROGUE_THRESHOLD, 256)
    }
}

impl RogueDetector {
    pub fn new(window: usize, threshold: f64, refractory: usize) -> Self {
        Self {
            window,
            threshold,
            sigma: 0.0,
            refractory,
            samples_since_peak: usize::MAX,
            peak_t: 0.0,
            peak_loc: 0.0,
            peak_amp: 0.0,
            peak_cf: 0.0,
            last_amp: 0.0,
            events: Vec::new(),
        }
    }

    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    pub fn events(&self) -> &[RogueEvent] {
        &self.events
    }

    pub fn event_count(&self) -> usize {
        self.events.len()
    }

    /// Feed a per-time-step field observation: the global envelope maximum
    /// amplitude and its location, plus the window RMS `sigma` (long-window
    /// running RMS of `|ψ|`). `sigma` should be supplied externally from the
    /// full field history; `None` falls back to an internal exponential RMS.
    pub fn observe(&mut self, t: f64, amp_max: f64, loc: f64, sigma: Option<f64>) {
        if let Some(s) = sigma {
            self.sigma = s;
        } else {
            // exponential running RMS (slow)
            let a = self.last_amp;
            self.sigma = (self.sigma * self.sigma * (self.window as f64 - 1.0)
                + a * a)
                / self.window as f64;
            self.sigma = self.sigma.sqrt();
        }
        self.last_amp = amp_max;
        let cf = if self.sigma > 0.0 {
            amp_max / self.sigma
        } else {
            0.0
        };

        if cf >= self.threshold {
            if cf > self.peak_cf {
                // new crest peak — track it
                self.peak_t = t;
                self.peak_loc = loc;
                self.peak_amp = amp_max;
                self.peak_cf = cf;
                self.samples_since_peak = 0;
            } else {
                // still above threshold but below the tracked peak: count
                // how long the crest has been receding before emitting.
                self.samples_since_peak = self.samples_since_peak.saturating_add(1);
                if self.samples_since_peak >= self.refractory {
                    self.emit();
                }
            }
        } else if self.peak_cf >= self.threshold {
            self.samples_since_peak = self.samples_since_peak.saturating_add(1);
            if self.samples_since_peak >= self.refractory {
                self.emit();
            }
        } else {
            self.peak_cf = 0.0;
        }
    }

    fn emit(&mut self) {
        if self.peak_cf < self.threshold {
            return;
        }
        self.events.push(RogueEvent {
            t: self.peak_t,
            location: self.peak_loc,
            amplitude: self.peak_amp,
            crest_factor: self.peak_cf,
            sigma: self.sigma,
            significant_height: 4.0 * self.sigma,
        });
        self.peak_cf = 0.0;
        self.peak_amp = 0.0;
        self.samples_since_peak = usize::MAX;
    }

    pub fn clear(&mut self) {
        self.events.clear();
        self.sigma = 0.0;
        self.samples_since_peak = usize::MAX;
        self.peak_cf = 0.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_rogue_spike_yields_single_event() {
        let mut d = RogueDetector::new(1024, 2.2, 64);
        for i in 0..300 {
            let amp = if i == 100 { 5.0 } else { 1.0 };
            d.observe(i as f64, amp, 0.0, Some(1.0));
        }
        assert_eq!(d.event_count(), 1);
        let e = &d.events[0];
        assert!((e.crest_factor - 5.0).abs() < 1e-9);
        assert!((e.t - 100.0).abs() < 1.0);
        assert!((e.significant_height - 4.0).abs() < 1e-9);
    }

    #[test]
    fn calm_sea_produces_no_events() {
        let mut d = RogueDetector::new(1024, 2.2, 64);
        for i in 0..500 {
            d.observe(i as f64, 1.0, 0.0, Some(1.0));
        }
        assert_eq!(d.event_count(), 0);
    }
}
