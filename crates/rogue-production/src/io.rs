//! CSV / TSV intake for real-world sea-state observations.
//!
//! Supported input schemes:
//! - one row per observation with headers `t,hs,tp` (seconds, metres, seconds)
//! - optional `gamma` column; optional `gamma` default.
//!
//! Example header:
//! ```text
//! t,hs,tp
//! 0.0,2.5,9.4
//! 300.0,3.2,10.1
//! ```

use crate::sea::SeaState;
use serde::{Deserialize, Serialize};
use std::fmt;

/// A tide/observation from a buoy or forecast source.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct Observation {
    /// Observation time in seconds from the start of the record.
    pub t: f64,
    /// Significant wave height Hs [m].
    pub hs: f64,
    /// Peak period Tp [s].
    pub tp: f64,
    /// Optional JONSWAP gamma (default 3.3 when absent).
    pub gamma: f64,
}

#[derive(Debug)]
pub enum IntakeError {
    Io(String),
    IoRead(std::io::Error),
    MissingHeader,
    MissingColumn(String),
    BadField(usize, String, String),
}

impl fmt::Display for IntakeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IntakeError::Io(s) => write!(f, "intake: {s}"),
            IntakeError::IoRead(e) => write!(f, "read error: {e}"),
            IntakeError::MissingHeader => write!(f, "missing header"),
            IntakeError::MissingColumn(c) => write!(f, "missing required column '{c}'"),
            IntakeError::BadField(r, c, v) => write!(f, "bad {c} at row {r}: '{v}'"),
        }
    }
}

impl std::error::Error for IntakeError {}

impl From<std::io::Error> for IntakeError {
    fn from(e: std::io::Error) -> Self {
        IntakeError::IoRead(e)
    }
}

/// Parse a CSV/TSV buffer of sea-state observations.
pub fn parse_observations(text: &str) -> Result<Vec<Observation>, IntakeError> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header = lines
        .next()
        .ok_or(IntakeError::MissingHeader)?;
    let split: Vec<&str> = header.trim().split([',', '\t']).map(|s| s.trim()).collect();
    let mut idx_t = None;
    let mut idx_hs = None;
    let mut idx_tp = None;
    let mut idx_gamma = None;
    for (i, name) in split.iter().enumerate() {
        match name.to_ascii_lowercase().as_str() {
            "t" | "time" | "time_s" | "timestamp" | "t_s" => idx_t = Some(i),
            "hs" | "sigwaveheight" | "significant_wave_height" => idx_hs = Some(i),
            "tp" | "tpeak" | "peakperiod" | "peak_period" => idx_tp = Some(i),
            "gamma" | "enh" | "peak_sharp" => idx_gamma = Some(i),
            _ => {}
        }
    }
    let t_idx = idx_t.ok_or(IntakeError::MissingColumn("t".into()))?;
    let hs_idx = idx_hs.ok_or(IntakeError::MissingColumn("hs".into()))?;
    let tp_idx = idx_tp.ok_or(IntakeError::MissingColumn("tp".into()))?;

    let mut out = Vec::new();
    for (r, line) in lines.enumerate() {
        let fields: Vec<&str> = line.trim().split([',', '\t']).map(|s| s.trim()).collect();
        let num = |i: usize, name: &str| -> Result<f64, IntakeError> {
            fields
                .get(i)
                .ok_or_else(|| IntakeError::BadField(r, name.into(), "missing".into()))?
                .parse::<f64>()
                .map_err(|_| IntakeError::BadField(r, name.into(), fields[i].into()))
        };
        let t = num(t_idx, "t")?;
        let hs = num(hs_idx, "hs")?;
        let tp = num(tp_idx, "tp")?;
        if hs <= 0.0 || tp <= 0.0 {
            return Err(IntakeError::BadField(r, "hs/tp".into(), format!("{hs}/{tp}")));
        }
        let gamma = match idx_gamma {
            Some(g) => num(g, "gamma")?,
            None => 3.3,
        };
        let obs = Observation {
            t,
            hs,
            tp,
            gamma,
        };
        out.push(obs);
    }
    Ok(out)
}

impl Observation {
    pub fn sea_state(&self) -> SeaState {
        SeaState {
            hs: self.hs,
            tp: self.tp,
            gamma: self.gamma,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_comma_csv() {
        let csv = "t,hs,tp,gamma\n0,2.5,9.4,3.3\n300.0,3.2,10.1,3.0\n";
        let obs = parse_observations(csv).unwrap();
        assert_eq!(obs.len(), 2);
        assert!((obs[0].hs - 2.5).abs() < 1e-9);
        assert!((obs[1].tp - 10.1).abs() < 1e-9);
        assert!((obs[1].gamma - 3.0).abs() < 1e-9);
    }

    #[test]
    fn aliases_are_accepted() {
        let csv = "timestamp,Hs,Tpeak\n0,2.5,9.4\n0.1,2.6,9.5\n";
        let obs = parse_observations(csv).unwrap();
        assert_eq!(obs.len(), 2);
        assert_eq!(obs[0].sea_state().k0() > 0.0, true);
    }

    #[test]
    fn rejects_bad_numeric() {
        let csv = "t,hs,tp\n0,abc,9.4\n";
        assert!(parse_observations(csv).is_err());
    }

    #[test]
    fn rejects_negative_heights() {
        let csv = "t,hs,tp\n0,-2.5,9.4\n";
        assert!(parse_observations(csv).is_err());
    }
}