//! Rogue-wave physics: exact breather solutions of the focusing cubic NLS,
//! Benjamin–Feir instability growth rates, and crest-factor detection.

pub mod benjamin_feir;
pub mod breathers;
pub mod detect;

pub use benjamin_feir::{bf_growth_rate, analyze_spectrum, BfAnalysis};
pub use breathers::{
    breather, BreatherKind, BreatherType, akhmediev, km_breather, peregrine,
};
pub use detect::{RogueDetector, RogueEvent};
