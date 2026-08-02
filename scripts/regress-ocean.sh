#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# JONSWAP random sea + Benjamin–Feir sideband growth law verification.
cargo test -p rogue-nls scenario::tests::jonswap_envelope_has_expected_rms
cargo test -p rogue-rogue --test exact_solutions benjamin_feir_growth_rate_measured
cargo run --release --example probe_bf
