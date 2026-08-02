#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# Blow-up detection: ETA recovers the exact 1/H² blow-up time, plus a live
# quintic run through the CLI.
cargo test -p rogue-blow-up
cargo run --release -- blowup 5.0 0.3 2>&1 | grep -E "eta|active|blow" | head
