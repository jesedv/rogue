#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# Replay the Peregrine soliton through the split-step solver.
cargo run --release -- peregrine 1.0 2>&1 | grep -E "crest|peak|mass|energy" | head
