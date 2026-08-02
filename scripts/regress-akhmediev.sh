#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# Replay the Akhmediev breather and check the solver reproduces it exactly.
cargo test -p rogue-rogue --test exact_solutions solver_reproduces_akhmediev_breather
