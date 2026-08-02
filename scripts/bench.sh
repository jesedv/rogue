#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# 2D NLS performance gate: 1024×1024 step must stay under 100 ms.
cargo run --release --bin bench_nls
