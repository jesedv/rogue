#!/usr/bin/env bash
set -euo pipefail

# Production forecast regression: parses a sea-state CSV, runs the physical
# bridge, and asserts the report summarizes as many observations as input.
cd "$(dirname "$0")/.."

CSV="${1:-data/sample-buoy.csv}"
N=$(wc -l < "$CSV" | tr -d ' ')
N=$((N - 1))   # header row

OUT=$(cargo run --release --quiet -- predict "$CSV" 2>&1)
echo "$OUT"

OBS=$(echo "$OUT" | grep -oP 'summary: \K[0-9]+' | head -1)
if [ "$OBS" != "$N" ]; then
  echo "FAIL: expected $N observations, got $OBS" >&2
  exit 1
fi

echo "PASS: production forecast summarized $N observation(s)"