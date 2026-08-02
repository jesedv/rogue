#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

for s in regress-akhmediev.sh regress-ocean.sh regress-blowup.sh; do
  echo "== $s =="
  bash "scripts/$s"
done
bash scripts/regress-peregrine.sh
