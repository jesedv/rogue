#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

cargo fetch
(cd education && npm install)
./dev wasm
