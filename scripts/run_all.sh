#!/usr/bin/env bash
# Reproduce every result in results/ from scratch (~12 min with 2 cores).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/data results/figures results/logs
pytest -q
spinners-sweep configs/hex.toml         2>&1 | tee results/logs/sweep_hex.log
spinners-sweep configs/scl_control.toml 2>&1 | tee results/logs/sweep_scl.log
python3 scripts/make_figures.py hex
PYTHONPATH=src python3 scripts/legacy_vs_corrected.py | tee results/logs/legacy_vs_corrected.log
