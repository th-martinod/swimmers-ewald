#!/usr/bin/env bash
# Reproduce every result in results/ from scratch (~55 min with 2 cores, mostly configs/paper.toml).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/data results/data_paper results/figures/paper results/logs
pytest -q
spinners-sweep configs/paper.toml --out results/data_paper 2>&1 | tee results/logs/sweep_paper.log
spinners-sweep configs/hex.toml         2>&1 | tee results/logs/sweep_hex.log
spinners-sweep configs/scl_control.toml 2>&1 | tee results/logs/sweep_scl.log
python3 scripts/paper_figures.py
PYTHONPATH=src python3 scripts/convergence.py | tee results/logs/convergence.log
python3 scripts/make_figures.py hex
PYTHONPATH=src python3 scripts/legacy_vs_corrected.py | tee results/logs/legacy_vs_corrected.log
