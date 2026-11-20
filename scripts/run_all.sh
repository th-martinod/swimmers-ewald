#!/usr/bin/env bash
# Reproduce every result in results/ from scratch (~35 min with 2 cores, mostly the sweep).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/data results/logs
export PYTHONPATH=src
python3 -m pytest -q
python3 -m spinners.sweep configs/paper.toml --out results/data 2>&1 | tee results/logs/sweep.log
python3 scripts/paper_figures.py 2>&1 | tee results/logs/paper_figures.log
python3 scripts/convergence.py 2>&1 | tee results/logs/convergence.log
python3 scripts/summary.py 2>&1 | tee results/logs/summary.log
