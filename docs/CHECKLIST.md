# Running and checking the results

## 1. Run

```bash
pip install -e ".[test]"
bash scripts/run_all.sh          # tests, sweep, figures, convergence, checks and summary
```

`scripts/check_results.py` summarizes the checks of the data in `results/logs/checks.log`. All of them must PASS
(INFO lines are for information only).

The sweep (`configs/paper.toml`, 240 tasks) takes most of the time. It can be sped up with
`python3 -m spinners.sweep configs/paper.toml --out results/data --workers N`.

## 2. Check the convergence plots (`results/figures/convergence/`)

| File | Expected |
|---|---|
| `conv_ewald` | error < 10⁻¹² at the production cutoff (cut = 36) for the three values of ξ |
| `conv_xi` | deviations from 10⁻¹⁶ to 10⁻¹¹ (the SCL with R = Q/4 lies higher because T is small) |
| `conv_direct` | T converges quickly. B_⊥ converges as 1/L (the pair is a force quadrupole) and the relative error of the production extrapolation lies between 10⁻⁵ and 10⁻³ |
| `conv_images` | 3D: the partial sum is stable to ~10⁻⁵ of B_⊥ from L ≈ 10h (production uses 24h) |
| `conv_ball` | the Fourier sum with the form factor of the ball approaches the mean-value evaluation as n grows |
| `conv_production` | maximum relative error of A_⊥ and B_⊥ with the sweep parameters. The SCL with R = Q/2 is left out for A_⊥ because it is exactly zero |

## 3. Claims of the article and where to check them

R = Q/2 unless stated otherwise. The blocks are those of `results/logs/summary.log` (written by
`scripts/summary.py`), and the numbered checks are those of `results/logs/checks.log`.

| # | Claim | Where to check it |
|---|---|---|
| 1 | With R = Q/4 the velocities are smaller by a factor of about 0.5, since B_⊥ is proportional to R | block "R = Q/4 vs R = Q/2" |
| 2 | With R = Q/4 the torques are smaller by a factor of about 0.25, since A_⊥ is proportional to R² | block "R = Q/4 vs R = Q/2" |
| 3 | With R = Q/4 the velocity in the 3D HEX and BCL lattices oscillates with θ by up to 3 % and 6 % at Q/h = 0.25 | check 4 (HEX) and `results/data/bcl_bulk3d_R0.25.csv` (BCL) |
| 4 | In 3D the SCL gives the largest velocities, closely followed by HEX and then by BCL. The FCL velocities are about three times smaller than those of SCL and HEX | block "B_perp/h averaged over theta" |
| 5 | Monolayer velocities are about 5 % smaller than in 3D for SCL and HEX, and about 10 % and 60 % larger for BCL and FCL | block "2D / 3D ratio of the theta-averaged velocity" |
| 6 | In monolayers SCL and HEX give almost identical velocities, while BCL and FCL are smaller by about 15 % and 50 % | block "B_perp/h averaged over theta" |
| 7 | The velocity does not depend on θ in HEX (3D and monolayer) and in the 3D BCL, and oscillates in the other cases | spread/\|mean\| in block "B_perp/h averaged over theta", check 4 |
| 8 | In the SCL, A_⊥ vanishes identically for R = Q/2 (7 to 14 orders of magnitude below the other lattices) | check 3, `figs_2/magstability_*` (Fig. S3) |
| 9 | BCL and HEX in 3D, and HEX in monolayers, are stable for \|θ\| < π/6 at every Q/h. The monolayer BCL has the same window except at Q/h = 0.25, where a narrow unstable band appears next to θ = 0 on the clockwise side | block "stable window around theta = 0", check 4 |
| 10 | In the FCL, θ = 0 is not an equilibrium orientation. A_⊥ > 0 for 0 < θ < θ*, with θ* ≈ 0.09 rad in 3D and ≈ 0.08 rad in monolayers, and the pairs rotate towards θ* from any θ between −0.43 and 0.62 rad | block "stable window around theta = 0" (lower limit) and `results/data/fcl_*_R0.5.csv` (θ* and upper limit) |
| 11 | The magnitude of A_⊥ grows with Q/h | block "growth of the mean \|A_perp\|/h with Q/h" |
| 12 | Within their stable window only BCL and HEX reach their maximum velocity in 3D, and HEX also in monolayers | block "velocity reachable inside the stable window" |
| 13 | Within its stable window the monolayer BCL reaches between about 0.7 and 1 times the velocity of the monolayer HEX, and between 0.73 and 0.80 times at Q/h = 0.25 | block "velocity reachable inside the stable window" |
| 14 | Fig. 7: the stabilizing torques of HEX are almost the same in 3D and in monolayers, and those of the 3D BCL are about 1.4 to 1.6 times larger | block "minima of A_perp/h" |
| 15 | Among the stable arrangements, the 3D HEX lattice gives the largest translational velocity | block "velocity reachable inside the stable window" |

## 4. Figures

The files in `results/figures/figs/` and `results/figures/figs_2/` have the same names as the figures in the
LaTeX source of the article (see `docs/FIGURE_MAP.md`). Previews are in `results/figures/preview/`.
