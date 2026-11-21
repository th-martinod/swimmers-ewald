# Correspondence with the figures of the article

The files in `results/figures/figs/` and `results/figures/figs_2/` have the same names as the files included by
the LaTeX source of the article, so they can be copied directly onto its `figs/` and `figs_2/` folders.
Each one exists as `.eps` and `.pdf`. Previews are in `results/figures/preview/`.
Figures 1 and 6 (diagrams) are not produced by this code.

## Main text (R = Q/2)

| Figure | Files | Content |
|---|---|---|
| 2 | `figs/CM_av_all_3D`, `figs/CM_vs_ang_all_3D` | B̄_⊥/h vs Q/h and B_⊥/B̄_⊥ vs θ, 3D SCL, BCL, FCL and HEX lattices |
| 3 | `figs/CM_av_all_2D`, `figs/CM_vs_ang_all_2D` | the same for monolayers |
| 4 | `figs/stability_{BCL,FCL,HEX}_3D` | sign of A_⊥ in the (θ, Q/h) plane, blue A_⊥ < 0, red A_⊥ > 0, thick black line A_⊥ = 0, thin labeled lines \|A_⊥\|/h = 10^k |
| 5 | `figs/stability_{BCL,FCL,HEX}_2D` | the same for monolayers |
| 7 | `figs/conf_min_all` | minimum of A_⊥/h for θ > 0 (counterclockwise, solid line and filled marker) and θ < 0 (clockwise, dashed line and small open marker), HEX 2D, HEX 3D, BCL 3D. Where the two coincide, a filled marker with a white center is seen |

## Supplementary Information

| Figure | Files | Content |
|---|---|---|
| S1 | `figs_2/CM_av_all_3D_R_p25`, `figs_2/CM_vs_ang_all_3D_R_p25` | 3D velocities, BCL and HEX, R = Q/4 |
| S2 | `figs_2/CM_av_all_2D_R_p25`, `figs_2/CM_vs_ang_all_2D_R_p25` | monolayer velocities, HEX, R = Q/4 |
| S3 | `figs_2/magstability_2D`, `figs_2/magstability_3D` | average of \|A_⊥\|/h over θ ∈ [−π/4, π/4], R = Q/2. SCL values that are exactly zero are drawn at 10⁻²² on the logarithmic axis |
| S4 | `figs_2/stability_HEX_3D_R_p25`, `figs_2/stability_BCL_3D_R_p25` | 3D maps, R = Q/4 |
| S5 | `figs_2/stability_HEX_2D_R_p25` | monolayer HEX map, R = Q/4 |
| S6 | `figs_2/conf_min_all_R_p25` | minima of A_⊥/h, R = Q/4 |

## Construction details

- B̄_⊥ is the average over the uniform grid of 80 angles in [−π, π).
- The maps use the grid of 101 angles in [−π/4, π/4] (excluding θ = 0, where A_⊥ = 0 trivially) and the 15 values
  of Q/h in `configs/paper.toml` (0.025 to 0.25). They leave the top 8 % of the canvas free for the panel letter
  that `overpic` adds at (57, 98).
- The minima are refined with a parabola through the three grid points around the minimum. If there is no stable
  orientation on one side, the point is omitted.
- To draw the contour lines, the maps interpolate A_⊥/(Q/h)⁵ (which varies slowly) between the 15 computed rows.
  The labels are placed in the column θ = −0.26. The ripples next to the narrow red band of the FCL (0 < θ < 0.09)
  come from the resolution in θ (0.016 rad), not from the physics.
- Every series has a line and small markers. In the curves against θ (one per value of Q/h) the markers are placed
  every 4 points and shifted between lattices, so that overlapping curves can be told apart (BCL and HEX in 3D with
  R = Q/2).
- Markers: SCL square, BCL triangle, FCL circle, HEX diamond. The color palette remains distinguishable for readers
  with color vision deficiency.

## Convergence plots (`results/figures/convergence/`)

| File | What it shows |
|---|---|
| `conv_ewald` | relative error of T and B_⊥ vs the Gaussian cutoff of the Ewald sum, three values of ξ, 3D lattices |
| `conv_xi` | independence of T and B_⊥ from ξ (at the 10⁻¹² level) |
| `conv_direct` | partial sums of the monolayers vs the cutoff radius L, and error of the production extrapolation |
| `conv_images` | image sum of ∇p used in the sub-volume average |
| `conv_ball` | independent Fourier sum with the form factor of the ball vs the mean-value evaluation |
| `conv_production` | production parameters vs a reference with larger cutoffs, all lattices, values of h and R |
