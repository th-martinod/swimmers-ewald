# Correspondencia con las figuras del manuscrito

Los archivos de `results/figures/figs/` y `results/figures/figs_2/` tienen el mismo nombre que los que
incluye el LaTeX, así que basta copiarlos sobre las carpetas `figs/` y `figs_2/` del manuscrito.
Cada uno existe en `.eps` y `.pdf`. Hay vistas previas en `results/figures/preview/`.
Las figuras 1 y 6 (diagramas) no salen de este código.

## Texto principal (R = Q/2)

| Figura | Archivos | Contenido |
|---|---|---|
| 2 | `figs/CM_av_all_3D`, `figs/CM_vs_ang_all_3D` | B̄_⊥/h vs Q/h y B_⊥/B̄_⊥ vs θ, redes 3D SCL, BCL, FCL, HEX |
| 3 | `figs/CM_av_all_2D`, `figs/CM_vs_ang_all_2D` | igual para monocapas |
| 4 | `figs/stability_{BCL,FCL,HEX}_3D` | signo de A_⊥ en el plano (θ, Q/h), azul A_⊥ < 0, rojo A_⊥ > 0, línea negra A_⊥ = 0 |
| 5 | `figs/stability_{BCL,FCL,HEX}_2D` | igual para monocapas |
| 7 | `figs/conf_min_all` | mínimo de A_⊥/h para θ > 0 (antihorario, marcador lleno) y θ < 0 (horario, vacío), HEX 2D, HEX 3D, BCL 3D |

## SI

| Figura | Archivos | Contenido |
|---|---|---|
| S1 | `figs_2/CM_av_all_3D_R_p25`, `figs_2/CM_vs_ang_all_3D_R_p25` | velocidades 3D, BCL y HEX, R = Q/4 |
| S2 | `figs_2/CM_av_all_2D_R_p25`, `figs_2/CM_vs_ang_all_2D_R_p25` | velocidades monocapa HEX, R = Q/4 |
| S3 | `figs_2/magstability_2D`, `figs_2/magstability_3D` | promedio de \|A_⊥\|/h sobre θ ∈ [−π/4, π/4], R = Q/2. Los ceros exactos de la SCL se dibujan en 10⁻²² |
| S4 | `figs_2/stability_HEX_3D_R_p25`, `figs_2/stability_BCL_3D_R_p25` | mapas 3D, R = Q/4 |
| S5 | `figs_2/stability_HEX_2D_R_p25` | mapa monocapa HEX, R = Q/4 |
| S6 | `figs_2/conf_min_all_R_p25` | mínimos de A_⊥/h, R = Q/4 |

## Detalles de construcción

- B̄_⊥ es el promedio sobre la malla uniforme de 80 ángulos en [−π, π).
- Los mapas usan la malla de 101 ángulos en [−π/4, π/4] (sin θ = 0, donde A_⊥ = 0 trivialmente) y los 15 valores
  de Q/h de `configs/paper.toml` (0.025 a 0.25). Dejan libre el 8 % superior del lienzo para la letra que agrega
  `overpic` en (57, 98).
- Los mínimos se refinan con una parábola por los tres puntos de la malla alrededor del mínimo. Si no hay
  orientación estable de un lado, el punto se omite.
- Colores y marcadores siguen a las figuras originales (SCL cuadrado, BCL triángulo, FCL círculo, HEX rombo).
  El verde pasó a verde azulado y el café a ocre para que la paleta sea distinguible con daltonismo.

## Gráficas de convergencia (`results/figures/convergence/`)

| Archivo | Qué muestra |
|---|---|
| `conv_ewald` | error relativo de T y B_⊥ vs el corte gaussiano de Ewald, tres valores de ξ, redes 3D |
| `conv_xi` | independencia de T y B_⊥ respecto a ξ (nivel 10⁻¹²) |
| `conv_direct` | sumas parciales de las monocapas vs el radio de corte L y error de la extrapolación de producción |
| `conv_images` | suma de imágenes de ∇p usada en el promedio del subvolumen |
| `conv_ball` | suma de Fourier independiente con factor de forma de la bola vs la evaluación por valor medio |
| `conv_production` | parámetros de producción vs referencia con cortes mayores, todas las redes, h y R |
