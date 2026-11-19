# Correspondencia entre las figuras del artículo y las versiones corregidas

Todas salen de `results/data_paper/` (barrido `configs/paper.toml`) y de `results/convergence/`.
Se generan con:

```bash
spinners-sweep configs/paper.toml --out results/data_paper      # barrido completo
python3 scripts/paper_figures.py                                 # figuras de datos
PYTHONPATH=src python3 scripts/convergence.py                    # convergencias
```

Los diagramas (Fig. 1 y la figura de cancelación) no se regeneran: no salen del código.

| Artículo / SI | Archivo corregido (`results/figures/paper/`) | Qué cambia respecto a la publicada |
|---|---|---|
| Fig. 2 (velocidad, 3D) | `fig2_velocity_3D.png` | Velocidad del par $U_\perp$ en lugar de $B_\perp$ con factor A. Se añade el panel de corrección de red. La dependencia con θ es pequeña y tiene período compatible con la simetría de cada red |
| Fig. 3 (velocidad, 2D) | `fig3_velocity_monolayer.png` | "2D" = monocapa en un fluido 3D |
| Fig. 4 (estabilidad, 3D) | `fig4_stability_3D.png` | Signo corregido. BCL y FCL con la geometría del artículo ($\hat n=a^{(3)}/|a^{(3)}|$) |
| Fig. 5 (estabilidad, 2D) | `fig5_stability_monolayer.png` | Idem, monocapa |
| Fig. 7 (mínimos de $A_\perp$) | `fig7_minima.png` | Misma definición que el notebook (mín. sobre θ>0 y θ<0 en $[-\pi/4,\pi/4]$), para todas las redes |
| Fig. S1 (velocidad 3D, R = Q/4) | `figS1_velocity_3D_R025.png` | Las cuatro redes |
| Fig. S2 (velocidad 2D, R = Q/4) | `figS2_velocity_monolayer_R025.png` | Las cuatro redes |
| Fig. S3 ($\overline{|A_\perp|}$) | `figS3_mean_abs_A.png` | Añade 2D puro. SCL con R = Q/2 en round-off |
| Fig. S4 (estabilidad 3D, R = Q/4) | `figS4_stability_3D_R025.png` | Incluye SCL, que con R = Q/4 ya no es neutral |
| Fig. S5 (estabilidad 2D, R = Q/4) | `figS5_stability_monolayer_R025.png` | Las cuatro redes |
| Fig. S6 (mínimos, R = Q/4) | `figS6_minima_R025.png` | Todas las redes |
| — (extra) | `figS7_velocity_pure2d.png`, `figS8_stability_pure2d.png`, `figS9_stability_pure2d_R025.png`, `figS10_minima_pure2d.png` | Hidrodinámica 2D pura, que es lo que el código viejo calculaba como "2D" |
| Figs. comentadas del SI (convergencia de Ewald, $B_\parallel$, $B_\perp$) | `conv_ewald.png`, `conv_xi.png`, `figS11_converged_velocity_3D.png` | Error contra número de imágenes para tres ξ, independencia de ξ, velocidades convergidas |
| — (diagnóstico) | `conv_production.png` | Error de truncamiento con los parámetros de producción (α = 1, n de `parametros_experimentos.csv`) |
| — (sumas directas) | `conv_direct.png` | Convergencia en el radio de corte para monocapa y 2D puro |

## Definiciones usadas

- $A_\perp/h$ = `A_conf` = sign(θ)·$(\hat Q\times V_{\rm conf})\cdot\hat n$ (ec. 5), con la fuerza sobre el fluido de la ec. (2).
- Velocidad = $U_\perp$, la velocidad del fluido promediada en los centros de los dos spinners, proyectada sobre $\hat Q\times\hat n$.
  El artículo usaba $B_\perp$ (promedio volumétrico con el factor A), que con la integral exacta es idénticamente nulo.
- BCL y FCL: vectores primitivos de Hasimoto con lado de cubo $h$. La base se rota para que $a^{(3)}$ quede a lo largo de $z$,
  y θ se mide desde la proyección de $a^{(1)}$ en el plano perpendicular a $a^{(3)}$.
  Monocapa = red $\{c_1a^{(1)}+c_2a^{(2)}\}$ (inclinada respecto al plano de giro). 2D puro = líneas a lo largo de $\hat n$ por esos puntos.
