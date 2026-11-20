# Lista para correr y revisar en local

## 1. Correr

```bash
pip install -e ".[test]"
bash scripts/run_all.sh          # tests, barrido, figuras, convergencia y resumen
```

El barrido (`configs/paper.toml`, 240 tareas) es lo que más tarda. Se puede acelerar con
`python3 -m spinners.sweep configs/paper.toml --out results/data --workers N`.

## 2. Revisar las gráficas de convergencia (`results/figures/convergence/`)

| Archivo | Lo esperado |
|---|---|
| `conv_ewald` | error < 10⁻¹² en el corte de producción (cut = 36) para los tres ξ |
| `conv_xi` | desviaciones de 10⁻¹⁶ a 10⁻¹¹ (la SCL con R = Q/4 queda más arriba porque T es pequeño) |
| `conv_direct` | T converge rápido. B_⊥ converge como 1/L (el par es un cuadrupolo de fuerzas) y la extrapolación de producción queda entre 10⁻⁵ y 10⁻³ relativo |
| `conv_images` | 3D: suma parcial estable en ~10⁻⁵ de B_⊥ desde L ≈ 10h (producción usa 24h) |
| `conv_ball` | la suma de Fourier con factor de forma de la bola se acerca a la evaluación por valor medio al crecer n |
| `conv_production` | error relativo máximo de A_⊥ y B_⊥ con los parámetros del barrido. La SCL con R = Q/2 se omite en A_⊥ porque es cero exacto |

## 3. Contrastar el texto con `results/logs/summary.log`

Afirmaciones del texto principal que dependen de los números y que pueden cambiar con la nueva definición de B_⊥
(promedio en la bola) y la convención aprobada de A_⊥:

| # | Afirmación | Bloque de `summary.py` |
|---|---|---|
| 1 | con R = Q/4 las velocidades bajan por un factor ≈ 0.7 | R = Q/4 vs R = Q/2 |
| 2 | con R = Q/4 los torques bajan por un factor ≈ 0.3 | R = Q/4 vs R = Q/2 |
| 3 | en 3D: SCL la mayor velocidad, luego HEX de cerca, luego BCL, FCL bastante menor | B_perp/h averaged |
| 4 | las velocidades en monocapa son menores que en 3D | 2D / 3D ratio |
| 5 | en monocapa SCL ≈ HEX, y BCL ≈ FCL bastante menores | B_perp/h averaged |
| 6 | en monocapa SCL y HEX casi no dependen de θ | spread/\|mean\| |
| 7 | FCL estable solo para θ > 0 (3D y monocapa) | stable window |
| 8 | BCL y HEX estables a ambos lados, ancho ≈ 1 rad en 3D, HEX 2D ≈ 1 rad, BCL 2D bastante más estrecho | stable window |
| 9 | la magnitud de A_⊥ crece con Q/h | growth of the mean \|A_perp\| |
| 10 | solo BCL y HEX alcanzan su velocidad máxima dentro de la región estable (3D), HEX 2D también | velocity reachable |
| 11 | BCL 2D solo alcanza su velocidad media, unas 2.5 veces menor que HEX 2D | velocity reachable |
| 12 | Fig. 7: HEX (2D y 3D) tiene los torques estabilizadores más grandes, BCL bastante menores | minima |
| 13 | Conclusiones: entre las redes estables, HEX 3D da la mayor velocidad | velocity reachable |

Si algún número cambia, se corrige en rojo en el texto (y en las figuras del SI los textos alternativos).

## 4. Copiar las figuras al manuscrito

`results/figures/figs/*` → `figs/` y `results/figures/figs_2/*` → `figs_2/`. Los nombres coinciden con los del
LaTeX (ver `docs/FIGURE_MAP.md`). Las vistas previas están en `results/figures/preview/`.
