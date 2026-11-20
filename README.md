# swimmers-ewald

Sumas hidrodinámicas en red para clusters de pares de spinners antiparalelos en flujo de Stokes.
Calcula los observables del artículo *Active motion and stability of a cluster of antiparallel spinners*
(Córdoba, Martinod, Schieber) con las convenciones aprobadas en la revisión y genera todas sus figuras de datos.

[![tests](https://github.com/th-martinod/swimmers-ewald/actions/workflows/tests.yml/badge.svg)](https://github.com/th-martinod/swimmers-ewald/actions/workflows/tests.yml)

- Convenciones: [`docs/CONVENTIONS.md`](docs/CONVENTIONS.md) (Figura 1 del artículo, ecs. 2 a 5).
- Correspondencia archivo ↔ figura del manuscrito: [`docs/FIGURE_MAP.md`](docs/FIGURE_MAP.md).
- Qué correr y qué revisar contra el texto del artículo: [`docs/CHECKLIST.md`](docs/CHECKLIST.md).

## Instalación y prueba rápida

Requiere Python ≥ 3.10.

```bash
git clone https://github.com/th-martinod/swimmers-ewald
cd swimmers-ewald
pip install -e ".[test]"
pytest -q                     # ~40 tests, unos segundos
```

Un cálculo suelto:

```python
from spinners.backends import make_field
from spinners.pair import observables

field = make_field("bulk3d", "hex", h=5.0)          # red hexagonal 3D, Q/h = 0.2
obs = observables(field, theta=0.2, Q=1.0, R=0.5)
print(obs["A_perp"], obs["B_perp"])                 # A_⊥/h y B_⊥/h en unidades de 1/Q
```

## Mapa del repositorio

```
src/spinners/
  lattices.py   bases de SCL, BCL, FCL y HEX (a⁽³⁾ sobre z), red recíproca, puntos de red ordenados por radio
  kernels.py    Stokeslet 3D y gradiente de su presión
  ewald.py      suma de Ewald (Hasimoto) para redes 3D, velocidades relativas a la media del fluido
  direct.py     sumas directas con corte radial y extrapolación L → ∞ (monocapas e imágenes de ∇p)
  backends.py   Field: geometrías bulk3d y monolayer, distancia a primeros vecinos d_nn
  pair.py       el par (ec. 1), puntos A, B, C, D y observables T, A_⊥, B_⊥
  sweep.py      barrido en paralelo definido por un .toml (comando spinners-sweep)
configs/paper.toml       el barrido de todas las figuras
scripts/paper_figures.py figuras del texto principal y del SI con los nombres del LaTeX
scripts/convergence.py   gráficas de convergencia
scripts/summary.py       números detrás de cada afirmación de la sección de resultados
scripts/run_all.sh       reproduce todo desde cero
results/                 se genera con scripts/run_all.sh (no se versiona en esta rama)
results/data/            CSV del barrido
results/figures/         figs/, figs_2/ (EPS y PDF), convergence/, preview/ (PNG)
results/convergence/     CSV de las pruebas de convergencia
tests/                   convenciones, simetrías, Ewald, promedio en el subvolumen
```

## Qué se calcula

Para cada red (SCL, BCL, FCL, HEX), geometría, radio R del spinner, parámetro de red h y ángulo θ:

| Columna | Significado |
|---|---|
| `T` | V_conf · e_⊥ con V_conf = [v(A) − v(B)] − [v(C) − v(D)] (ec. 4) |
| `A_perp` | A_⊥/h = −T sin θ (ec. 5). Negativo es estable |
| `B_perp` | B_⊥/h = V^f · e_⊥, con V^f el promedio de la velocidad del fluido en la bola de radio d_nn/2 centrada en el par (ec. 3) |
| `err_A`, `err_B` | distancia entre la suma parcial de mayor corte y el valor extrapolado (0 en Ewald para A) |

Unidades: μ = 1, ζV^S = 1, Q = 1.

Geometrías:

- `bulk3d`: red 3D de pares en un fluido 3D. Suma de Ewald con el desdoblamiento de Hasimoto. El corte gaussiano
  `cut = 36` deja errores relativos por debajo de 10⁻¹² (ver `conv_ewald`).
- `monolayer`: red 2D {c₁a⁽¹⁾ + c₂a⁽²⁾} en un fluido 3D sin fronteras, en reposo lejos. Las combinaciones del par
  son libres de fuerza, así que las sumas directas convergen y se extrapolan en el radio de corte con a + b/L + c/L².

El promedio en la bola se evalúa con el teorema del valor medio para campos biarmónicos (detalles en
`docs/CONVENTIONS.md`) y se valida contra una suma de Fourier independiente con el factor de forma de la bola.

## Tests

| Archivo | Garantiza |
|---|---|
| `test_symmetry.py` | sentido de giro de la Figura 1, v(r₁) = v(r₂), campo par respecto al centro, V_conf = 2[v(A) − v(B)], HEX estable para 0 < \|θ\| < π/6 e inestable en (π/6, π/4), periodicidad π/3, amplitud sin 6θ de la Sección S3, cancelación exacta en la SCL solo para R = Q/2, orientación de BCL y FCL |
| `test_ewald.py` | independencia respecto a ξ, Ewald contra suma directa 3D |
| `test_average.py` | promedio en la bola contra la suma de Fourier en las 8 combinaciones red × geometría, promedio cerrado del par aislado |

## Reproducir todo

```bash
bash scripts/run_all.sh        # unos 35 min con 2 núcleos, casi todo es el barrido
```

Para regenerar solo las figuras cuando ya existen los CSV de `results/data`:

```bash
python3 scripts/paper_figures.py
PYTHONPATH=src python3 scripts/convergence.py
```

## Extender

- Otra red: agregar su base en `lattices.basis` con a⁽³⁾ sobre z y la proyección de a⁽¹⁾ sobre +x.
- Otro observable: calcularlo en `pair._scalars` a partir de las velocidades en `eval_points` y agregarlo a
  `sweep.FIELDS`.
- Otro subvolumen: cambiar `self_ball_average` y el término de valor medio en `pair._scalars`, y el factor de forma
  en `tests/test_average.py`.

## Limitaciones

- Fuerzas puntuales. A_⊥ es un proxy cinemático del torque orientador, porque en una red de Bravais los centros
  de los dos spinners se advectan con la misma velocidad.
- El valor de B_⊥ depende del subvolumen de referencia. La bola de radio d_nn/2 es la mayor que no se solapa con
  las de los vecinos y no tiene orientación propia.
