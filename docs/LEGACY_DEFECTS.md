# Defectos de `base_sums.rkt` (versión de producción)

Verificados con el port `src/spinners/legacy.py`. Reproducibles con
`PYTHONPATH=src python3 scripts/legacy_vs_corrected.py` y cubiertos por `tests/test_legacy.py`.
Parámetros de producción (`parametros_experimentos.csv`): Q = 1, α = 1, h = 4…40, na = 10…50.

1. **Q distinta en las partes de Fourier y real.** Para `hex`, `randdirf` gira Q +30° y la proyección usa $b^{(1)}$.
   Como con α = 1 la parte real vale ~$e^{-\pi h^2}$, el resultado es la parte de Fourier en θ+30°, que es
   exactamente $-A_\perp(\theta)$. Los mapas HEX publicados tienen el signo invertido.
   Para `bcl` y `fcl`, Q de Fourier y Q real difieren 29–35° y Q no es perpendicular a $\hat n=\hat z$.
2. **Factor de red A en la velocidad.** Solo deja los $k$ con algún índice nulo (3D). Produce anisotropía
   artificial (HEX 3D: 59 % de variación con θ y período π/2, cuando la simetría exige ~0). También se aplica en la suma real,
   donde no tiene sentido.
3. **Prefactor espurio $h^3/\pi^3$ (y falta un 4).** legacy/correcto $=h^3/(4\pi^3)$ exactamente. Con h = 4…40,
   varía 1000 veces a lo largo del eje Q/h.
4. **Truncamiento.** n/h baja de 2.5 a 1.25. Para Q/h ≤ 0.0625 el resultado es ruido de truncamiento
   (error 340 %–14 000 %, signo invertido).
5. **Prefactor de la parte real 4/π demasiado grande.** Oculto porque con α = 1 la parte real es despreciable.
6. **"2D"** = Fourier con $k_z=0$ (fuerzas lineales) más una suma real 3D en el plano. No es consistente (depende de α)
   y no es una monocapa en un fluido 3D.
7. Desplazamiento de las fuerzas en Fourier `beadr·(Q×n)` sin dividir por |Q|: solo es correcto si |Q| = 1.
8. Corridas 3D de velocidad de SCL, BCL y FCL con `rotplane` vacío: si llegó así, las fuerzas apuntan a lo largo de $\hat n$.
9. Posprocesamiento: `Last[...]` sin criterio de convergencia, pares Vper/VperAv emparejados por posición,
   caso h = 28 duplicado (exp 12 y 13).
