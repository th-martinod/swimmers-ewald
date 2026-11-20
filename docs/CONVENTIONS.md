# Convenciones (las aprobadas para el artículo)

Todas coinciden con la Figura 1 y la Figura 6 del artículo corregido.

| Elemento | Convención |
|---|---|
| Unidades | μ = 1, ζV^S = 1, \|Q\| = 1. Así `A_perp` es A_⊥/h y `B_perp` es B_⊥/h, ambos en unidades de 1/Q |
| Marco | n̂ = ẑ, Q̂ = (cos θ, sin θ, 0), e_⊥ = Q̂ × n̂ |
| Spinners | r₁ = −Q/2 gira antihorario visto desde la punta de n̂, r₂ = +Q/2 gira horario, centro r_c = 0 |
| Ec. (2) | μ∇²v = ∇p + ΣF, con F la fuerza que el fluido ejerce sobre los spinners. Sobre el fluido actúa −F (`FLUID_SIGN = -1`) |
| Puntos de fuerza (Fig. 1) | A = Q/2 + R e_⊥, B = Q/2 − R e_⊥, C = −Q/2 + R e_⊥, D = −Q/2 − R e_⊥, signos s_A = s_D = +1, s_B = s_C = −1 |
| Ec. (4) | V_conf = [v(A) − v(B)] − [v(C) − v(D)] = Σ s_P v(P) |
| Escalar T | T = V_conf · e_⊥ (columna `T`) |
| Ec. (5) | A_⊥ = (μh/ζV^S)(Q̂ × V_conf)·(â₁ × Q̂) = −(μh/ζV^S) T sin θ. Se usa sin θ, no sign(θ) |
| Estabilidad | A_⊥ < 0 (azul) empuja Q hacia a⁽¹⁾. En HEX: estable para 0 < \|θ\| < π/6, porque T ≃ c₆ · 126 K(h) sin 6θ con c₆ > 0 |
| Ec. (3) | [V^f]_⊥ = (ζV^S/μ)(B_⊥/h) e_⊥, con V^f el promedio de la velocidad del fluido en el subvolumen de referencia |
| Subvolumen | bola de radio r_f = √2h/4 ≈ 0.354h centrada en el par, la misma en todas las redes y geometrías. Es la bola más grande que no se solapa con las de los primeros vecinos en ninguna de las cuatro redes (la mitad de la distancia entre vecinos de la FCL). Contiene las cuatro fuerzas del par y ninguna de los vecinos mientras Q/h ≤ 0.25 |
| Redes BCL y FCL | se rotan para que a⁽³⁾ quede sobre ẑ. θ se mide desde la proyección de a⁽¹⁾ sobre el plano ⊥ a⁽³⁾ |
| Geometrías | `bulk3d`: red 3D en fluido 3D (Ewald, velocidades relativas a la velocidad media del fluido). `monolayer`: red 2D {c₁a⁽¹⁾ + c₂a⁽²⁾} en un fluido 3D sin fronteras, en reposo lejos (sumas directas extrapoladas) |

## El promedio sobre el subvolumen

El promedio sobre la celda completa se anula para k ≠ 0 y el modo k = 0 queda indeterminado para nadadores
libres de fuerza, por eso V^f se define como promedio local. Se evalúa sin sumas de Fourier:

V^f = ⟨v_par⟩_bola + [v_per(0) − v_par(0)] + (a²/10μ) ∇p_img(0),  a = r_f

En el CSV: `B_perp = −QR/(10π r_f³) + B_img + r_f² gradp_img/10`. Las columnas `B_img` y `gradp_img` permiten
recalcular B_⊥ para cualquier otro radio sin repetir las sumas.

- ⟨v_par⟩_bola es el promedio exacto del flujo del propio par:
  (1/8πμa³) Σ_P [(2a² − 4|P|²/5) F_P + (2/5) P (P·F_P)] = −ζV^S Q R e_⊥/(10πμa³).
- El flujo de las imágenes es biarmónico dentro de la bola, así que su promedio es u(0) + a²Δu(0)/10
  (teorema del valor medio), con Δu = ∇p/μ.
- ∇p_img(0) suma las imágenes del par (decae como d⁻⁵, convergencia absoluta).

`tests/test_average.py` compara con una suma de Fourier independiente con factor de forma de la bola,
3(sin x − x cos x)/x³ con x = 2π|k|a. El acuerdo es mejor que 2·10⁻³ relativo en todas las redes y geometrías.

### Por qué un radio común

El término propio del par, −ζV^S QR/(10πμ r_f³), domina B_⊥. Con r_f = d_nn/2 (la primera versión), r_f³ cambiaba
2.8 veces entre la SCL y la FCL, y la comparación entre redes reflejaba sobre todo el tamaño de la bola. Con el mismo
radio en todas las redes, a igual Q/h el término propio es el mismo y las diferencias vienen solo de la red.
A_⊥ no depende de r_f.

## Advertencia física que ya está en el texto

Por simetría de inversión v(r₂) = v(r₁) en toda red de Bravais (`test_centres_are_advected_with_the_same_velocity`).
V_conf mide la diferencia de velocidad a través del diámetro 2R y [V_conf]_⊥ se usa como proxy cinemático del
torque orientador.
