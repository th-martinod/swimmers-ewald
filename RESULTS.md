# Resultados (red hexagonal, R = Q/2 y Q/4, Q/h = 0.025…0.25)

Datos en `results/data/hex_*.csv`, figuras en `results/figures/`. Todo con el código corregido y validado (25 tests).

## 1. El criterio de advección da rotación nula, exactamente

Para toda red de Bravais el campo de velocidades total es **par** respecto al centro del par:

- la red es centrosimétrica, $\{r_n\}=\{-r_n\}$,
- el campo de un par aislado es par, $u(-x)=u(x)$: girar el par 180° equivale a cambiar el signo de sus cuatro fuerzas,
  y la rotación cambia el signo de los vectores en el plano.

Entonces $v(r_1)=v(-r_2)=v(r_2)$ y $\dot Q=v(r_2)-v(r_1)=0$ para todo θ, R, Q/h y geometría.
Lo mismo vale para cualquier corrección de tamaño finito evaluada en los centros (Faxén), porque $\nabla^{2m}v$ también es par.
Numéricamente $|\omega_{\dot Q}|\le 10^{-14}$ en todo el barrido (`tests/test_symmetry.py`).

**Consecuencia:** con partículas puntuales en flujo de Stokes, una red perfecta no ejerce ningún torque orientador sobre el par.

## 2. Qué mide realmente $V_{\rm conf}$

Por la misma paridad, $(\hat Q\times V_{\rm conf})\cdot\hat n=-2\,[v(p_1)-v(p_2)]\cdot\hat e_\perp$ (verificado a 1e-10).
Es la diferencia de la velocidad normal entre los dos puntos de fuerza de un spinner, es decir, un **estiramiento local a lo ancho
del spinner**. No es una velocidad de rotación de Q.

## 3. Signo corregido del criterio del artículo

En las tres geometrías, para los dos valores de R y para todos los Q/h:

- $A_\perp>0$ para $0<|\theta|<\pi/6$: θ = 0 es **inestable**,
- $A_\perp<0$ para $\pi/6<|\theta|<\pi/4$: los puntos fijos estables están en $\theta=\pm\pi/6$ (Q apunta entre dos vecinos).

Es el resultado opuesto al publicado. La forma angular es la misma en las tres geometrías ($\propto\sin6\theta$).
La amplitud escala como $(Q/h)^5$ en bulk 3D y monocapa, y como $(Q/h)^4$ en 2D puro. Es proporcional a $R^2$
(R = Q/4 da exactamente 1/4).

## 4. Velocidad de traslación

- El par se mueve a lo largo de $-\hat Q\times\hat n$, que coincide con la flecha de la Fig. 1.
- La velocidad la domina el **par aislado** (interacción entre sus dos spinners):
  $U_{\rm iso}=0.02847$ (3D, R = Q/2), 0.01817 (3D, R = Q/4), $1/(5\pi)=0.06366$ (2D, R = Q/2), 0.03745 (2D, R = Q/4), en unidades ζV^S/(μQ).
- La red solo la corrige: −1.8 % (bulk), −2.5 % (monocapa) y +0.8 % (2D) en Q/h = 0.25, y tiende a 0 al diluir.
- La dependencia con θ tiene período π/3, como exige la simetría C6. Su amplitud es 6–8 % en Q/h = 0.25, 0.06–0.2 % en
  Q/h = 0.1 y < 1e-5 en Q/h = 0.025. Nada que ver con las oscilaciones de período π/2 y ~60 % de la Fig. 2, que venían del factor A.

El artículo afirma lo contrario en tres puntos: velocidad que crece fuerte con la concentración, que depende de θ y que
"surge solo de las interacciones". Las dos primeras eran artefactos del factor A y del prefactor $h^3/\pi^3$. La tercera es una cuestión de redacción: el par aislado ya se traslada
(Pingali y Luijten, Yeo et al.), y la red solo lo corrige.

## 5. Control SCL

Con R = Q/2 la tasa de $V_{\rm conf}$ se anula a precisión de máquina en las tres geometrías (máx. 5e-16, 1e-14 y 7e-13).
Con R = Q/4 no se anula (máx. 1.7e-4, 1.7e-4 y 4.6e-4). $\dot Q\equiv 0$ también en la SCL.
Ver `results/data/scl_*.csv`.
