# Convenciones

- $\hat n=\hat z$. $Q=(\cos\theta,\sin\theta,0)$ con $|Q|=1$. $\hat e_\perp=\hat Q\times\hat n$.
- El spinner antihorario está en $r_1=-Q/2$ y el horario en $r_2=+Q/2$ (artículo: $Q=r_2-r_1$, $r_c=0$).
- Fuerzas puntuales de la ec. (1): $s_P\,\hat Q$ en $P=\pm Q/2\pm R\,\hat e_\perp$ con $s=(+,-,-,+)$ para
  $(Q/2+R\hat e,\ Q/2-R\hat e,\ -Q/2+R\hat e,\ -Q/2-R\hat e)$.
- Con la ec. (2) escrita como $\mu\nabla^2v=\nabla p+F$, la fuerza **sobre el fluido** es $-F$. Esa es la que hace girar
  $r_1$ en sentido antihorario y la que usa `base_sums.rkt` (`FLUID_SIGN = -1` en `pair.py`).
  La Sección S3 del SI usa $+F$ y por eso tiene el signo contrario.
- $V_{\rm conf}$ (ec. 5): $[v(Q/2-R\hat e)-v(Q/2+R\hat e)]-[v(-Q/2-R\hat e)-v(-Q/2+R\hat e)]$.
- Tasa de rotación de un vector $W$: $(\hat Q\times W)\cdot\hat n$, positiva antihoraria.
- $A_\perp=\mathrm{sign}(\theta)\,(\hat Q\times W)\cdot\hat n$. Es lo que hace `base_sums.rkt` (normaliza $a^{(1)}\times Q$).
  $A_\perp<0$ empuja $Q$ hacia $a^{(1)}$. En $\theta=0$ vale 0 por construcción.
- `bulk3d`: velocidades relativas a la velocidad media del fluido en la celda (modo $k=0$ excluido).
  `monolayer` y `pure2d`: relativas al fluido en reposo lejos de la capa.
- Red hexagonal 3D: hexagonal simple, capas apiladas a lo largo de $a^{(3)}$ con separación $h$.
