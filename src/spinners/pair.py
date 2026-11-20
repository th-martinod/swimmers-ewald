"""One antiparallel spinner pair (Eq. 1 of the paper) and the observables of the article.

Conventions (docs/CONVENTIONS.md, Figure 1 of the article):
- n = z, Qhat = (cos th, sin th, 0), e = Qhat x n. The counterclockwise spinner is at r1 = -Q/2 and the
  clockwise spinner at r2 = +Q/2 (Q = r2 - r1, pair centre r_c = 0).
- Eq. (1) gives the force density F that the fluid exerts on the spinners: s_P zeta V^S Qhat at
  A = Q/2 + R e (s = +1), B = Q/2 - R e (-1), C = -Q/2 + R e (-1), D = -Q/2 - R e (+1).
  The spinners exert -F on the fluid (Eq. 2): FLUID_SIGN = -1.
- Eq. (4): V_conf = [v(A) - v(B)] - [v(C) - v(D)].
- A_perp = (mu h / zeta V^S) (Qhat x V_conf) . (a1hat x Qhat) = -(mu h / zeta V^S) T sin(theta),
  with T = V_conf . e. A_perp < 0 pushes Q back towards a1 (stabilizing).
- B_perp: [V^f]_perp = (zeta V^S / mu) (B_perp / h) e, with V^f the average of the fluid velocity over
  the ball of radius r_f = sqrt(2) h / 4 centred on the pair, the same for every lattice (SI Section S1).
  B_perp = B_self + B_img + r_f^2 gradp_img / 10, with B_self = -Q R / (10 pi r_f^3) the exact average of the
  flow of the pair itself, B_img the velocity of the images at the pair centre and gradp_img the pressure
  gradient of the images there (both along e). The last two are exported so that B_perp can be re-evaluated
  for any other radius without repeating the lattice sums.
Units: mu = 1, zeta V^S = 1, Q = 1, so A_perp/h and B_perp/h are given in units of 1/Q.
"""
import numpy as np

from .direct import extrapolate
from .kernels import stokeslet_dot

FLUID_SIGN = -1.0
NZ = np.array([0.0, 0.0, 1.0])


def frame(theta):
    qh = np.array([np.cos(theta), np.sin(theta), 0.0])
    return qh, np.cross(qh, NZ)


def sources(theta, Q=1.0, R=0.5):
    """Positions of the four point forces (A, B, C, D) and the forces exerted on the fluid."""
    qh, e = frame(theta)
    pos = np.array([Q / 2 * qh + R * e, Q / 2 * qh - R * e, -Q / 2 * qh + R * e, -Q / 2 * qh - R * e])
    s = np.array([1.0, -1.0, -1.0, 1.0])
    return pos, FLUID_SIGN * s[:, None] * qh[None, :]


def eval_points(theta, Q=1.0, R=0.5):
    """[A, B, C, D, r1, r2, r_c]."""
    P, _ = sources(theta, Q, R)
    qh, _ = frame(theta)
    return np.vstack([P, [-Q / 2 * qh, Q / 2 * qh, np.zeros(3)]])


def self_ball_average(P, F, a):
    """Exact average over the ball |r| < a of the free-space flow of forces F at P (all |P| < a), mu = 1:
    <G(r - P)> = [(2 a^2 - 4|P|^2/5) I + 2 P P / 5] / (8 pi a^3)."""
    P, F = np.asarray(P, float), np.asarray(F, float)
    p2 = np.sum(P * P, axis=1)
    pf = np.sum(P * F, axis=1)
    return ((2 * a**2 - 0.8 * p2)[:, None] * F + 0.4 * P * pf[:, None]).sum(0) / (8 * np.pi * a**3)


def _scalars(v, gp, theta, qh, e, P, F, a):
    """Observables from the velocities v[..., 7, 3] at eval_points and the image pressure gradient
    gp[..., 3] at the centre (any leading shape, e.g. one entry per cutoff radius)."""
    Vconf = (v[..., 0, :] - v[..., 1, :]) - (v[..., 2, :] - v[..., 3, :])
    T = Vconf @ e
    # average over the ball |r| < a: own pair exactly, images by the biharmonic mean-value theorem
    v_self0 = sum(stokeslet_dot(-P[s], F[s]) for s in range(4))
    v_img = (v[..., 6, :] - v_self0) @ e
    gp_img = gp @ e
    B = self_ball_average(P, F, a) @ e + v_img + a**2 / 10 * gp_img
    return {"T": T, "A_perp": -T * np.sin(theta), "B_perp": B, "B_img": v_img, "gradp_img": gp_img,
            "Qdot_perp": (v[..., 5, :] - v[..., 4, :]) @ e}


def observables(field, theta, Q=1.0, R=0.5):
    """field: a Field from backends.make_field. Returns a dict of scalar observables.

    For direct lattice sums every observable is extrapolated in the cutoff radius on its own, after
    the cancellations, which gives a realistic error estimate."""
    qh, e = frame(theta)
    P, F = sources(theta, Q, R)
    X = eval_points(theta, Q, R)
    a = field.r_f
    if field.direct:
        v = field.velocity_vs_L(X, P, F)
        gp = field.grad_p_vs_L(P, F)
        series = _scalars(v, gp, theta, qh, e, P, F, a)
        out, err = {}, {}
        for k, y in series.items():
            val, er = extrapolate(field.L, y[:, None])
            out[k], err[k] = float(val[0]), float(er[0])
    else:
        # Ewald velocities. The image sum of grad p decays as d^-5 (absolutely convergent, tail ~ L^-2 with
        # shell-to-shell fluctuations), so the largest cutoff is used as it is and the change from the previous
        # cutoff is the error estimate (no 1/L extrapolation, see scripts/convergence.py, conv_images)
        v = field.velocity(X, P, F)
        gps = field.grad_p_vs_L(P, F)
        series = _scalars(np.broadcast_to(v, (len(gps),) + v.shape), gps, theta, qh, e, P, F, a)
        out = {k: float(y[-1]) for k, y in series.items()}
        err = {k: float(abs(y[-1] - y[-2])) for k, y in series.items()}
    return {"T": out["T"], "A_perp": out["A_perp"], "B_perp": out["B_perp"], "B_img": out["B_img"],
            "gradp_img": out["gradp_img"], "Qdot_perp": out["Qdot_perp"],
            "err_A": err["A_perp"], "err_B": err["B_perp"]}
