"""Geometry of one antiparallel spinner pair (Eq. 1 of the paper) and the observables.

Conventions (see docs/CONVENTIONS.md):
- n = z. Q = Q (cos th, sin th, 0). e_perp = Qhat x n. The CCW spinner sits at r1 = -Q/2,
  the CW spinner at r2 = +Q/2 (paper: Q = r2 - r1, r_c = 0).
- Eq. (1) lists point forces s_P * Qhat at P = +-Q/2 +- R e_perp. With Eq. (2) written as
  mu lap v = grad p + F, the force ON THE FLUID is -F, which is the one that spins r1 CCW.
  FLUID_SIGN = -1 implements that (it is also what base_sums.rkt does).
- Rotation rate of a velocity W:  rate(W) = (Qhat x W) . n   (> 0 counterclockwise).
- Paper criterion: A_perp = sign(theta) * rate(W) (base_sums.rkt normalizes a1 x Q, hence sign);
  A_perp < 0 pushes Q back towards a1 (stabilizing).
"""
import numpy as np

from .kernels import stokeslet3d_dot, stokeslet2d_dot

FLUID_SIGN = -1.0
NZ = np.array([0.0, 0.0, 1.0])


def frame(theta, Q=1.0):
    qh = np.array([np.cos(theta), np.sin(theta), 0.0])
    e = np.cross(qh, NZ)
    return qh, e


def sources(theta, Q=1.0, R=0.5, fluid_sign=FLUID_SIGN):
    qh, e = frame(theta, Q)
    pos = np.array([Q / 2 * qh + R * e, Q / 2 * qh - R * e, -Q / 2 * qh + R * e, -Q / 2 * qh - R * e])
    s = np.array([1.0, -1.0, -1.0, 1.0])
    return pos, fluid_sign * s[:, None] * qh[None, :]


def eval_points(theta, Q=1.0, R=0.5):
    """[p1, p2, p3, p4, r1, r2] with V_conf = (v(p1)-v(p2)) - (v(p3)-v(p4))  (Eq. 5)."""
    qh, e = frame(theta, Q)
    return np.array([Q / 2 * qh - R * e, Q / 2 * qh + R * e, -Q / 2 * qh - R * e, -Q / 2 * qh + R * e,
                     -Q / 2 * qh, Q / 2 * qh])


def isolated_pair_velocity(theta, Q, R, X, kernel="3d"):
    """Free-space velocity at X from the pair itself (no images)."""
    P, F = sources(theta, Q, R)
    kern = stokeslet3d_dot if kernel == "3d" else stokeslet2d_dot
    return sum(kern(X - P[s], F[s]) for s in range(4))


def observables(field, theta, Q=1.0, R=0.5, kernel="3d"):
    """field(X, P, F) -> (v (E,3), err (E,3)). Returns a dict of scalar observables."""
    qh, e = frame(theta, Q)
    P, F = sources(theta, Q, R)
    X = eval_points(theta, Q, R)
    v, err = field(X, P, F)
    Vconf = (v[0] - v[1]) - (v[2] - v[3])
    Qdot = v[5] - v[4]
    U = 0.5 * (v[4] + v[5])
    v0 = isolated_pair_velocity(theta, Q, R, X[4:], kernel)
    U_iso = 0.5 * (v0[0] + v0[1])
    rate = lambda W: float(np.dot(np.cross(qh, W), NZ))
    sg = float(np.sign(round(theta, 14)))
    econf = float(np.linalg.norm(err[0] + err[1] + err[2] + err[3]))
    eq = float(np.linalg.norm(err[4] + err[5]))
    return {
        "rate_conf": rate(Vconf),           # (Qhat x V_conf).n   [velocity units, as in the paper]
        "omega_qdot": rate(Qdot) / Q,       # angular velocity of Q from Qdot = v(r2) - v(r1)
        "A_conf": sg * rate(Vconf),         # paper's A_perp/h (V_conf criterion), <0 = stabilizing
        "A_qdot": sg * rate(Qdot),          # same sign convention with Qdot instead of V_conf
        "Vconf_par": float(np.dot(Vconf, qh)),
        "U_perp": float(np.dot(U, e)),       # translation of the pair centre along Qhat x n
        "U_par": float(np.dot(U, qh)),
        "U_perp_isolated": float(np.dot(U_iso, e)),
        "err_conf": econf,
        "err_qdot": eq,
    }
