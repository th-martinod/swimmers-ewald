"""Triply periodic Stokeslet via the Hasimoto/Beenakker Ewald split.

G_per(x) = 1/(4 pi^2 V0) sum_{k!=0} phi(k) (I - kk/k^2)/k^2 cos(2 pi k.x)
          + c sum_n A_xi(x - r_n),     phi = (1 + pi xi k^2) exp(-pi xi k^2),  c = -1/(4 pi)

with A_xi the real-space kernel of base_sums.rkt (matA). The coefficient c was derived from
the short-distance limit (A_xi -> -(I/r + rr/r^3)/2) and is verified by xi-independence
(tests/test_ewald.py). Velocities are relative to the volume-mean fluid velocity (k = 0 excluded).
"""
import numpy as np
from scipy.special import erfc

from .lattices import reciprocal, points

PI = np.pi
C_REAL = -1.0 / (4.0 * PI)


def _phi1d2(y):
    return np.exp(-y) / y + np.sqrt(PI) / (2.0 * y**1.5) * erfc(np.sqrt(y))


def _matA_dot(xi, D, f):
    """matA(xi, D) . f for D (...,3), f broadcastable. Coincident points return 0."""
    r2 = np.sum(D * D, axis=-1)
    safe = np.where(r2 < 1e-24, 1.0, r2)
    t1 = 2.0 * xi**-0.5 * np.exp(-PI * safe / xi)
    t2 = PI * xi**-1.5 * _phi1d2(PI * safe / xi)
    rf = np.sum(D * f, axis=-1)
    out = ((t1 - t2 * safe)[..., None] * f) - (t2 * rf)[..., None] * D
    out[r2 < 1e-24] = 0.0
    return out


class Ewald3D:
    """Periodic Stokeslet sums for a 3D Bravais lattice with basis A (rows)."""

    def __init__(self, A, xi=None, cut=36.0):
        self.A = np.asarray(A, float)
        self.B, self.V0 = reciprocal(self.A)
        self.xi = float(xi) if xi is not None else self.V0 ** (2.0 / 3.0)
        kmax = np.sqrt(cut / (PI * self.xi))
        self.rmax = np.sqrt(cut * self.xi / PI)
        k = points(self.B, kmax, dims=3)[1:]  # drop k = 0
        k2 = np.sum(k * k, axis=1)
        self.k = k
        self.khat = k / np.sqrt(k2)[:, None]
        self.wk = (1 + PI * self.xi * k2) * np.exp(-PI * self.xi * k2) / k2 / (4 * PI**2 * self.V0)
        # self term: lim_{r->0} [c A_xi(r) - S(r)] = kappa I with kappa = -1/(2 pi sqrt(xi))
        # (from c A_xx = -exp(-y)/(4 pi sqrt(xi)) + erfc(sqrt(y))/(8 pi r), y = pi r^2/xi)
        self.kappa = -1.0 / (2.0 * PI * np.sqrt(self.xi))

    def velocity(self, X, P, F):
        """Velocity at points X (E,3) induced by periodic point forces F (S,3) at P (S,3).
        Coincident source/target pairs include the regular self term kappa*F."""
        X, P, F = (np.atleast_2d(np.asarray(a, float)) for a in (X, P, F))
        d = X[:, None, :] - P[None, :, :]  # (E,S,3)
        # Fourier part
        g = F[None, :, :] - self.khat[:, None, :] * np.einsum("kj,sj->ks", self.khat, F)[:, :, None]
        cosk = np.cos(2 * PI * np.einsum("kj,esj->kes", self.k, d))
        vF = np.einsum("k,kes,ksj->ej", self.wk, cosk, g)
        # real part
        rmax = self.rmax + np.max(np.linalg.norm(d, axis=-1))
        rn = points(self.A, rmax, dims=3)
        vR = np.zeros_like(vF)
        for s in range(P.shape[0]):
            D = d[:, s, None, :] - rn[None, :, :]  # (E,N,3)
            vR += C_REAL * _matA_dot(self.xi, D, F[s]).sum(axis=1)
            coinc = np.linalg.norm(d[:, s, :], axis=-1) < 1e-12
            vR[coinc] += self.kappa * F[s]
        return vF + vR
