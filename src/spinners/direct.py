"""Direct real-space lattice sums for 2D lattices (monolayer in 3D fluid, or pure 2D).

All quantities used here are built from force-free combinations whose lattice sums converge
absolutely once shells are summed symmetrically (cutoff by radius on a centrosymmetric lattice).
The cutoff is extrapolated to infinity with a least-squares fit a + b/L + c/L^2.
"""
import numpy as np

from .kernels import KERNELS
from .lattices import points


class Direct2DLattice:
    def __init__(self, A, kernel="3d", L_over_h=(40, 60, 90, 120, 160)):
        self.A = np.asarray(A, float)
        self.kernel = KERNELS[kernel]
        h = np.cbrt(abs(np.linalg.det(self.A)))  # lattice length scale
        self.L = np.array(L_over_h, float) * h
        # pure 2D: line forces along z through the {a1, a2} net -> distances measured in the xy plane
        metric = (lambda r: r * np.array([1.0, 1.0, 0.0])) if kernel == "2d" else None
        self.rn = points(self.A, self.L.max(), dims=2, metric=metric)
        self.rad = np.linalg.norm(self.rn if metric is None else metric(self.rn), axis=1)
        self.cuts = np.searchsorted(self.rad, self.L * (1 + 1e-12), side="right")

    def velocity_vs_L(self, X, P, F):
        """(nL, E, 3) partial sums for each cutoff L (own cell always included)."""
        X, P, F = (np.atleast_2d(np.asarray(a, float)) for a in (X, P, F))
        E = X.shape[0]
        out = np.zeros((len(self.L), E, 3))
        for e in range(E):
            contrib = np.zeros((len(self.rn), 3))
            for s in range(P.shape[0]):
                D = X[e] - P[s] - self.rn
                contrib += self.kernel(D, F[s])
            cs = np.cumsum(contrib, axis=0)
            out[:, e, :] = cs[self.cuts - 1]
        return out

    def velocity(self, X, P, F):
        """Extrapolated velocity (E,3) and an error estimate (E,3)."""
        vs = self.velocity_vs_L(X, P, F)
        return extrapolate(self.L, vs)


def extrapolate(L, vals):
    """Fit vals(L) = a + b/L + c/L^2 per component. Returns (a, |vals[-1] - a|)."""
    M = np.stack([np.ones_like(L), 1 / L, 1 / L**2], 1)
    flat = vals.reshape(len(L), -1)
    coef, *_ = np.linalg.lstsq(M, flat, rcond=None)
    a = coef[0].reshape(vals.shape[1:])
    return a, np.abs(vals[-1] - a)
