"""Direct real-space lattice sums with extrapolation of the cutoff radius.

Used for the monolayer (a two-dimensional lattice of pairs in an unbounded 3D fluid) and for the
image sums of the pressure gradient. Every quantity is built from force-free combinations whose
lattice sums converge absolutely once shells are summed symmetrically (cutoff by radius on a
centrosymmetric lattice). The cutoff is extrapolated to infinity with a fit a + b/L + c/L^2.
"""
import numpy as np

from .kernels import stokeslet_dot, pressure_gradient_dot
from .lattices import points


class LatticeSum:
    """Partial sums over the lattice points r_n (sorted by radius) for several cutoff radii L.

    dims = 2: monolayer {c1 a1 + c2 a2}.  dims = 3: bulk lattice.
    include_origin = False drops r_n = 0 (sum over the images of the central pair only)."""

    def __init__(self, A, dims=2, L_over_h=(40, 60, 90, 120, 160), include_origin=True):
        self.A = np.asarray(A, float)
        h = np.cbrt(abs(np.linalg.det(self.A)))  # lattice length scale
        self.L = np.array(L_over_h, float) * h
        rn = points(self.A, self.L.max(), dims=dims)
        self.rn = rn if include_origin else rn[1:]
        rad = np.linalg.norm(self.rn, axis=1)
        self.cuts = np.searchsorted(rad, self.L * (1 + 1e-12), side="right")

    def _vs_L(self, kernel, X, P, F):
        X, P, F = (np.atleast_2d(np.asarray(a, float)) for a in (X, P, F))
        out = np.zeros((len(self.L), X.shape[0], 3))
        for e in range(X.shape[0]):
            contrib = np.zeros((len(self.rn), 3))
            for s in range(P.shape[0]):
                contrib += kernel(X[e] - P[s] - self.rn, F[s])
            out[:, e, :] = np.cumsum(contrib, axis=0)[self.cuts - 1]
        return out

    def velocity_vs_L(self, X, P, F):
        """(nL, E, 3) partial sums of the velocity at X for each cutoff L."""
        return self._vs_L(stokeslet_dot, X, P, F)

    def grad_p_vs_L(self, X, P, F):
        """(nL, E, 3) partial sums of the pressure gradient at X for each cutoff L."""
        return self._vs_L(pressure_gradient_dot, X, P, F)


def extrapolate(L, vals):
    """Fit vals(L) = a + b/L + c/L^2 per component. Returns (a, |vals[-1] - a|)."""
    M = np.stack([np.ones_like(L), 1 / L, 1 / L**2], 1)
    flat = vals.reshape(len(L), -1)
    coef, *_ = np.linalg.lstsq(M, flat, rcond=None)
    a = coef[0].reshape(vals.shape[1:])
    return a, np.abs(vals[-1] - a)
