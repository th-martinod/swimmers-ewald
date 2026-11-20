"""Periodic flow of the pairs for the two geometries of the article.

bulk3d     3D lattice of pairs in a 3D fluid. Ewald sums (velocities relative to the mean fluid velocity).
monolayer  2D lattice {c1 a1 + c2 a2} of pairs in an unbounded 3D fluid. Direct sums (fluid at rest far away).

Both expose grad_p_vs_L(P, F): the pressure gradient at the pair centre produced by the image pairs
(absolutely convergent, ~ d^-5), needed for the average over the reference sub-volume. For bulk3d the
partial sums at the cutoffs GRADP_L_3D (in units of h) are used directly, the last one as the value.
"""
import numpy as np

from .direct import LatticeSum
from .ewald import Ewald3D
from .lattices import basis, points

GEOMETRIES = ("bulk3d", "monolayer")
GRADP_L_3D = (16, 24)
# Radius of the reference ball over which the fluid velocity is averaged (Eq. 3), in units of h. The same for
# every lattice and geometry: the largest ball centred on the pair that does not overlap with the balls of the
# nearest neighbours in any of the four lattices (half the FCL nearest-neighbour distance, sqrt(2) h / 2).
R_F_OVER_H = np.sqrt(2.0) / 4


class Field:
    def __init__(self, geometry, lattice, h, xi=None, L_over_h=(40, 60, 90, 120, 160), gradp_L_3d=GRADP_L_3D):
        A = basis(lattice, h)
        self.geometry, self.lattice, self.h, self.A = geometry, lattice, float(h), A
        if geometry == "bulk3d":
            self.direct = False
            self.ewald = Ewald3D(A, xi=xi)
            self.images = LatticeSum(A, dims=3, L_over_h=gradp_L_3d, include_origin=False)
            self.d_nn = np.linalg.norm(points(A, 2 * np.linalg.norm(A, axis=1).max(), dims=3)[1])
        elif geometry == "monolayer":
            self.direct = True
            self.sums = LatticeSum(A, dims=2, L_over_h=L_over_h)
            self.images = LatticeSum(A, dims=2, L_over_h=L_over_h, include_origin=False)
            self.d_nn = np.linalg.norm(points(A, 2 * np.linalg.norm(A[:2], axis=1).max(), dims=2)[1])
        else:
            raise ValueError(f"geometry '{geometry}' not implemented")
        self.L = self.images.L
        self.r_f = R_F_OVER_H * self.h

    def velocity(self, X, P, F):
        """Velocity (E,3) at X. For the monolayer, the extrapolated direct sum."""
        if self.direct:
            from .direct import extrapolate
            return extrapolate(self.sums.L, self.sums.velocity_vs_L(X, P, F))[0]
        return self.ewald.velocity(X, P, F)

    def velocity_vs_L(self, X, P, F):
        return self.sums.velocity_vs_L(X, P, F)

    def grad_p_vs_L(self, P, F):
        """(nL, 3) image sums of the pressure gradient at the pair centre for each cutoff radius."""
        return self.images.grad_p_vs_L(np.zeros((1, 3)), P, F)[:, 0, :]


def make_field(geometry, lattice, h, **kw):
    return Field(geometry, lattice, h, **kw)
