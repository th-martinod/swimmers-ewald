"""Field back-ends for the three geometries."""
import numpy as np

from .ewald import Ewald3D
from .direct import Direct2DLattice
from .lattices import basis

GEOMETRIES = ("bulk3d", "monolayer", "pure2d")


def make_field(geometry, lattice, h, xi=None, L_over_h=(40, 60, 90, 120, 160)):
    A = basis(lattice, h)
    if geometry == "bulk3d":
        ew = Ewald3D(A, xi=xi)
        return (lambda X, P, F: (ew.velocity(X, P, F), np.zeros((len(X), 3)))), "3d"
    if geometry == "monolayer":
        d = Direct2DLattice(A, kernel="3d", L_over_h=L_over_h)
        return d.velocity, "3d"
    if geometry == "pure2d":
        d = Direct2DLattice(A, kernel="2d", L_over_h=L_over_h)
        return d.velocity, "2d"
    raise ValueError(geometry)
