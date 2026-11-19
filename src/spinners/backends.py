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
    if geometry in ("monolayer", "pure2d"):
        kern = "3d" if geometry == "monolayer" else "2d"
        d = Direct2DLattice(A, kernel=kern, L_over_h=L_over_h)

        def field(X, P, F):
            return d.velocity(X, P, F)

        field.vs_L = d.velocity_vs_L  # lets pair.observables extrapolate each observable directly
        field.L = d.L
        return field, kern
    raise ValueError(geometry)
