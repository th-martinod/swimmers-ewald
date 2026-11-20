import numpy as np
import pytest

from spinners.ewald import Ewald3D
from spinners.kernels import stokeslet_dot
from spinners.lattices import basis, points
from spinners.pair import sources, eval_points, frame


@pytest.mark.parametrize("lattice", ["hex", "scl", "bcl", "fcl"])
def test_xi_independence(lattice):
    """The Ewald sum must not depend on the splitting parameter."""
    h = 5.0
    A = basis(lattice, h)
    P, F = sources(0.3, 1.0, 0.5)
    X = np.vstack([eval_points(0.3, 1.0, 0.5), [[1.3, -0.7, 0.4]]])
    V0 = abs(np.linalg.det(A))
    ref = Ewald3D(A, xi=V0 ** (2 / 3)).velocity(X, P, F)
    for s in (0.4, 2.5):
        v = Ewald3D(A, xi=s * V0 ** (2 / 3)).velocity(X, P, F)
        assert np.allclose(v, ref, rtol=0, atol=1e-12 * np.abs(ref).max())


def test_ewald_matches_direct_3d_sum():
    """T = V_conf . e from Ewald equals the absolutely convergent direct 3D lattice sum (~ d^-5)."""
    h, th, R = 5.0, 0.6, 0.5
    A = basis("hex", h)
    P, F = sources(th, 1.0, R)
    X = eval_points(th, 1.0, R)[:4]
    _, e = frame(th)
    v = Ewald3D(A).velocity(X, P, F)
    T_ewald = np.dot((v[0] - v[1]) - (v[2] - v[3]), e)
    rn = points(A, 14 * h, dims=3)
    Ls, vals = [], []
    for L in (8, 10, 12, 14):
        r = rn[np.linalg.norm(rn, axis=1) <= L * h + 1e-9]
        w = np.array([sum(stokeslet_dot(x - P[s] - r, F[s]).sum(0) for s in range(4)) for x in X])
        vals.append(np.dot((w[0] - w[1]) - (w[2] - w[3]), e))
        Ls.append(L * h)
    Ls, vals = np.array(Ls), np.array(vals)
    M = np.stack([np.ones_like(Ls), Ls**-2, Ls**-3], 1)
    a = np.linalg.lstsq(M, vals, rcond=None)[0][0]
    assert abs(T_ewald - a) < 2e-4 * abs(T_ewald)
