"""The average over the reference sub-volume (ball of radius sqrt(2) h / 4) against an independent Fourier sum."""
import numpy as np
import pytest

from spinners.lattices import reciprocal
from spinners.pair import observables, sources, frame, self_ball_average

PI = np.pi


def _index_cube(n, dims):
    m = np.arange(-n, n + 1)
    g = np.meshgrid(*([m] * dims), indexing="ij")
    return np.stack([x.ravel() for x in g], 1).astype(float)


def _ball_weight(k, a):
    """Average of exp(-2 pi i k.r) over the ball |r| < a: 3 (sin x - x cos x)/x^3, x = 2 pi |k| a."""
    x = 2 * PI * np.linalg.norm(k, axis=-1) * a
    return 3 * (np.sin(x) - x * np.cos(x)) / x**3


def _fourier_average(k, w, P, F):
    k2 = np.sum(k * k, axis=1)
    kh = k / np.sqrt(k2)[:, None]
    g = np.cos(2 * PI * k @ P.T) @ F
    g -= kh * np.sum(kh * g, axis=1)[:, None]
    return (w / k2) @ g


def fourier_ball_bulk(A, a, P, F, n=40):
    B, V0 = reciprocal(A)
    k = _index_cube(n, 3) @ B
    k = k[np.any(k != 0, axis=1)]
    return _fourier_average(k, _ball_weight(k, a) / (4 * PI**2 * V0), P, F)


def fourier_ball_monolayer(A, a, P, F, n=36, nq=240):
    a1, a2 = A[0], A[1]
    m = np.cross(a1, a2)
    A0 = np.linalg.norm(m)
    m /= A0
    G = np.array([a1, a2])
    kp = _index_cube(n, 2) @ np.linalg.solve(G @ G.T, G)
    x, wx = np.polynomial.legendre.leggauss(nq)
    phi, wphi = x * PI / 2, wx * PI / 2
    s = np.maximum(np.linalg.norm(kp, axis=1), 1 / np.linalg.norm(a1))
    kap = s[:, None] * np.tan(phi)[None]
    jac = s[:, None] / np.cos(phi)[None] ** 2
    k = (kp[:, None, :] + kap[..., None] * m).reshape(-1, 3)
    w = (jac * wphi[None]).reshape(-1) * _ball_weight(k, a) / (4 * PI**2 * A0)
    return _fourier_average(k, w, P, F)


@pytest.mark.parametrize("geo", ["bulk3d", "monolayer"])
@pytest.mark.parametrize("lattice", ["scl", "hex", "bcl", "fcl"])
def test_ball_average_matches_fourier_sum(fields, geo, lattice):
    fld = fields(geo, lattice, 5.0)
    th = 0.3
    P, F = sources(th, 1.0, 0.5)
    _, e = frame(th)
    ref = (fourier_ball_bulk if geo == "bulk3d" else fourier_ball_monolayer)(fld.A, fld.r_f, P, F) @ e
    assert np.isclose(observables(fld, th, 1.0, 0.5)["B_perp"], ref, rtol=2e-3)


def test_self_average_of_the_isolated_pair():
    """Own pair: <v>_ball = -zeta V^S Q R e / (10 pi mu a^3) (forces on the fluid, Eq. 2)."""
    for th, R, a in ((0.0, 0.5, 2.5), (1.1, 0.25, 1.7)):
        P, F = sources(th, 1.0, R)
        _, e = frame(th)
        assert np.allclose(self_ball_average(P, F, a), -R * e / (10 * PI * a**3), atol=1e-16)


@pytest.mark.parametrize("geo", ["bulk3d", "monolayer"])
@pytest.mark.parametrize("lattice", ["scl", "hex", "bcl", "fcl"])
def test_reference_ball_is_common_and_isolates_the_pair(fields, geo, lattice):
    """Same radius sqrt(2) h / 4 in every lattice. For Q/h <= 0.25 it contains the four point forces of the
    pair and none of the neighbours (needed by the mean-value evaluation), and it does not overlap with the
    balls of the nearest neighbours."""
    h = 4.0                                            # Q/h = 0.25, the largest value of the article
    fld = fields(geo, lattice, h)
    assert np.isclose(fld.r_f, np.sqrt(2) * h / 4)
    assert 2 * fld.r_f <= fld.d_nn * (1 + 1e-12)
    for th in (0.0, 0.4, 1.3):
        for R in (0.5, 0.25):
            P, _ = sources(th, 1.0, R)
            reach = np.linalg.norm(P, axis=1).max()
            assert reach < fld.r_f
            assert fld.d_nn - reach > fld.r_f
