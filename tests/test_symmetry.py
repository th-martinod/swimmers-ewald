"""Conventions, symmetries and analytical results of the article (Section III and SI Section S3)."""
import numpy as np
import pytest

from spinners.lattices import basis
from spinners.pair import observables, sources, eval_points, frame

GEOS = ["bulk3d", "monolayer"]
LATS = ["hex", "scl", "bcl", "fcl"]


def test_forces_on_the_fluid_spin_r1_counterclockwise():
    """Convention of Fig. 1: the forces exerted on the fluid turn the fluid around r1 = -Q/2
    counterclockwise and around r2 = +Q/2 clockwise (torque along +n and -n)."""
    th = 0.37
    qh, _ = frame(th)
    P, F = sources(th, 1.0, 0.5)
    tz1 = sum(np.cross(P[s] + qh / 2, F[s])[2] for s in (2, 3))
    tz2 = sum(np.cross(P[s] - qh / 2, F[s])[2] for s in (0, 1))
    assert tz1 > 0 and tz2 < 0


@pytest.mark.parametrize("geo", GEOS)
@pytest.mark.parametrize("lattice", LATS)
def test_centres_are_advected_with_the_same_velocity(fields, geo, lattice):
    """v(r1) = v(r2) on any Bravais lattice (inversion symmetry about the pair centre)."""
    fld = fields(geo, lattice, 6.0)
    for th in (0.11, -0.37, 0.52):
        for R in (0.5, 0.25):
            o = observables(fld, th, 1.0, R)
            assert abs(o["Qdot_perp"]) < 1e-12 * max(1.0, abs(o["B_perp"]))


@pytest.mark.parametrize("geo", GEOS)
def test_field_is_even_about_pair_centre(fields, geo):
    fld = fields(geo, "hex", 6.0)
    P, F = sources(0.23, 1.0, 0.5)
    X = np.array([[0.31, 0.17, 0.0], [0.8, -0.4, 0.0], [1.7, 0.9, 0.0]])
    v = fld.velocity(np.vstack([X, -X]), P, F)
    assert np.abs(v[:3] - v[3:]).max() < 1e-12 * np.abs(v).max()


@pytest.mark.parametrize("geo", GEOS)
def test_vconf_is_twice_the_velocity_difference_across_one_spinner(fields, geo):
    """By inversion symmetry, V_conf = 2 [v(A) - v(B)] (SI, Section S2)."""
    fld = fields(geo, "hex", 6.0)
    th, R = 0.29, 0.5
    _, e = frame(th)
    P, F = sources(th, 1.0, R)
    v = fld.velocity(eval_points(th, 1.0, R), P, F)
    o = observables(fld, th, 1.0, R)
    assert np.isclose(o["T"], 2 * np.dot(v[0] - v[1], e), rtol=1e-8, atol=0)


@pytest.mark.parametrize("geo", GEOS)
def test_hexagonal_lattice_is_stable_around_theta_zero(fields, geo):
    """A_perp < 0 for 0 < |theta| < pi/6 and A_perp > 0 for pi/6 < |theta| < pi/4 (SI, Section S3)."""
    fld = fields(geo, "hex", 5.0)
    for th in (0.05, 0.2, 0.4, 0.5):
        assert observables(fld, th, 1.0, 0.5)["A_perp"] < 0
        assert observables(fld, -th, 1.0, 0.5)["A_perp"] < 0
    for th in (0.56, 0.7):
        assert observables(fld, th, 1.0, 0.5)["A_perp"] > 0


@pytest.mark.parametrize("geo", GEOS)
def test_hex_symmetries(fields, geo):
    """A_perp is even in theta, T and B_perp have period pi/3."""
    fld = fields(geo, "hex", 5.0)
    for th in (0.1, 0.33):
        a, b, c = (observables(fld, t, 1.0, 0.5) for t in (th, -th, th + np.pi / 3))
        assert np.isclose(a["A_perp"], b["A_perp"], rtol=1e-8)
        assert np.isclose(a["T"], c["T"], rtol=1e-8)
        assert np.isclose(a["B_perp"], c["B_perp"], rtol=1e-6)


@pytest.mark.parametrize("geo,c6", [("monolayer", 0.927), ("bulk3d", 0.955)])
def test_six_fold_amplitude_of_section_S3(fields, geo, c6):
    """T_hex = c6 * 126 K(h) sin(6 theta) with K = 45 F Q^2 R^2 / (64 pi mu h^5), Q/h = 0.2."""
    h, R = 5.0, 0.5
    K = 45 * R**2 / (64 * np.pi * h**5)
    fld = fields(geo, "hex", h)
    for th in (0.1, 0.25):
        assert np.isclose(observables(fld, th, 1.0, R)["T"] / (126 * K * np.sin(6 * th)), c6, atol=2e-3)


@pytest.mark.parametrize("geo", GEOS)
def test_square_cancellation_only_for_touching_spinners(fields, geo):
    fld = fields(geo, "scl", 5.0)
    ref = abs(observables(fields(geo, "hex", 5.0), 0.2, 1.0, 0.5)["T"])
    for th in (0.1, 0.3):
        assert abs(observables(fld, th, 1.0, 0.5)["T"]) < 1e-10 * ref
    assert abs(observables(fld, 0.3, 1.0, 0.25)["T"]) > 1e-4 * ref


@pytest.mark.parametrize("lattice", ["bcl", "fcl"])
def test_tilted_lattices_have_a3_along_z_and_theta_from_projected_a1(lattice):
    A = basis(lattice, 5.0)
    assert np.allclose(A[2][:2], 0) and A[2][2] > 0
    assert abs(A[0][1]) < 1e-12 and A[0][0] > 0        # projection of a1 along +x
    assert abs(np.dot(A[0], A[2])) > 1e-6               # a1 is not perpendicular to a3
