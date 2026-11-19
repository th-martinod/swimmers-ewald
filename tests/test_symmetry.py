import numpy as np
import pytest

from spinners.pair import observables, sources, eval_points, frame, isolated_pair_velocity

GEOS = ["bulk3d", "monolayer", "pure2d"]


@pytest.mark.parametrize("geo", GEOS)
@pytest.mark.parametrize("lattice", ["hex", "scl"])
def test_qdot_vanishes_on_any_bravais_lattice(fields, geo, lattice):
    """v(r1) = v(r2) exactly: the field is even about the pair centre (lattice inversion
    symmetry + evenness of the pair's own field), so Qdot has no perpendicular part."""
    fld, k = fields(geo, lattice, 6.0)
    for th in (0.11, -0.37, 0.52):
        for R in (0.5, 0.25):
            o = observables(fld, th, 1.0, R, k)
            assert abs(o["omega_qdot"]) < 1e-12 * max(1.0, abs(o["U_perp"]))


@pytest.mark.parametrize("geo", GEOS)
def test_field_is_even_about_pair_centre(fields, geo):
    fld, _ = fields(geo, "hex", 6.0)
    P, F = sources(0.23, 1.0, 0.5)
    X = np.array([[0.31, 0.17, 0.0], [0.8, -0.4, 0.0], [1.7, 0.9, 0.0]])
    v, _ = fld(np.vstack([X, -X]), P, F)
    assert np.abs(v[:3] - v[3:]).max() < 1e-12 * np.abs(v).max()


@pytest.mark.parametrize("geo", GEOS)
def test_vconf_is_normal_strain_across_spinner(fields, geo):
    """rate(V_conf) = -2 (v(p1) - v(p2)).e_perp : an extension rate, not a rotation of Q."""
    fld, k = fields(geo, "hex", 6.0)
    th, R = 0.29, 0.5
    qh, e = frame(th)
    P, F = sources(th, 1.0, R)
    v, _ = fld(eval_points(th, 1.0, R), P, F)
    o = observables(fld, th, 1.0, R, k)
    assert np.isclose(o["rate_conf"], -2 * np.dot(v[0] - v[1], e), rtol=1e-10, atol=0)


@pytest.mark.parametrize("geo", GEOS)
def test_hex_symmetries(fields, geo):
    fld, k = fields(geo, "hex", 5.0)
    for th in (0.1, 0.33):
        a = observables(fld, th, 1.0, 0.5, k)
        b = observables(fld, -th, 1.0, 0.5, k)
        c = observables(fld, th + np.pi / 3, 1.0, 0.5, k)
        assert np.isclose(a["A_conf"], b["A_conf"], rtol=1e-8)
        assert np.isclose(a["rate_conf"], c["rate_conf"], rtol=1e-8)
        assert np.isclose(a["U_perp"], c["U_perp"], rtol=1e-8)


@pytest.mark.parametrize("geo", GEOS)
def test_square_cancellation_only_for_touching_spinners(fields, geo):
    fld, k = fields(geo, "scl", 5.0)
    fh, kh = fields(geo, "hex", 5.0)
    ref = abs(observables(fh, 0.2, 1.0, 0.5, kh)["rate_conf"])
    for th in (0.1, 0.3):
        assert abs(observables(fld, th, 1.0, 0.5, k)["rate_conf"]) < 1e-10 * ref
    assert abs(observables(fld, 0.3, 1.0, 0.25, k)["rate_conf"]) > 1e-4 * ref


def test_own_pair_gives_no_rotation():
    th, R = 0.4, 0.5
    qh, e = frame(th)
    X = eval_points(th, 1.0, R)
    for kern in ("3d", "2d"):
        v = isolated_pair_velocity(th, 1.0, R, X, kern)
        Vc = (v[0] - v[1]) - (v[2] - v[3])
        assert abs(np.dot(np.cross(qh, Vc), [0, 0, 1])) < 1e-14
        assert abs(np.dot(np.cross(qh, v[5] - v[4]), [0, 0, 1])) < 1e-14
