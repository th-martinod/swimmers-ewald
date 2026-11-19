"""Regression against the numpy port of base_sums.rkt (production version)."""
import numpy as np
import pytest

from spinners import legacy
from spinners.backends import make_field
from spinners.pair import observables


@pytest.mark.parametrize("dim", ["2D", "3D"])
def test_legacy_hex_torque_has_inverted_sign(dim):
    nF = 40 if dim == "2D" else 20
    bug = legacy.conf_total(0.2, "hex", dim, 5.0, 1.0, nF, 4, bug=True)[0]
    ok = legacy.conf_total(0.2, "hex", dim, 5.0, 1.0, nF, 4, bug=False)[0]
    assert np.isclose(bug, -ok, rtol=1e-8)


def test_legacy_consistent_equals_corrected_bulk_times_spurious_factor():
    h = 5.0
    ok = legacy.conf_total(0.6, "hex", "3D", h, 1.0, 24, 4, bug=False)[0]
    fld, k = make_field("bulk3d", "hex", h)
    new = observables(fld, 0.6, 1.0, 0.5, k)["A_conf"]
    assert np.isclose(ok, new * h**3 / (4 * np.pi**3), rtol=1e-6)
