import pytest

from spinners.backends import make_field

SMALL_L = (30, 45, 60, 90)


@pytest.fixture(scope="session")
def fields():
    """Cached fields: fields(geometry, lattice, h). Smaller cutoff radii than production, for speed."""
    cache = {}

    def get(geometry, lattice, h):
        key = (geometry, lattice, h)
        if key not in cache:
            cache[key] = make_field(geometry, lattice, h, L_over_h=SMALL_L)
        return cache[key]

    return get
