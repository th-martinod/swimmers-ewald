"""Bravais lattices used in the paper (basis vectors as rows)."""
import numpy as np

SQ3 = np.sqrt(3.0)


def _rotate_a3_to_z(A: np.ndarray) -> np.ndarray:
    """Rigidly rotate a basis so that a3 -> z and the projection of a1 on the plane
    perpendicular to a3 -> +x. Then n = a3/|a3| = z and theta is measured from that projection."""
    z = A[2] / np.linalg.norm(A[2])
    x = A[0] - np.dot(A[0], z) * z
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return A @ np.array([x, y, z]).T


def basis(lattice: str, h: float) -> np.ndarray:
    """Real-space basis a1, a2, a3 (rows), oriented so that a3 (the rotation axis n) is along z.

    hex: simple hexagonal, layers stacked along a3 with spacing h (h = nearest-neighbour distance).
    scl: simple cubic of side h.
    bcl, fcl: primitive vectors of Hasimoto (1959) / base_sums.rkt with cube side h, so V0 = h^3/2 and h^3/4.
      For these a1, a2 are not perpendicular to a3, so the basis is rotated (a3 -> z) and theta is
      measured from the projection of a1 on the plane perpendicular to a3 (article definition)."""
    if lattice == "hex":
        return np.array([[h, 0.0, 0.0], [-h / 2, SQ3 * h / 2, 0.0], [0.0, 0.0, h]])
    if lattice == "scl":
        return np.eye(3) * h
    if lattice == "bcl":
        return _rotate_a3_to_z(np.array([[1.0, 1.0, -1.0], [-1.0, 1.0, 1.0], [1.0, -1.0, 1.0]]) * h / 2)
    if lattice == "fcl":
        return _rotate_a3_to_z(np.array([[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]]) * h / 2)
    raise ValueError(f"lattice '{lattice}' not implemented")


def reciprocal(A: np.ndarray):
    """Reciprocal basis (rows, a_i . b_j = delta_ij) and cell volume."""
    V0 = abs(np.linalg.det(A))
    return np.linalg.inv(A).T, V0


def points(A: np.ndarray, rmax: float, dims: int = 3, metric=None) -> np.ndarray:
    """Lattice points with |metric(r)| <= rmax, sorted by that radius (origin included).

    dims=2 uses only a1, a2 (a monolayer, tilted with respect to z for bcl and fcl).
    metric: optional function r -> r' used to measure distances (e.g. projection on the xy plane)."""
    a = A[:dims]
    if dims == 3:
        V = abs(np.linalg.det(A))
        spac = [V / np.linalg.norm(np.cross(A[(i + 1) % 3], A[(i + 2) % 3])) for i in range(3)]
    else:
        area = np.linalg.norm(np.cross(A[0], A[1]))
        spac = [area / np.linalg.norm(A[1]), area / np.linalg.norm(A[0])]
    pad = 1.0 if metric is None else 3.0
    rng = [np.arange(-int(np.ceil(pad * rmax / s)) - 1, int(np.ceil(pad * rmax / s)) + 2) for s in spac]
    grids = np.meshgrid(*rng, indexing="ij")
    m = np.stack([g.ravel() for g in grids], 1)
    r = m @ a
    d = np.linalg.norm(r if metric is None else metric(r), axis=1)
    keep = d <= rmax * (1 + 1e-12)
    r, d = r[keep], d[keep]
    o = np.argsort(d, kind="stable")
    return r[o]
