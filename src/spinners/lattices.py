"""Bravais lattices used in the paper (basis vectors as rows)."""
import numpy as np

SQ3 = np.sqrt(3.0)


def basis(lattice: str, h: float) -> np.ndarray:
    """Real-space basis a1, a2, a3 (rows). a3 is the rotation axis n."""
    if lattice == "hex":  # simple hexagonal, layers stacked along a3, c = h
        return np.array([[h, 0.0, 0.0], [-h / 2, SQ3 * h / 2, 0.0], [0.0, 0.0, h]])
    if lattice == "scl":
        return np.eye(3) * h
    raise ValueError(f"lattice '{lattice}' not implemented (only hex and scl)")


def reciprocal(A: np.ndarray):
    """Reciprocal basis (rows, a_i . b_j = delta_ij) and cell volume."""
    V0 = abs(np.linalg.det(A))
    return np.linalg.inv(A).T, V0


def points(A: np.ndarray, rmax: float, dims: int = 3) -> np.ndarray:
    """Lattice points with |r| <= rmax, sorted by |r| (origin included).

    dims=2 uses only a1, a2 (an in-plane monolayer)."""
    a = A[:dims]
    if dims == 3:
        V = abs(np.linalg.det(A))
        spac = [V / np.linalg.norm(np.cross(A[(i + 1) % 3], A[(i + 2) % 3])) for i in range(3)]
    else:
        area = np.linalg.norm(np.cross(A[0], A[1]))
        spac = [area / np.linalg.norm(A[1]), area / np.linalg.norm(A[0])]
    rng = [np.arange(-int(np.ceil(rmax / s)) - 1, int(np.ceil(rmax / s)) + 2) for s in spac]
    grids = np.meshgrid(*rng, indexing="ij")
    m = np.stack([g.ravel() for g in grids], 1)
    r = m @ a
    d = np.linalg.norm(r, axis=1)
    keep = d <= rmax * (1 + 1e-12)
    r, d = r[keep], d[keep]
    o = np.argsort(d, kind="stable")
    return r[o]
