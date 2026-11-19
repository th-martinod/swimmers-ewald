"""Free-space Green's functions applied to a force: S(r).f (mu = 1)."""
import numpy as np


def stokeslet3d_dot(D: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Oseen tensor (I/r + rr/r^3)/(8 pi) applied to f. D: (...,3), f: (3,) or (...,3).
    Coincident points (|D| = 0) give 0 (their in-plane contribution is along f only)."""
    r = np.linalg.norm(D, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rf = np.sum(D * f, axis=-1)
        out = (f / r[..., None] + D * (rf / r**3)[..., None]) / (8 * np.pi)
    out[r < 1e-12] = 0.0
    return out


def stokeslet2d_dot(D: np.ndarray, f: np.ndarray) -> np.ndarray:
    """2D Stokeslet (-ln r I + rr/r^2)/(4 pi) for line forces (per unit length), in-plane."""
    D2 = D.copy()
    D2[..., 2] = 0.0
    r = np.linalg.norm(D2, axis=-1)
    f2 = np.array(f, dtype=float, copy=True)
    f2[..., 2] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        rf = np.sum(D2 * f2, axis=-1)
        out = (-np.log(r)[..., None] * f2 + D2 * (rf / r**2)[..., None]) / (4 * np.pi)
    out[r < 1e-12] = 0.0
    return out


KERNELS = {"3d": stokeslet3d_dot, "2d": stokeslet2d_dot}
