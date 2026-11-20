"""Free-space Stokeslet (mu = 1) applied to a force f exerted on the fluid at the origin."""
import numpy as np


def stokeslet_dot(D: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Velocity (I/r + rr/r^3)/(8 pi) . f at D (...,3). Coincident points (|D| = 0) give 0."""
    r = np.linalg.norm(D, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rf = np.sum(D * f, axis=-1)
        out = (f / r[..., None] + D * (rf / r**3)[..., None]) / (8 * np.pi)
    out[r < 1e-12] = 0.0
    return out


def pressure_gradient_dot(D: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Gradient of the Stokeslet pressure p = (D.f)/(4 pi r^3):  [f/r^3 - 3 (D.f) D/r^5]/(4 pi).
    Away from the force, lap(v) = grad(p)/mu."""
    r = np.linalg.norm(D, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rf = np.sum(D * f, axis=-1)
        out = (f / r[..., None] ** 3 - 3 * D * (rf / r**5)[..., None]) / (4 * np.pi)
    out[r < 1e-12] = 0.0
    return out
