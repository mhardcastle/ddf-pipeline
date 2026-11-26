"""Known-answer visibility/gain model for the calibration harness (layer B).

Pure numpy: a point source at the phase centre gives V_ij = g_i conj(g_j) S, so
the true gains (and the gain products, up to a global degeneracy) are known.
The casacore injection into a real MS lives in inject_gains.py (run on Jean Zay).
"""
import numpy as np


def random_gains(n_ant, seed=0, amp_spread=0.1):
    """Complex gains with |g| ~ 1 +/- amp_spread and random phases."""
    rng = np.random.default_rng(seed)
    amp = 1.0 + amp_spread * rng.normal(size=n_ant)
    phase = 2j * np.pi * rng.random(n_ant)
    return amp * np.exp(phase)


def point_source_visibilities(gains, flux=1.0):
    """V_ij = g_i conj(g_j) S for a point source at the phase centre."""
    g = np.asarray(gains, dtype=complex)
    return flux * np.outer(g, np.conj(g))


def normalised_gain_products(V):
    """V_ij / V_00: the observable gain products, up to the global scale."""
    return V / V[0, 0]


def closure_phase(V, i, j, k):
    return np.angle(V[i, j] * V[j, k] * np.conj(V[i, k]))
