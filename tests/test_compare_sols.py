"""Unit tests for the phase-invariant solution comparison (no casacore)."""
import numpy as np

from scripts.compare_sols import gains_from_G, phase_invariant_diff


def _g(seed, nt=2, nch=1, na=5):
    rng = np.random.default_rng(seed)
    return np.exp(2j * np.pi * rng.random((nt, nch, na)))


def test_phase_invariant_diff_removes_global_phase():
    g = _g(0)
    phase = np.exp(2j * np.pi * np.random.default_rng(1).random((2, 1, 1)))
    assert phase_invariant_diff(g, g * phase) < 1e-12


def test_phase_invariant_diff_detects_a_real_change():
    g = _g(0)
    g2 = g.copy()
    g2[..., 2] *= 1.5  # a real per-antenna change, not a global degeneracy
    assert phase_invariant_diff(g, g2) > 0.05


def test_gains_from_G_scalar():
    g = _g(2)
    nt, nch, na = g.shape
    G = np.zeros((nt, nch, na, 1, 2, 2), dtype=complex)
    G[..., 0, 0] = g[..., None]
    assert np.allclose(gains_from_G(G), g)
