"""Unit tests for the phase-invariant, direction-aligned solution comparison."""
import numpy as np

from scripts.compare_sols import direction_order, phase_invariant_diff, scalar_gains


def _g(seed, nt=2, nch=1, na=5, nd=3):
    rng = np.random.default_rng(seed)
    return np.exp(2j * np.pi * rng.random((nt, nch, na, nd)))


def test_phase_invariant_diff_removes_global_factor():
    g = _g(0)
    # a different global complex factor per (time, channel, direction)
    factor = np.exp(2j * np.pi * np.random.default_rng(1).random((2, 1, 1, 3)))
    assert phase_invariant_diff(g, g * factor) < 1e-12


def test_phase_invariant_diff_detects_a_real_change():
    g = _g(0)
    g2 = g.copy()
    g2[..., 2, :] *= 1.5  # a real per-antenna change, not a global degeneracy
    assert phase_invariant_diff(g, g2) > 0.05


def test_direction_order_matches_reordered_runs():
    ra = np.array([0.30, 0.10, 0.20])
    dec = np.array([0.00, 0.00, 0.00])
    perm = direction_order(ra, dec)
    # after ordering, ra is sorted, so two permutations of the same set align
    ra2 = ra[[2, 0, 1]]
    dec2 = dec[[2, 0, 1]]
    assert np.allclose(ra[perm], ra2[direction_order(ra2, dec2)])


def test_scalar_gains_from_6d():
    g = _g(2, nt=2, nch=1, na=4, nd=3)
    G = np.zeros((2, 1, 4, 3, 2, 2), dtype=complex)
    G[..., 0, 0] = g
    assert np.allclose(scalar_gains(G), g)
