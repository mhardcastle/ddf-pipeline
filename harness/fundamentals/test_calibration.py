"""Fundamentals section 8: 1GC observables and the gain-product rank-1 structure."""
import numpy as np


def test_gain_products_are_rank_one():
    # For one point source, V_ij = g_i g_j* S, so V/S = g g^H is rank 1.
    S = 2.0
    g = np.exp(1j * np.array([0.1, -0.3, 0.7, 1.1, -0.5]))
    V = np.outer(g, np.conj(g)) * S
    assert np.linalg.matrix_rank(V, tol=1e-9) == 1
    U, s, _ = np.linalg.svd(V)
    g_hat = U[:, 0] * np.sqrt(s[0])
    ratios = g_hat / g
    assert np.allclose(ratios / ratios[0], 1.0, atol=1e-9)  # recovered up to a global phase


def test_closure_phase_of_a_single_source_is_zero():
    # V_ij = g_i g_j* on every baseline has zero closure phase.
    g = np.exp(1j * np.array([0.2, -0.4, 0.9, 1.3]))
    V = np.outer(g, np.conj(g))
    c = np.angle(V[0, 1] * V[1, 2] * np.conj(V[0, 2]))
    assert np.isclose(c, 0.0, atol=1e-12)
