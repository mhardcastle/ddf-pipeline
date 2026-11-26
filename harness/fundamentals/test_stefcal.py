"""Paper-derived: StEFCal / complex gain calibration (Smirnov & Tasse 2015).

We pin the StEFCal update itself: starting from the true gains, one sweep returns
them (the true solution is a fixed point), and the solution is degenerate in a
global phase. Convergence rate is numerical and is not asserted.
"""
import numpy as np


def _hermitian_model(n, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    m = a + a.conj().T
    np.fill_diagonal(m, np.abs(np.diag(m)))
    return m


def stefcal_sweep(V, M, g):
    """One StEFCal sweep: g_i = sum_j V_ij g_j conj(M_ij) / sum_j |g_j|^2 |M_ij|^2."""
    n = V.shape[0]
    g = np.asarray(g, dtype=complex).copy()
    for i in range(n):
        js = np.arange(n) != i
        num = np.sum(V[i, js] * g[js] * np.conj(M[i, js]))
        den = np.sum(np.abs(g[js]) ** 2 * np.abs(M[i, js]) ** 2)
        if den:
            g[i] = num / den
    return g


def test_stefcal_update_fixes_the_true_solution():
    n = 6
    rng = np.random.default_rng(1)
    g_true = (1.0 + 0.1 * rng.normal(size=n)) * np.exp(2j * np.pi * rng.random(n))
    M = _hermitian_model(n)
    V = g_true[:, None] * np.conj(g_true)[None, :] * M

    g_next = stefcal_sweep(V, M, g_true)

    # the true gains (up to a global phase) are a fixed point of the update
    assert np.allclose(g_next / g_next[0], g_true / g_true[0], atol=1e-10)


def test_gain_solution_is_degenerate_in_global_phase():
    # g and g e^{i phi} give exactly the same visibilities, so the global phase is free.
    n = 5
    rng = np.random.default_rng(2)
    g = np.exp(2j * np.pi * rng.random(n))
    M = _hermitian_model(n, seed=2)
    V1 = g[:, None] * np.conj(g)[None, :] * M
    g2 = g * np.exp(0.7j)
    V2 = g2[:, None] * np.conj(g2)[None, :] * M
    assert np.allclose(V1, V2)
