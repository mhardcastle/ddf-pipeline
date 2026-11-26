"""Fundamentals section 7: Jones algebra, the RIME, and closure invariants."""
import numpy as np


def _scalar_jones(g):
    return np.eye(2, dtype=complex) * g


def test_rime_with_scalar_gains():
    g1 = 1.2 * np.exp(0.3j)
    g2 = 0.8 * np.exp(-0.5j)
    B = 2.0
    V = (_scalar_jones(g1) @ (np.eye(2) * B) @ _scalar_jones(g2).conj().T)[0, 0]
    assert np.isclose(V, g1 * np.conj(g2) * B)


def test_closure_phase_is_invariant_to_antenna_gains():
    rng = np.random.default_rng(0)
    n_ant = 4
    Vtrue = rng.normal(size=(n_ant, n_ant)) + 1j * rng.normal(size=(n_ant, n_ant))
    gains = np.exp(2j * np.pi * rng.random(n_ant))
    Vcal = Vtrue * gains[:, None] * np.conj(gains[None, :])
    for _ in range(10):
        i, j, k = rng.choice(n_ant, size=3, replace=False)
        c_true = np.angle(Vtrue[i, j] * Vtrue[j, k] * np.conj(Vtrue[i, k]))
        c_cal = np.angle(Vcal[i, j] * Vcal[j, k] * np.conj(Vcal[i, k]))
        assert np.isclose(c_true, c_cal)


def test_closure_amplitude_is_invariant_to_antenna_gains():
    rng = np.random.default_rng(1)
    n_ant = 4
    Vtrue = rng.random((n_ant, n_ant)) + 0.1
    gains = 1.0 + rng.random(n_ant)
    Vcal = Vtrue * gains[:, None] * gains[None, :]
    for _ in range(10):
        i, j, k, l = rng.choice(n_ant, size=4, replace=False)
        a_true = (Vtrue[i, j] * Vtrue[k, l]) / (Vtrue[i, k] * Vtrue[j, l])
        a_cal = (Vcal[i, j] * Vcal[k, l]) / (Vcal[i, k] * Vcal[j, l])
        assert np.isclose(a_true, a_cal)
