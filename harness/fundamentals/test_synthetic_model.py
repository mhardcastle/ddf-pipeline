"""Layer B model and sols-recovery helpers (exact, no casacore)."""
import numpy as np

from harness.synthetic import visibility_model as vm
from harness.synthetic.check_sols import compare_gains, gains_from_G


def test_point_source_products_recover_gains():
    gains = vm.random_gains(6, seed=3)
    V = vm.point_source_visibilities(gains, flux=2.0)
    products = vm.normalised_gain_products(V)
    expected = np.outer(gains, np.conj(gains)) / (gains[0] * np.conj(gains[0]))
    assert np.allclose(products, expected)


def test_closure_phase_is_zero_for_point_source():
    gains = vm.random_gains(5, seed=4)
    V = vm.point_source_visibilities(gains)
    for i, j, k in [(0, 1, 2), (1, 3, 4), (0, 2, 4)]:
        assert np.isclose(vm.closure_phase(V, i, j, k), 0.0, atol=1e-12)


def test_gains_from_G_recovers_scalar_gains():
    g = np.exp(1j * np.array([0.1, -0.3, 0.7, 1.1]))
    nt, nch, nd = 3, 2, 1
    G = np.zeros((nt, nch, len(g), nd, 2, 2), dtype=complex)
    G[..., 0, 0] = g[None, None, :, None]
    G[..., 1, 1] = g[None, None, :, None]
    assert np.allclose(gains_from_G(G), g)


def test_compare_gains_removes_global_phase():
    g = vm.random_gains(5, seed=5)
    g_shifted = g * np.exp(0.7j)
    assert compare_gains(g_shifted, g) < 1e-12
