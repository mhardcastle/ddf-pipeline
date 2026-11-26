"""Fundamentals section 4: the 2-element interferometer and Van Cittert-Zernike."""
import numpy as np


def test_point_source_lm_from_visibility_phase():
    # V(u, v) = S exp(-2 pi i (u l + v m)); recover (l, m) from the phase.
    u = np.array([10.0, -5.0, 8.0, 2.0])
    v = np.array([0.0, 6.0, -9.0, 3.0])
    l, m = 0.01, -0.02
    phi = -2 * np.pi * (u * l + v * m)
    assert np.all(np.abs(phi) < np.pi)  # chosen so there is no phase wrapping
    V = np.exp(1j * phi)
    A = np.column_stack([u, v])
    lm, *_ = np.linalg.lstsq(A, -np.angle(V) / (2 * np.pi), rcond=None)
    assert np.allclose(lm, [l, m], atol=1e-9)


def test_phase_centre_source_has_constant_real_visibility():
    S = 2.0
    u = np.linspace(-100, 100, 11)
    V = np.full_like(u, S, dtype=complex)  # l = m = 0
    assert np.allclose(V.real, S)
    assert np.allclose(V.imag, 0.0)


def test_baseline_to_uv():
    lam = 3.0
    baseline_east = 30.0
    assert np.isclose(baseline_east / lam, 10.0)
