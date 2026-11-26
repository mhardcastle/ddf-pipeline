"""Fundamentals section 2: Fourier, convolution, sampling (exact assertions)."""
import numpy as np


def test_dft_of_complex_exponential_is_a_delta():
    n, k0 = 64, 5
    x = np.exp(2j * np.pi * k0 * np.arange(n) / n)
    X = np.fft.fft(x)
    assert np.argmax(np.abs(X)) == k0
    assert np.isclose(np.abs(X[k0]), n)
    assert np.allclose(np.delete(np.abs(X), k0), 0.0, atol=1e-9)


def test_convolution_theorem_circular():
    rng = np.random.default_rng(0)
    a = rng.normal(size=32) + 1j * rng.normal(size=32)
    b = rng.normal(size=32) + 1j * rng.normal(size=32)
    via_fft = np.fft.ifft(np.fft.fft(a) * np.fft.fft(b))
    direct = np.array([sum(a[j] * b[(i - j) % len(a)] for j in range(len(a)))
                       for i in range(len(a))])
    assert np.allclose(via_fft, direct, atol=1e-9)


def test_parseval():
    rng = np.random.default_rng(1)
    x = rng.normal(size=64) + 1j * rng.normal(size=64)
    assert np.isclose(np.sum(np.abs(x) ** 2),
                      np.sum(np.abs(np.fft.fft(x)) ** 2) / len(x))


def test_resolvable_frequency_is_recovered():
    n, f = 128, 0.1  # cycles/sample, below Nyquist (0.5)
    x = np.exp(2j * np.pi * f * np.arange(n))
    peak = np.argmax(np.abs(np.fft.fft(x)))
    assert np.isclose(peak / n, f, atol=1.0 / n)
