"""Fundamentals section 5: sampling, the PSF, and weighting."""
import numpy as np


def test_psf_of_complete_uv_is_a_delta():
    n = 32
    psf = np.fft.fftshift(np.fft.ifft2(np.ones((n, n)))).real
    centre = np.unravel_index(np.argmax(psf), psf.shape)
    assert np.isclose(psf[centre], 1.0)
    mask = np.ones_like(psf, dtype=bool)
    mask[centre] = False
    assert np.allclose(psf[mask], 0.0, atol=1e-9)


def test_uniform_weighting_flattens_uv_density():
    counts = np.array([1.0, 2.0, 4.0])  # visibilities per uv cell
    uniform = 1.0 / counts
    assert np.allclose(counts * uniform, 1.0)  # equal total weight per cell
    natural_total = counts.copy()               # equal per-visibility weight
    assert not np.allclose(natural_total, natural_total[0])
