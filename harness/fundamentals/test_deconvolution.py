"""Fundamentals section 6: deconvolution (a Hogbom CLEAN on a synthetic dirty image)."""
import numpy as np


def _gaussian_psf(n, sigma=2.0):
    """A unit-peak Gaussian dirty beam centred at (n//2, n//2)."""
    y, x = np.mgrid[0:n, 0:n]
    cy = cx = n // 2
    psf = np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * sigma ** 2))
    return psf / psf.max()


def _place(psf, iy, ix):
    """Shift a PSF centred at (n//2, n//2) so its centre lands on (iy, ix)."""
    cy = cx = psf.shape[0] // 2
    return np.roll(psf, (iy - cy, ix - cx), axis=(0, 1))


def _hogbom(dirty, psf, gain=0.1, niter=500, threshold=0.005):
    residual = dirty.copy()
    model = np.zeros_like(dirty)
    for _ in range(niter):
        iy, ix = np.unravel_index(np.argmax(residual), residual.shape)
        peak = residual[iy, ix]
        if peak <= threshold:
            break
        model[iy, ix] += gain * peak
        residual -= gain * peak * _place(psf, iy, ix)
    return model, residual


def test_hogbom_recovers_two_point_sources():
    n = 64
    psf = _gaussian_psf(n, sigma=2.0)
    sources = [((20, 20), 1.0), ((42, 46), 0.5)]
    dirty = sum(flux * _place(psf, iy, ix) for (iy, ix), flux in sources)

    model, residual = _hogbom(dirty, psf)

    for (iy, ix), flux in sources:
        assert np.isclose(model[iy, ix], flux, atol=0.05)
    # the residual is cleaned down to the threshold
    assert np.max(np.abs(residual)) < 0.02
