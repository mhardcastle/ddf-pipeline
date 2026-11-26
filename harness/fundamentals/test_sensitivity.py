"""Paper-derived: weighting and image sensitivity (Xie et al. 2024).

The image noise for a set of visibility weights is sqrt(sum w^2) / sum w. Natural
weighting (equal weights) minimises it; uniform weighting trades sensitivity for
resolution; the noise falls as 1/sqrt(N).
"""
import numpy as np


def image_noise(weights):
    w = np.asarray(weights, dtype=float)
    return np.sqrt(np.sum(w ** 2)) / np.sum(w)


def test_natural_weighting_is_most_sensitive():
    counts = np.array([1.0, 2.0])  # two uv cells, unequally populated
    w_nat = np.ones(int(counts.sum()))                       # every visibility equal
    w_uni = np.concatenate([np.full(int(c), 1.0 / c) for c in counts])  # equal per cell
    assert image_noise(w_nat) <= image_noise(w_uni) + 1e-12


def test_noise_scales_inverse_sqrt_n():
    n = 100
    assert np.isclose(image_noise(np.ones(n)), 1.0 / np.sqrt(n))
    assert np.isclose(image_noise(np.ones(4 * n)) * 2.0, image_noise(np.ones(n)))
