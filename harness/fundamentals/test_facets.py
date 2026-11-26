"""Paper-derived: direction cosines and facet coverage (Tasse et al. 2017; van Weeren 2016).

The full facet tesselation is exercised by the acceptance/e2e runs; here we pin the
geometry it relies on: the direction-cosine transform and the field coverage.
"""
import numpy as np


def lm_from_radec(ra, dec, ra0, dec0):
    dra = ra - ra0
    l = np.cos(dec) * np.sin(dra)
    m = np.sin(dec) * np.cos(dec0) - np.cos(dec) * np.sin(dec0) * np.cos(dra)
    return l, m


def radec_from_lm(l, m, ra0, dec0):
    n = np.sqrt(1.0 - l * l - m * m)
    dec = np.arcsin(m * np.cos(dec0) + n * np.sin(dec0))
    ra = ra0 + np.arctan2(l, n * np.cos(dec0) - m * np.sin(dec0))
    return ra, dec


def test_lm_radec_roundtrip():
    ra0, dec0 = 0.3, 0.9
    for ra, dec in [(0.34, 0.93), (0.27, 0.88), (0.30, 0.90)]:
        l, m = lm_from_radec(ra, dec, ra0, dec0)
        ra2, dec2 = radec_from_lm(l, m, ra0, dec0)
        assert np.isclose(ra2, ra, atol=1e-12)
        assert np.isclose(dec2, dec, atol=1e-12)


def test_regular_facet_grid_tiles_the_field():
    n_side = 4
    lmax = mmax = 0.2
    edges_l = np.linspace(-lmax, lmax, n_side + 1)
    edges_m = np.linspace(-mmax, mmax, n_side + 1)
    areas = np.diff(edges_l)[:, None] * np.diff(edges_m)[None, :]
    # no gaps, no overlap: the facets tile exactly the field
    assert np.isclose(areas.sum(), (2 * lmax) * (2 * mmax))
    # every facet centre lies inside the field
    centres_l = 0.5 * (edges_l[:-1] + edges_l[1:])
    centres_m = 0.5 * (edges_m[:-1] + edges_m[1:])
    assert centres_l.min() > -lmax and centres_l.max() < lmax
    assert centres_m.min() > -mmax and centres_m.max() < mmax
