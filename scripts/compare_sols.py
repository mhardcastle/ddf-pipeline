#!/usr/bin/env python3
"""Compare two killMS .sols.npz up to the self-calibration degeneracy.

The gains are defined only up to a global complex factor per solution slot
(Smirnov & Tasse 2015). Also, the facet directions can be produced in a different
order between runs, so comparing ``direction 0`` element-wise compares different
sky directions. This loads the direction positions (ClusterCat/SkyModel), orders
both runs by (ra, dec), checks the positions agree, and then reports the
phase/amplitude-invariant gain difference max |g_a/g_a[0] - g_b/g_b[0]|.

If the directions (or their count) differ beyond the tolerance, the solutions are
not comparable and that is reported instead of a misleading number.
"""
import argparse
import sys

import numpy as np


def load(path):
    d = np.load(path, allow_pickle=False)
    sols = d["Sols"]
    names = sols.dtype.names or ()
    G = np.asarray(sols["G"]) if "G" in names else np.asarray(d["Gains"])
    cat = d["ClusterCat"] if "ClusterCat" in d else d["SkyModel"]
    ra = np.asarray(cat["ra"], dtype=float)
    dec = np.asarray(cat["dec"], dtype=float)
    return G, ra, dec


def scalar_gains(G):
    """(nt, nch, na, nd) complex scalar gains from Sols.G (or already 4-D)."""
    G = np.asarray(G, dtype=complex)
    if G.ndim == 6:
        return G[..., 0, 0]
    if G.ndim == 4:
        return G
    raise ValueError("unexpected Sols.G shape %r" % (G.shape,))


def direction_order(ra, dec):
    return np.lexsort((dec, ra))


def phase_invariant_diff(g_a, g_b):
    """g_a, g_b: (nt, nch, na, nd). Removes the global complex factor per slot."""
    a = g_a / g_a[..., :1, :]
    b = g_b / g_b[..., :1, :]
    return float(np.nanmax(np.abs(a - b)))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sols_a")
    p.add_argument("sols_b")
    p.add_argument("--tol", type=float, default=0.05,
                   help="warn if the invariant gain difference exceeds this (default 0.05)")
    p.add_argument("--tol-pos", type=float, default=1e-3,
                   help="max direction position difference in radians (default 1e-3)")
    args = p.parse_args(argv)

    Ga, ra_a, dec_a = load(args.sols_a)
    Gb, ra_b, dec_b = load(args.sols_b)
    ga, gb = scalar_gains(Ga), scalar_gains(Gb)
    if ga.shape != gb.shape:
        print("incomparable: Sols.G shapes differ (%r vs %r)" % (ga.shape, gb.shape))
        return 1

    ia = direction_order(ra_a, dec_a)
    ib = direction_order(ra_b, dec_b)
    sep = float(np.max(np.maximum(np.abs(ra_a[ia] - ra_b[ib]), np.abs(dec_a[ia] - dec_b[ib]))))
    print("max direction position difference = %.3e rad" % sep)
    if sep > args.tol_pos:
        print("incomparable: facet directions differ (the clustering differs between runs)")
        return 1

    d = phase_invariant_diff(ga[..., ia], gb[..., ib])
    verdict = "ok" if d <= args.tol else "DIFFERS"
    print("max phase-invariant gain diff = %.4f (tol %.3f) %s" % (d, args.tol, verdict))
    return 0 if verdict == "ok" else 1


if __name__ == "__main__":
    sys.exit(main())
