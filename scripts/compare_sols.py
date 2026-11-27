#!/usr/bin/env python3
"""Compare two killMS .sols.npz up to the self-calibration degeneracy.

The gains are defined only up to a global complex factor per solution slot (the
self-cal degeneracy, Smirnov & Tasse 2015), so a raw element-wise diff of Sols.G
overstates the difference. This reports the invariant difference
max |g_a/g_a[0] - g_b/g_b[0]| over all (time, channel) slots and antennas.

If the two runs solved a different number of directions (facets), the shapes
differ and the gains are not comparable; that is reported instead.
"""
import argparse
import sys

import numpy as np


def gains_from_G(G):
    G = np.asarray(G, dtype=complex)
    if G.ndim == 6:
        return G[..., 0, 0][:, :, :, 0]  # (nt, nch, na), direction 0
    if G.ndim == 4:
        return G[..., 0]
    raise ValueError("unexpected Sols.G shape %r" % (G.shape,))


def load(path):
    d = np.load(path, allow_pickle=False)
    return gains_from_G(d["Sols"]["G"])


def phase_invariant_diff(g_a, g_b):
    a = g_a / g_a[..., :1]
    b = g_b / g_b[..., :1]
    return float(np.nanmax(np.abs(a - b)))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sols_a")
    p.add_argument("sols_b")
    p.add_argument("--tol", type=float, default=0.05,
                   help="warn if the invariant gain difference exceeds this (default 0.05)")
    args = p.parse_args(argv)

    g_a, g_b = load(args.sols_a), load(args.sols_b)
    if g_a.shape != g_b.shape:
        print("incomparable: Sols.G shapes differ (%r vs %r) — the facet directions differ"
              % (g_a.shape, g_b.shape))
        return 1
    d = phase_invariant_diff(g_a, g_b)
    verdict = "ok" if d <= args.tol else "DIFFERS"
    print("max phase-invariant gain diff = %.4f (tol %.3f) %s" % (d, args.tol, verdict))
    return 0 if verdict == "ok" else 1


if __name__ == "__main__":
    sys.exit(main())
