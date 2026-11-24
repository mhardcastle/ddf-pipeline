#!/usr/bin/env python3
"""Data-level comparison of two pipeline products (.fits, .npy, .npz).

Usage:
  compare_products.py A B [--rtol 1e-6] [--atol 0.0]

Compares arrays instead of bytes, so floating-point/threading differences show up as a
magnitude rather than a false mismatch. For FITS it compares each HDU's data and prints
the header keys that differ (ignoring volatile date/checksum keys).
"""
import argparse
import os
import sys

import numpy as np

VOLATILE = {"DATE", "DATE-OBS", "DATE-MAP", "TIMESYS", "ORIGIN", "HISTORY", "COMMENT",
            "CHECKSUM", "DATASUM"}


def compare_arrays(a, b, rtol, atol):
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return f"shape differs {a.shape} vs {b.shape}"
    if a.size == 0:
        return "empty"
    a = a.astype(np.float64, copy=False)
    b = b.astype(np.float64, copy=False)
    finite_a, finite_b = np.isfinite(a), np.isfinite(b)
    if not np.array_equal(finite_a, finite_b):
        return f"non-finite mask differs (A={finite_a.sum()} finite, B={finite_b.sum()})"
    if not finite_a.any():
        return "identical (all non-finite)"
    fa, fb = a[finite_a], b[finite_b]
    diff = np.abs(fa - fb)
    max_abs = float(diff.max())
    max_rel = float((diff / np.maximum(np.abs(fa), 1e-30)).max())
    ok = np.allclose(fa, fb, rtol=rtol, atol=atol)
    return f"{'equal' if ok else 'DIFFERS'} max_abs={max_abs:.3e} max_rel={max_rel:.3e}"


def compare_fits(x, y, rtol, atol):
    from astropy.io import fits
    with fits.open(x) as ha, fits.open(y) as hb:
        if len(ha) != len(hb):
            print(f"  HDU count differs: {len(ha)} vs {len(hb)}")
        for i, (ea, eb) in enumerate(zip(ha, hb)):
            label = f"HDU {i} ({type(ea).__name__})"
            if ea.data is None and eb.data is None:
                print(f"  {label}: no data")
            elif ea.data is None or eb.data is None:
                print(f"  {label}: only one HDU has data")
            else:
                print(f"  {label}: {compare_arrays(ea.data, eb.data, rtol, atol)}")
            ka = {k: v for k, v in ea.header.items() if k not in VOLATILE}
            kb = {k: v for k, v in eb.header.items() if k not in VOLATILE}
            for k in sorted(set(ka) | set(kb)):
                if ka.get(k) != kb.get(k):
                    print(f"    header {k}: {ka.get(k)!r} vs {kb.get(k)!r}")


def compare_npz(x, y, rtol, atol):
    a = np.load(x, allow_pickle=True)
    b = np.load(y, allow_pickle=True)
    fa, fb = set(a.files), set(b.files)
    if fa != fb:
        print(f"  keys differ: only A={sorted(fa - fb)} only B={sorted(fb - fa)}")
    for k in sorted(fa & fb):
        try:
            print(f"  {k}: {compare_arrays(a[k], b[k], rtol, atol)}")
        except Exception as exc:  # non-numeric content
            print(f"  {k}: not numeric ({exc})")


def compare_npy(x, y, rtol, atol):
    print(f"  data: {compare_arrays(np.load(x), np.load(y), rtol, atol)}")


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("a")
    parser.add_argument("b")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=0.0)
    args = parser.parse_args(argv[1:])

    base = args.a[:-3] if args.a.endswith(".gz") else args.a
    ext = os.path.splitext(base)[1].lower()
    print(f"{args.a}\n  vs {args.b}")
    if ext in (".fits", ".fit"):
        compare_fits(args.a, args.b, args.rtol, args.atol)
    elif ext == ".npz":
        compare_npz(args.a, args.b, args.rtol, args.atol)
    elif ext == ".npy":
        compare_npy(args.a, args.b, args.rtol, args.atol)
    else:
        print("  not a fits/npy/npz file")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
