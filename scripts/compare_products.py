#!/usr/bin/env python3
"""Data-level comparison of two pipeline products (.fits, .npy, .npz).

Usage:
  compare_products.py A B [--rtol 1e-6] [--atol 0.0]
  compare_products.py stats A

Compares arrays instead of bytes, so floating-point/threading differences show up as a
magnitude rather than a false mismatch. Handles structured and complex arrays (calibration
solutions). For FITS it compares each HDU's data, reports the data scale (ref_absmax) and
the differing header keys (ignoring volatile dates/checksums).
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

    if a.dtype.names:  # structured
        if a.dtype.names != b.dtype.names:
            return f"fields differ {a.dtype.names} vs {b.dtype.names}"
        parts = [f"{f}: {compare_arrays(a[f], b[f], rtol, atol)}" for f in a.dtype.names]
        return "; ".join(parts)

    kind = a.dtype.kind
    if kind in ("U", "S", "O"):
        return f"{'equal' if np.array_equal(a, b) else 'DIFFERS'} ({a.size} items)"
    if a.size == 0:
        return "empty"

    if kind == "c":  # complex
        diff = np.abs(a - b)
        rms_a = float(np.sqrt(np.mean(np.abs(a) ** 2))) or 1.0
        rms_diff = float(np.sqrt(np.mean(diff ** 2)))
        ok = np.allclose(a, b, rtol=rtol, atol=atol)
        return (f"{'equal' if ok else 'DIFFERS'} max_abs={float(diff.max()):.3e} "
                f"rms_ratio={rms_diff / rms_a:.3e} ref_absmax={float(np.abs(a).max()):.3e}")

    af = a.astype(np.float64, copy=False)
    bf = b.astype(np.float64, copy=False)
    finite_a, finite_b = np.isfinite(af), np.isfinite(bf)
    if not np.array_equal(finite_a, finite_b):
        return f"non-finite mask differs (A={finite_a.sum()} finite, B={finite_b.sum()})"
    if not finite_a.any():
        return "identical (all non-finite)"
    fa, fb = af[finite_a], bf[finite_b]
    diff = np.abs(fa - fb)
    ref = float(np.abs(fa).max()) or 1.0
    rms_a = float(np.sqrt(np.mean(fa * fa))) or 1.0
    rms_diff = float(np.sqrt(np.mean(diff * diff)))
    max_abs = float(diff.max())
    max_rel = float((diff / np.maximum(np.abs(fa), 1e-30)).max())
    ok = np.allclose(fa, fb, rtol=rtol, atol=atol)
    return (f"{'equal' if ok else 'DIFFERS'} max_abs={max_abs:.3e} max_rel={max_rel:.3e} "
            f"rms_diff={rms_diff:.3e} rms_ratio={rms_diff / rms_a:.3e} ref_absmax={ref:.3e}")


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
        except Exception as exc:
            print(f"  {k}: not comparable ({exc})")


def stats(path):
    from astropy.io import fits
    with fits.open(path) as hdus:
        for i, hdu in enumerate(hdus):
            if hdu.data is None:
                continue
            d = np.asarray(hdu.data)
            if d.dtype.kind == "c":
                d = np.abs(d)
            if d.dtype.kind not in "biuf":
                continue
            d = d.astype(np.float64, copy=False)
            finite = d[np.isfinite(d)]
            if not finite.size:
                print(f"HDU {i}: all non-finite")
                continue
            print(f"HDU {i}: shape={d.shape} min={finite.min():.4g} max={finite.max():.4g} "
                  f"mean={finite.mean():.4g} std={finite.std():.4g} median={np.median(finite):.4g} "
                  f"p1={np.percentile(finite, 1):.4g} p99={np.percentile(finite, 99):.4g} "
                  f"nonfinite={d.size - finite.size}")


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("a")
    parser.add_argument("b", nargs="?")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=0.0)
    args = parser.parse_args(argv[1:])

    if args.a == "stats":
        if not args.b:
            parser.error("stats needs a file")
        stats(args.b)
        return 0

    if not args.b:
        parser.error("compare needs two files")
    base = args.a[:-3] if args.a.endswith(".gz") else args.a
    ext = os.path.splitext(base)[1].lower()
    print(f"{args.a}\n  vs {args.b}")
    if ext in (".fits", ".fit"):
        compare_fits(args.a, args.b, args.rtol, args.atol)
    elif ext == ".npz":
        compare_npz(args.a, args.b, args.rtol, args.atol)
    elif ext == ".npy":
        print(f"  data: {compare_arrays(np.load(args.a), np.load(args.b), args.rtol, args.atol)}")
    else:
        print("  not a fits/npy/npz file")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
