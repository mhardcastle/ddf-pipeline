#!/usr/bin/env python3
"""Smoke check for a pipeline run directory.

Usage:
  smoke_check.py RUN_DIR [--reference tests/data/smoke_reference.json]

The pipeline is not pixel-reproducible (threading/FFT/iterative solving), but its
restored-image statistics are stable run-to-run (~1% mean/std, ~5% tails). This check
asserts the named final artefacts exist and that the image's global statistics match the
stored reference within per-statistic tolerances.
"""
import argparse
import json
import os
import sys

import numpy as np
from astropy.io import fits

DEFAULT_REF = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "tests", "data", "smoke_reference.json")


def image_stats(path):
    with fits.open(path) as hdus:
        data = np.asarray(hdus[0].data, dtype=np.float64)
    finite = data[np.isfinite(data)]
    return {
        "mean": float(finite.mean()),
        "std": float(finite.std()),
        "p1": float(np.percentile(finite, 1)),
        "p99": float(np.percentile(finite, 99)),
        "max": float(finite.max()),
    }


def scan_tracebacks(run_dir):
    hits = []
    logs = os.path.join(run_dir, "logs")
    for root, _dirs, files in os.walk(logs):
        for name in files:
            path = os.path.join(root, name)
            try:
                with open(path, errors="ignore") as fh:
                    if "Traceback" in fh.read():
                        hits.append(os.path.relpath(path, run_dir))
            except OSError:
                pass
    return hits


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir")
    parser.add_argument("--reference", default=DEFAULT_REF)
    args = parser.parse_args(argv[1:])

    with open(args.reference) as fh:
        ref = json.load(fh)

    ok = True
    for name in ref["artifacts"]:
        exists = os.path.isfile(os.path.join(args.run_dir, name))
        print(f"[{'ok' if exists else 'MISSING'}] artefact {name}")
        ok = ok and exists

    tracebacks = scan_tracebacks(args.run_dir)
    if tracebacks:
        for hit in tracebacks:
            print(f"[FAIL] Traceback in {hit}")
        ok = False
    else:
        print("[ok] no Traceback in logs")

    image = os.path.join(args.run_dir, ref["image"])
    if not os.path.isfile(image):
        print(f"[MISSING] image {ref['image']}")
        return 1

    stats = image_stats(image)
    for key, value in stats.items():
        expected = ref["stats"][key]
        tol = ref["tolerance"][key]
        rel = abs(value - expected) / (abs(expected) or 1.0)
        passed = rel <= tol
        print(f"[{'ok' if passed else 'FAIL'}] {key}: {value:.5g} vs ref {expected:.5g} "
              f"(rel {rel:.2%}, tol {tol:.0%})")
        ok = ok and passed

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
