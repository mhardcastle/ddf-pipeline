#!/usr/bin/env python3
"""Manifest and compare the files of two pipeline run directories.

Usage:
  compare_runs.py manifest DIR OUT        # write sorted '<sha256>  <relpath>' lines
  compare_runs.py compare A.manifest B.manifest

The manifest hashes every file under DIR (byte SHA256). Note that byte hashes
differ for files carrying timestamps (logs, FITS headers); the compare step tells
you *which* files differ, then you decide the oracle from the stable ones.
"""
import hashlib
import os
import sys


def manifest(directory, out):
    rows = []
    for root, _dirs, files in os.walk(directory):
        for name in files:
            path = os.path.join(root, name)
            rel = os.path.relpath(path, directory)
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            rows.append((rel, h.hexdigest()))
    rows.sort(key=lambda r: r[0])
    with open(out, "w") as fh:
        fh.write("\n".join(f"{digest}  {rel}" for rel, digest in rows) + "\n")
    print(f"{len(rows)} files -> {out}")


def load(path):
    table = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            digest, rel = line.split("  ", 1)
            table[rel] = digest
    return table


def compare(a, b):
    left, right = load(a), load(b)
    same = [p for p in left if p in right and left[p] == right[p]]
    changed = [p for p in left if p in right and left[p] != right[p]]
    only_left = [p for p in left if p not in right]
    only_right = [p for p in right if p not in left]

    print(f"identical : {len(same)}")
    print(f"changed   : {len(changed)}")
    print(f"only in A : {len(only_left)}")
    print(f"only in B : {len(only_right)}")
    for label, items in (("CHANGED", changed), ("ONLY_A", only_left), ("ONLY_B", only_right)):
        if items:
            print(f"\n{label}:")
            for p in sorted(items):
                print(f"  {p}")


def main(argv):
    if len(argv) == 4 and argv[1] == "manifest":
        manifest(argv[2], argv[3])
    elif len(argv) == 4 and argv[1] == "compare":
        compare(argv[2], argv[3])
    else:
        print(__doc__)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
