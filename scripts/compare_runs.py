#!/usr/bin/env python3
"""Manifest and compare the files of two pipeline run directories.

Usage:
  compare_runs.py manifest DIR OUT [--exclude PAT,PAT]
  compare_runs.py compare A.manifest B.manifest

manifest hashes every file under DIR (byte SHA256) into a sorted manifest, skipping
files whose relative path matches any --exclude glob. Use --exclude to keep noisy
monitoring artefacts (mem, iomonitor*) out of the science comparison.
"""
import argparse
import fnmatch
import hashlib
import os
import sys


def excluded(rel, patterns):
    base = os.path.basename(rel)
    return any(fnmatch.fnmatch(rel, p) or fnmatch.fnmatch(base, p) for p in patterns)


def manifest(directory, out, excludes):
    rows = []
    skipped = 0
    for root, _dirs, files in os.walk(directory):
        for name in files:
            path = os.path.join(root, name)
            rel = os.path.relpath(path, directory)
            if excluded(rel, excludes):
                skipped += 1
                continue
            digest = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    digest.update(chunk)
            rows.append((rel, digest.hexdigest()))
    rows.sort(key=lambda r: r[0])
    with open(out, "w") as fh:
        fh.write("\n".join(f"{d}  {r}" for r, d in rows) + "\n")
    print(f"{len(rows)} files -> {out} ({skipped} excluded)")


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
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("manifest", help="hash a run directory into a manifest")
    m.add_argument("directory")
    m.add_argument("out")
    m.add_argument("--exclude", default="",
                   help="comma-separated globs to skip (e.g. 'mem,iomonitor*,summary.txt')")

    c = sub.add_parser("compare", help="compare two manifests")
    c.add_argument("a")
    c.add_argument("b")

    args = parser.parse_args(argv[1:])
    if args.cmd == "manifest":
        patterns = [x.strip() for x in args.exclude.split(",") if x.strip()]
        manifest(args.directory, args.out, patterns)
    else:
        compare(args.a, args.b)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
