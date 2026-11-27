#!/bin/bash
# Repackage a Measurement-Set archive without the stale casacore table locks.
#
# The archived datasets (Sim.tgz, and the LoTSS "big" folder share) carry
# `table.lock` files left behind by the process that generated them. casacore
# treats those locks inconsistently on Lustre: one run "steals" them, the next
# deadlocks with "Resource deadlock avoided". A clean archive has none.
#
# This only removes `table.lock`. The `table.f0/f1/f0i/f1i` files are the
# StandardStMan/TiledStMan *data*, not locks, and are left untouched.
#
# Usage:
#   repackage_dataset.sh <archive.tgz> [dest.tgz]
#
# Defaults: reads ./Sim.tgz, writes ./Sim_clean.tgz.
set -euo pipefail

SRC="${1:-Sim.tgz}"
DST="${2:-${SRC%.tgz}_clean.tgz}"

die() { echo "error: $*" >&2; exit 1; }

[ -f "$SRC" ] || die "archive not found: $SRC"
[ "$DST" != "$SRC" ] || die "destination must differ from source"

work="$(mktemp -d "${TMPDIR:-/tmp}/repackage.XXXXXX")"
trap 'rm -rf "$work"' EXIT

echo "extract $SRC"
tar -xzf "$SRC" -C "$work"

before=$(find "$work" -name table.lock | wc -l)
echo "stale table.lock files: $before"
find "$work" -name table.lock -delete
after=$(find "$work" -name table.lock | wc -l)
[ "$after" -eq 0 ] || die "could not remove all table.lock files"

# One entry is expected at the archive root (e.g. Sim/, or a bare MS dir).
top="$(find "$work" -mindepth 1 -maxdepth 1 -print -quit)"
echo "re-tar $DST (from $(basename "$top"))"
tar --format=posix --pax-option=delete=atime,delete=ctime \
    --sort=name --numeric-owner --owner=0 --group=0 --mtime='@0' \
    -cf - -C "$work" "$(basename "$top")" | gzip -n > "$DST"

echo "removed=$before remaining=$after"
echo "sha256 $(sha256sum "$DST")"
