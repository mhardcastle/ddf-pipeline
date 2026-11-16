#!/usr/bin/env python
"""
Split each Measurement Set listed in a text file into nMS time chunks.

For every MS in the list, creates <MSName>.tChunk_<i>.ms next to the original,
for i = 0 .. nMS-1, and writes the list of all chunks to <mslist>.chunked.txt
(e.g. mslist.txt -> mslist.chunked.txt). By default the chunks are casacore *reference tables*
(virtual: no data copied). Use --deep to write standalone copies instead.

Reference tables have no subtable directories of their own: the subtables are
reached through table keywords, e.g. pt.table("chunk.ms::SPECTRAL_WINDOW").
For tools that open "<ms>/SPECTRAL_WINDOW" as a path instead, --link-subtables
creates a relative symlink to each subtable of the original MS inside every
chunk (existing chunks also get the links on a re-run).

Fast path (time-ordered MS with the same number of rows per timestamp, which is
the usual layout): only a handful of TIME values are read, and chunks are made
with a row range (OFFSET/LIMIT), so nothing scans the full table.
Fallback (irregular layout): chunks are selected with a TIME range computed
from the first/last time and the integration interval.

Usage:
    python chunk_ms_time.py mslist.txt --nMS 4
    python chunk_ms_time.py mslist.txt --nMS 4 --overwrite
    python chunk_ms_time.py mslist.txt --nMS 4 --deep
    python chunk_ms_time.py mslist.txt --nMS 4 --dry-run
    python chunk_ms_time.py mslist.txt --nMS 4 --link-subtables
"""

import argparse
import logging
import os
import sys
import time

import numpy as np
import casacore.tables as pt

log = logging.getLogger("chunk_ms")


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
class Timer:
    """Context manager logging how long a step took."""

    def __init__(self, what):
        self.what = what

    def __enter__(self):
        log.info("  %s ...", self.what)
        self.t0 = time.time()
        return self

    def __exit__(self, *exc):
        log.info("  %s done in %.1f s", self.what, time.time() - self.t0)


def read_mslist(path):
    """Return MS paths from a text file, skipping blank lines and # comments."""
    mss = []
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                mss.append(os.path.abspath(line.rstrip("/")))
    return mss


def chunk_name(ms, i):
    """/path/obs.ms -> /path/obs.tChunk_<i>.ms"""
    base = ms[:-3] if ms.lower().endswith(".ms") else ms
    return f"{base}.tChunk_{i}.ms"


def ffmt(x):
    """Plain float literal for TaQL (numpy 2 repr gives 'np.float64(...)')."""
    return f"{float(x):.6f}"


def rows_per_timestamp(t, nrows, probe=200000):
    """
    Count how many rows share the first timestamp, reading TIME only in
    growing blocks from the start of the table (never the whole column).
    """
    t0 = t.getcell("TIME", 0)
    n = min(probe, nrows)
    while True:
        times = t.getcol("TIME", 0, n)
        diff = np.nonzero(times != t0)[0]
        if diff.size:
            return int(diff[0])
        if n == nrows:
            return nrows  # single timestamp
        n = min(2 * n, nrows)


def check_regular_layout(t, nrows, nbl):
    """
    Verify the MS is time ordered with nbl rows per timestamp by sampling:
    the last block, and a few block boundaries. Cheap: only getcell calls.
    """
    if nrows % nbl:
        return False, f"nrows ({nrows}) not a multiple of rows/timestamp ({nbl})"
    ntimes = nrows // nbl
    last = t.getcol("TIME", nrows - nbl, nbl)
    if not np.all(last == last[0]):
        return False, "last block of rows does not share one timestamp"
    for k in np.unique(np.linspace(1, ntimes - 1, min(ntimes - 1, 20)).astype(int)):
        r = k * nbl
        if not (t.getcell("TIME", r - 1) < t.getcell("TIME", r) == t.getcell("TIME", r + nbl - 1)):
            return False, f"irregular block around row {r}"
    return True, ""


def link_subtables(out, subtables):
    """
    Make <out>/<SUBTABLE> a relative symlink to each subtable of the original
    MS, unless something with that name already exists (e.g. a deep copy).
    """
    made = []
    for sub in subtables:
        dst = os.path.join(out, os.path.basename(sub.rstrip("/")))
        if os.path.lexists(dst):
            continue
        os.symlink(os.path.relpath(sub, out), dst)
        made.append(os.path.basename(dst))
    if made:
        log.info("    linked subtables: %s", ", ".join(made))


# ----------------------------------------------------------------------------
# main work
# ----------------------------------------------------------------------------
def write_chunk(t, out, deep, offset=None, limit=None, cond=None):
    """Create one chunk, either by row range or by TaQL condition."""
    if cond is None:
        sel = t.query(offset=offset, limit=limit, name="" if deep else out)
    else:
        sel = t.query(cond, name="" if deep else out)
    nrows = sel.nrows()
    if deep:
        with Timer(f"deep copy of {nrows} rows to {os.path.basename(out)}"):
            sel.copy(out, deep=True)
    sel.close()
    return nrows


def chunk_ms(ms, nms, overwrite=False, deep=False, dry_run=False, link=False):
    log.info("=" * 70)
    log.info("MS: %s", ms)
    if not os.path.isdir(ms):
        log.warning("  not found, skipping")
        return None

    with Timer("opening table"):
        t = pt.table(ms, ack=False)

    try:
        subtables = [os.path.abspath(x) for x in t.getsubtables()]
        log.info("  subtables: %s", ", ".join(os.path.basename(x) for x in subtables))
        nrows = t.nrows()
        tfirst = t.getcell("TIME", 0)
        tlast = t.getcell("TIME", nrows - 1)
        interval = t.getcell("INTERVAL", 0)
        n_est = int(round((tlast - tfirst) / interval)) + 1
        log.info("  rows: %d", nrows)
        log.info("  first/last TIME: %s / %s  (span %.1f s)", ffmt(tfirst), ffmt(tlast), tlast - tfirst)
        log.info("  INTERVAL: %.3f s -> ~%d timestamps expected", interval, n_est)

        with Timer("finding rows per timestamp"):
            nbl = rows_per_timestamp(t, nrows)
        log.info("  rows per timestamp: %d", nbl)

        with Timer("checking time-ordered regular layout"):
            regular, why = check_regular_layout(t, nrows, nbl)

        if regular:
            ntimes = nrows // nbl
            log.info("  layout regular: %d timestamps -> using row ranges (fast)", ntimes)
            if ntimes != n_est:
                log.warning("  %d timestamps found vs %d from INTERVAL: gaps in the data?",
                            ntimes, n_est)
        else:
            ntimes = n_est
            log.warning("  layout irregular (%s) -> falling back to TIME selection", why)

        if nms > ntimes:
            log.warning("  nMS=%d > %d timestamps, skipping", nms, ntimes)
            return None

        groups = np.array_split(np.arange(ntimes), nms)
        chunks = []
        for i, g in enumerate(groups):
            out = chunk_name(ms, i)
            name = os.path.basename(out)

            if regular:
                offset, limit = int(g[0]) * nbl, len(g) * nbl
                cond = None
                desc = f"rows {offset}..{offset + limit - 1} ({len(g)} timestamps)"
            else:
                # cut halfway between expected timestamps to avoid float equality
                lo = tfirst + (g[0] - 0.5) * interval
                hi = tfirst + (g[-1] + 0.5) * interval
                parts = []
                if i > 0:
                    parts.append(f"TIME > {ffmt(lo)}")
                if i < nms - 1:
                    parts.append(f"TIME < {ffmt(hi)}")
                cond = " && ".join(parts) or "True"
                offset = limit = None
                desc = f"WHERE {cond}"

            log.info("  chunk %d/%d -> %s: %s", i + 1, nms, name, desc)

            if os.path.exists(out):
                if not overwrite:
                    log.info("    exists, skipping (use --overwrite)")
                    if link and not dry_run:
                        link_subtables(out, subtables)
                    chunks.append(out)
                    continue
                if not dry_run:
                    log.info("    removing existing %s", name)
                    pt.tabledelete(out)

            if dry_run:
                log.info("    [dry-run] not written")
                chunks.append(out)
                continue

            t0 = time.time()
            n = write_chunk(t, out, deep, offset, limit, cond)
            if link:
                link_subtables(out, subtables)
            log.info("    written: %d rows%s in %.1f s",
                     n, " (deep copy)" if deep else " (reference table)", time.time() - t0)
            chunks.append(out)
    finally:
        t.close()
    return chunks


def main():
    p = argparse.ArgumentParser(description="Split MSs into nMS time chunks.")
    p.add_argument("mslist", help="Text file with one MS path per line")
    p.add_argument("--nMS", type=int, required=True, help="Number of time chunks per MS")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing chunks")
    p.add_argument("--deep", action="store_true",
                   help="Write standalone copies instead of reference tables")
    p.add_argument("--dry-run", action="store_true", help="Only print what would be done")
    p.add_argument("--link-subtables", action="store_true",
                   help="Symlink the original subtables into each chunk, for tools "
                        "that open <ms>/SUBTABLE as a path instead of <ms>::SUBTABLE")
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                        format="%(asctime)s %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")

    if args.nMS < 1:
        sys.exit("--nMS must be >= 1")

    mss = read_mslist(args.mslist)
    if not mss:
        sys.exit(f"No MS found in {args.mslist}")
    log.info("%d MS(s) in %s, nMS=%d%s%s%s", len(mss), args.mslist, args.nMS,
             ", deep copy" if args.deep else "",
             ", link subtables" if args.link_subtables else "",
             ", DRY RUN" if args.dry_run else "")

    t0 = time.time()
    n_ok = 0
    all_chunks = []
    for k, ms in enumerate(mss, 1):
        log.info("[%d/%d]", k, len(mss))
        chunks = chunk_ms(ms, args.nMS, args.overwrite, args.deep, args.dry_run,
                          args.link_subtables)
        if chunks is not None:
            n_ok += 1
            all_chunks.extend(chunks)
    log.info("=" * 70)

    root, ext = os.path.splitext(args.mslist)
    out_list = f"{root}.chunked{ext or '.txt'}"
    if args.dry_run:
        log.info("[dry-run] would write %d chunk paths to %s", len(all_chunks), out_list)
    else:
        with open(out_list, "w") as f:
            f.write("\n".join(all_chunks) + ("\n" if all_chunks else ""))
        log.info("Wrote %d chunk paths to %s", len(all_chunks), out_list)

    log.info("Done: %d/%d MSs processed in %.1f s", n_ok, len(mss), time.time() - t0)


if __name__ == "__main__":
    main()

    
