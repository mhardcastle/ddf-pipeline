#!/usr/bin/env python3
"""SLURM performance report: time, CPU efficiency, memory and energy per job.

Reads the pipeline step (``<jobid>.0``) from ``sacct`` and prints one line per job:

    job  nodes  elapsed_s  cpu_s  cpu_eff  maxrss_GiB  maxvmsize_GiB  energy_J

Usage:
  perf_report.py -j 839738,841042
  NODES="12" ... perf_report.py --from-csv "$SCRATCH/ddf-sim-dataset/scaling/scaling.csv"
"""
import argparse
import subprocess
import sys

FIELDS = ["JobID", "State", "Elapsed", "TotalCPU", "AllocCPUS", "NNodes",
          "MaxRSS", "MaxVMSize", "ConsumedEnergyRaw"]


def hms_to_s(value):
    """'00:14:13' or '1-02:03:04' -> seconds; 0 for Unknown/empty."""
    value = (value or "").strip()
    if not value or value == "Unknown":
        return 0.0
    days = 0
    if "-" in value:
        d, value = value.split("-", 1)
        days = int(d)
    parts = [float(x) for x in value.split(":")]
    if len(parts) == 3:
        h, m, sec = parts
    elif len(parts) == 2:
        h, m, sec = 0.0, parts[0], parts[1]
    else:
        h, m, sec = 0.0, 0.0, parts[0]
    return days * 86400 + h * 3600 + m * 60 + sec


def size_to_bytes(value):
    """'14628K', '1.5G', '512M', '0' -> bytes."""
    value = (value or "").strip()
    if not value or value in ("0", "Unknown"):
        return 0
    units = {"K": 1024, "M": 1024 ** 2, "G": 1024 ** 3, "T": 1024 ** 4}
    if value[-1] in units:
        return int(float(value[:-1]) * units[value[-1]])
    return int(float(value))


def sacct_rows(jobids):
    cmd = ["sacct", "-P", "-n", "-j", str(jobids), "--format=" + ",".join(FIELDS)]
    out = subprocess.check_output(cmd, text=True)
    return [line.split("|") for line in out.splitlines() if line.strip()]


def metrics(row):
    jid, state, elapsed, totalcpu, cpus, nodes, maxrss, maxvmsize, energy = row
    elapsed_s = hms_to_s(elapsed)
    cpu_s = hms_to_s(totalcpu)
    alloc = int(cpus) if cpus.strip().isdigit() else 0
    nnodes = int(nodes) if nodes.strip().isdigit() else 0
    eff = cpu_s / (elapsed_s * alloc) if elapsed_s and alloc else 0.0
    return {
        "job": jid, "state": state, "nodes": nnodes,
        "elapsed_s": elapsed_s, "cpu_s": cpu_s, "cpu_eff": eff,
        "maxrss_GiB": size_to_bytes(maxrss) / 1024 ** 3,
        "maxvmsize_GiB": size_to_bytes(maxvmsize) / 1024 ** 3,
        "energy_J": int(energy) if energy.strip().isdigit() else 0,
    }


def jobids_from_csv(path):
    ids = []
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.strip().split(",")
            if len(parts) >= 3 and parts[2].isdigit():
                ids.append(parts[2])
    return ids


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-j", "--jobs", default=None, help="comma-separated job IDs")
    p.add_argument("--from-csv", default=None, help="read job IDs (column 3) from a scaling.csv")
    args = p.parse_args(argv)

    ids = args.jobs
    if args.from_csv:
        ids = ",".join(jobids_from_csv(args.from_csv))
    if not ids:
        p.error("give -j or --from-csv")

    rows = sacct_rows(ids)
    print("%-12s %-10s %5s %9s %9s %7s %11s %13s %10s"
          % ("job", "state", "nodes", "elapsed_s", "cpu_s", "cpu_eff",
             "maxrss_GiB", "maxvmsize_GiB", "energy_J"))
    for row in rows:
        if not row[0].endswith(".0"):
            continue  # the pipeline step, not .batch/.extern
        m = metrics(row)
        print("%-12s %-10s %5d %9.0f %9.0f %6.1f%% %11.2f %13.2f %10d"
              % (m["job"], m["state"], m["nodes"], m["elapsed_s"], m["cpu_s"],
                 100 * m["cpu_eff"], m["maxrss_GiB"], m["maxvmsize_GiB"], m["energy_J"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
