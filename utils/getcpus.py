import os
import re
import psutil


def get_physical_cpus():
    return len(psutil.Process().cpu_affinity())


def _affinity_count():
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return None


def _parse_slurm_cpus(s):
    """Parse SLURM_JOB_CPUS_PER_NODE, e.g. '5(x2)' or '72(x2),36'.
    Returns a list with one CPU count per node."""
    counts = []
    for part in s.split(','):
        m = re.fullmatch(r'\s*(\d+)(?:\(x(\d+)\))?\s*', part)
        if m:
            counts.extend([int(m.group(1))] * int(m.group(2) or 1))
    return counts


def getcpus():
    aff = _affinity_count()

    nodefile = os.getenv('PBS_NODEFILE')
    if nodefile and os.path.isfile(nodefile):
        with open(nodefile) as f:
            return sum(1 for line in f if line.strip())

    # CPUs allocated on this node
    n = 0
    oncpus = os.getenv('SLURM_CPUS_ON_NODE')
    if oncpus and oncpus.isdigit():
        n = int(oncpus)
    else:
        slurmcpus = os.getenv('SLURM_JOB_CPUS_PER_NODE')
        if slurmcpus:
            counts = _parse_slurm_cpus(slurmcpus)
            if counts:
                n = counts[0]

    if n:
        if aff and aff < n:
            n = aff
        return n

    return aff or get_physical_cpus()
