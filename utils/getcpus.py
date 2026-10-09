import os
import psutil


def get_physical_cpus():
    return psutil.cpu_count(logical=False)


def _affinity_count():
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return None


def getcpus():
    aff = _affinity_count()
    nodefile = os.getenv('PBS_NODEFILE')
    if nodefile and os.path.isfile(nodefile):
        return len(open(nodefile).readlines())
    slurmcpus = os.getenv('SLURM_JOB_CPUS_PER_NODE')
    if slurmcpus:
        n = int(slurmcpus)
        if aff and aff < n:
            n = aff
        return n
    return aff or get_physical_cpus()
