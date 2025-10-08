import os
import psutil

def get_physical_cpus():
    return len(psutil.Process().cpu_affinity())

def getcpus():
    nodefile=os.getenv('PBS_NODEFILE')
    slurmcpus=os.getenv('SLURM_NTASKS_PER_NODE')
    if nodefile and os.path.isfile(nodefile):
        lines=len(open(nodefile).readlines())
        return lines
    elif slurmcpus:
        # if multiple node: 256(x2)
        return int(slurmcpus.split("(")[0])
    else:
        return get_physical_cpus()
