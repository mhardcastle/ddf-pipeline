import os
import psutil

def get_physical_cpus():
    return len(psutil.Process().cpu_affinity())

def getcpus():
    nodefile=os.getenv('PBS_NODEFILE')
    if nodefile and os.path.isfile(nodefile):
        lines=len(open(nodefile).readlines())
        return lines
    else:
        return get_physical_cpus()
