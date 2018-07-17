import os
import psutil

def get_physical_cpus():
    return psutil.cpu_count(logical=False)

def getcpus():
    nodefile=os.getenv('PBS_NODEFILE')
    if nodefile:
        lines=len(open(nodefile).readlines())
        return lines
    else:
        return get_physical_cpus()
