import os
import psutil

def get_physical_cpus():
    return psutil.cpu_count()
    # find total number of physical cores, ignoring hyperthreading
    lines=open('/proc/cpuinfo').readlines()
    cpus=[]
    for l in lines:
        bits=l.split(':')
        if 'physical id' in l:
            phys=int(bits[1])
        if 'core id' in l:
            core=int(bits[1])
        if l=='\n':
            cpus.append((phys,core))
    cpus=set(cpus)
    return len(cpus)

def getcpus():
    nodefile=os.getenv('PBS_NODEFILE')
    if nodefile:
        lines=len(open(nodefile).readlines())
        return lines
    else:
        return get_physical_cpus()
