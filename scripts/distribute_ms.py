# Distribute measurement sets across hosts
# Inputs:
# -- mpirun-style host list
# -- location on these hosts to work in
# Assumes we are in a working directory where make_mslists.py has already run

import sys
import os
import glob

hostlist=sys.argv[1]
location=sys.argv[2]

hosts=[l.split()[0] for l in open(hostlist).readlines()]
for host in hosts:
    os.system(f'ssh {host} mkdir {location}')

# we rely on make_mslists.py to sort out the mslist.txt and
# big-mslist.txt.  so then the algorithm is to split first the
# mslist.txt and then the big-mslist.txt equally between the nodes --
# taking the mslist.txt first avoids imbalances in the early pipeline runs.

mslist=[l.rstrip() for l in open('mslist.txt').readlines()]
hostd={}
counter=0
outfile=open('dist_mslist.txt','w')
for ms in mslist:
    host=hosts[counter % len(hosts)]
    outfile.write(f'{host}:{location}/{ms}\n')
    os.system(f'rsync -av {ms} {host}:{location}')
    hostd[ms]=host
    counter+=1
outfile.close()

bigmslist=[l.rstrip() for l in open('big-mslist.txt').readlines()]
outfile=open('dist_big-mslist.txt','w')
for ms in bigmslist:
    if ms in hostd:
        host=hostd[ms]
    else:
        host=hosts[counter % len(hosts)]
        counter+=1
    outfile.write(f'{host}:{location}/{ms}\n')
    os.system(f'rsync -av {ms} {host}:{location}')
    hostd[ms]=host
outfile.close()

    
    
