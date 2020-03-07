#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from past.utils import old_div
import numpy as np
import matplotlib.pyplot as plt
from surveys_db import SurveysDB
import subprocess
csize=2.5

def qstat():
    # check the torque queue for lofar jobs
    # necessarily only works on Herts systems
    try:
        results=subprocess.check_output(['qstat', '-a','lofar']).split('\n')
    except OSError:
        results=[]
    jobs={}
    for l in results:
        bits=l.split()
        if len(bits)>10:
            jobname=bits[3]
            status=bits[9]
            if 'ddf-' in jobname:
                jobs[jobname[4:]]=status
            if 'ddfp-' in jobname:
                jobs[jobname[5:]]=status
    return jobs

def plotcircle(name,ra,dec,color,pcolor='black'):
    circle1=plt.Circle((ra,dec),csize,color=color,alpha=0.2)
    plt.gcf().gca().add_artist(circle1)
    plt.scatter(ra,dec,color=pcolor)
    plt.text(ra,dec,name)

jobs=qstat()
print(jobs)
circles=[]

s_colours={'D/L failed':'black','List failed':'black','Downloading':'red','Downloaded':'orange','Unpacking':'orange','Unpacked':'orange','Ready':'yellow','Queued':'blue','Running':'cyan','Complete':'green','Archiving':'green','Archived':'olive','Stopped':'magenta','Failed':'magenta','DI_started':'white'}
sdb=SurveysDB(readonly=True)
sdb.cur.execute('select * from fields where status!="Not started"')
fields=sdb.cur.fetchall()
sdb.close()

for f in fields:
    status=f['status']
    shortname=f['id']
    job_status = ''
    d_colour='black'
    if status=='Running' or status=='Ready' or status=='Failed':
        if len(jobs)>0:
            if shortname not in jobs:
                job_status = 'not queued'
                d_colour = 'red'
            elif jobs[shortname] == 'Q':
                job_status = 'queued'
                d_colour = 'orange'
            elif jobs[shortname] == 'R':
                job_status = 'running'
                d_colour = 'green'
            else:
                job_status = 'other ('+jobs[shortname]+')'
            
    circles.append((shortname,f['ra'],f['decl'],s_colours[status],d_colour))
    print("%-15s %8.3f %8.3f %-12s %-20s %-10s" % (shortname,f['ra'],f['decl'],status,f['nodename'],job_status))


ras=np.array([c[1] for c in circles])
decs=np.array([c[2] for c in circles])
rarange=max(ras)-min(ras)
decrange=max(decs)-min(decs)
yfigsize=1
xfigsize=(old_div(rarange,decrange))*np.cos(np.mean(decs)*np.pi/180.0)
maxs=max((xfigsize,yfigsize))
xfigsize*=old_div(16,maxs)
yfigsize*=old_div(12,maxs)
print(xfigsize,yfigsize)

plt.figure(figsize=(xfigsize,yfigsize))

plt.xlim(np.min(ras)-csize,np.max(ras)+csize)
plt.ylim(np.min(decs)-csize,np.max(decs)+csize)


for c in circles:
    plotcircle(*c)

plt.gca().invert_xaxis()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.tight_layout()
plt.show()
