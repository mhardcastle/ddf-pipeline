#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from past.utils import old_div
import sys
import pyrap.tables as pt
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from astropy.io import fits
from auxcodes import getpos,getposim
csize=2.5
import subprocess

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
            if 'ddfp' in jobname:
                jobs[jobname[5:]]=status
    return jobs

def plotcircle(name,ra,dec,color,pcolor='black'):
    circle1=plt.Circle((ra,dec),csize,color=color,alpha=0.2)
    plt.gcf().gca().add_artist(circle1)
    plt.scatter(ra,dec,color=pcolor)
    plt.text(ra,dec,name)

jobs=qstat()
circles=[]
s_colours={'unpacking':'red','downloaded':'orange','started':'yellow','complete':'green'}
s_files=['dirin_SSD_init','dirin_SSD','phase1','ampphase1','full_low','full_low_m','full_ampphase1','full_ampphase1m']
for d in sys.argv[1:]:
    os.chdir(d)
    d_colour = 'black'
    shortname=d.split('/')[-1]
    mss=glob.glob('*.ms')
    name=None
    if len(mss)>0:
        status=''
        secondary_status=''
        job_status=''
        name,ra,dec=getpos(mss[1])
        if os.path.isfile('big-mslist.txt'):
            status='downloaded'
        else:
            status='unpacking'
        if os.path.isfile('image_dirin_SSD_init.tessel.reg'):
            status='started'
        if os.path.isfile('summary.txt'):
            status='complete'
        if status == 'started' or status == 'downloaded':
            job_status = ''
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
            
            for ft in s_files:
                if os.path.isfile('image_'+ft+'.app.restored.fits'):
                    secondary_status = ft
        circles.append((name,ra,dec,s_colours[status],d_colour))
        print("%-45s %-15s %8.3f %8.3f %-12s %-10s %s" % (d,name,ra,dec,status,job_status,secondary_status))
        continue
    # else check for image only
    ims=glob.glob('image_full_ampphase1m_shift.int.facetRestored.fits')
    if len(ims)>0:
        ra,dec=getposim(ims[0])
        if name is None:
            name=d.split('/')[-1]
            name=name.split('_')[0]
        print("%-45s %-15s %8.3f %8.3f %-12s" % (d,name,ra,dec,'final'))
        circles.append((name,ra,dec,'blue'))

ras=np.array([c[1] for c in circles])
decs=np.array([c[2] for c in circles])
rarange=max(ras)-min(ras)
decrange=max(decs)-min(decs)
yfigsize=1
xfigsize=(old_div(rarange,decrange))*np.cos(np.mean(decs)*np.pi/180.0)
maxs=max((xfigsize,yfigsize))
xfigsize*=old_div(18,maxs)
yfigsize*=old_div(18,maxs)
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
