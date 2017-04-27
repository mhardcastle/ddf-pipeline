#!/usr/bin/env python

import sys
import pyrap.tables as pt
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from astropy.io import fits
csize=2.5

def getpos(ms):
    t = pt.table(ms+ '/OBSERVATION', readonly=True, ack=False)
    name=t[0]['LOFAR_TARGET']

    t = pt.table(ms+'/FIELD', readonly=True, ack=False)

    direction = t[0]['PHASE_DIR']
    ra, dec = direction[0]

    if (ra<0):
        ra+=2*np.pi;

    return name[0],ra*(180/np.pi),dec*(180/np.pi)

def plotcircle(name,ra,dec,color):
    circle1=plt.Circle((ra,dec),csize,color=color,alpha=0.2)
    plt.gcf().gca().add_artist(circle1)
    plt.scatter(ra,dec)
    plt.text(ra,dec,name)

def getposim(image):
    hdus=fits.open(image)
    ra=hdus[0].header['CRVAL1']
    dec=hdus[0].header['CRVAL2']
    return ra,dec

circles=[]
s_colours={'downloading':'red','downloaded':'orange','started':'yellow','complete':'green'}
for d in sys.argv[1:]:
    os.chdir(d)
    mss=glob.glob('*.ms')
    name=None
    if len(mss)>0:
        name,ra,dec=getpos(mss[0])
        if os.path.isfile('big-mslist.txt'):
            status='downloaded'
        else:
            status='downloading'
        if os.path.isfile('image_dirin_SSD_init.tessel.reg'):
            status='started'
        if os.path.isfile('summary.txt'):
            status='complete'
        circles.append((name,ra,dec,s_colours[status]))
        print "%-30s %-15s %8.3f %8.3f %s" % (d,name,ra,dec,status)
        continue
    # else check for image only
    ims=glob.glob('image_full_ampphase1m.int.restored.fits')
    if len(ims)>0:
        ra,dec=getposim(ims[0])
        if name is None:
            name=d.split('/')[-1]
        circles.append((name,ra,dec,'magenta'))

ras=np.array([ra for _,ra,_,_ in circles])
decs=np.array([dec for _,_,dec,_ in circles])
rarange=max(ras)-min(ras)
decrange=max(decs)-min(decs)
yfigsize=6
xfigsize=yfigsize*(rarange/decrange)*np.cos(np.mean(decs)*np.pi/180.0)
print xfigsize,yfigsize

plt.figure(figsize=(xfigsize,yfigsize))

plt.xlim(np.min(ras)-csize,np.max(ras)+csize)
plt.ylim(np.min(decs)-csize,np.max(decs)+csize)


for name,ra,dec,colour in circles:
    plotcircle(name,ra,dec,colour)

plt.gca().invert_xaxis()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.show()
