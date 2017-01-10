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
for d in sys.argv[1:]:
    os.chdir(d)
    mss=glob.glob('*.ms')
    name=None
    if len(mss)>0:
        name,ra,dec=getpos(mss[0])
        print d,name,ra,dec
        circles.append((name,ra,dec,'blue'))
    
    ims=glob.glob('image_full_ampphase1m.smooth.int.restored.fits')
    if len(ims)>0:
        ra,dec=getposim(ims[0])
        if name is None:
            name=d.split('/')[-1]
        circles.append((name,ra,dec,'red'))

ras=np.array([ra for _,ra,_,_ in circles])
decs=np.array([dec for _,_,dec,_ in circles])

plt.figure(figsize=(10/np.cos(np.mean(decs)*np.pi/180.0),10))

plt.xlim(np.min(ras)-csize,np.max(ras)+csize)
plt.ylim(np.min(decs)-csize,np.max(decs)+csize)


for name,ra,dec,colour in circles:
    plotcircle(name,ra,dec,colour)

plt.gca().invert_xaxis()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.show()
