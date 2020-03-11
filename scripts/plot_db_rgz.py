#!/usr/bin/env python

from __future__ import print_function
from builtins import str
from builtins import range
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from surveys_db import SurveysDB
fontsize=16 # adjust to taste
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times'],'size':fontsize})
rc('text', usetex=True)
import sys
import datetime

org=180

def cc(ra,dec):
    ra=np.array(ra)

    x = np.remainder(ra+360-org,360) # shift RA values
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left

    return np.radians(x),np.radians(dec)

def plot_select(r,sf,label,**kwargs):

    ra=[]
    dec=[]
    r_in=[]
    r_out=[]
    for f in r:
        if sf(f):
            ra.append(f['ra'])
            dec.append(f['decl'])
            r_in.append(f)
        else:
            r_out.append(f)

    ra_r,dec_r=cc(ra,dec)
    plt.scatter(ra_r,dec_r,label=label,**kwargs)
    print("%-20s : %i" % (label,len(r_in)))
    return r_in,r_out
    
with SurveysDB(readonly=True) as sdb:
    sdb.cur.execute('select * from fields where dr2')
    #sdb.cur.execute('select * from fields where status!="Not started"')
    results=sdb.cur.fetchall()
    
print(len(results),'fields are in dr2')
        
fig = plt.figure(figsize=(16, 8))
fig.add_subplot(111, projection='aitoff')

# GP

for b in [-10,0,10]:

    lon=np.linspace(-180,180,1000)
    lat=b*np.ones_like(lon)

    sc=SkyCoord(l=lon,b=lat,unit=(u.deg,u.deg),frame='galactic')

    ra=np.array(sc.icrs.ra)
    dec=np.array(sc.icrs.dec)

    ra_r,dec_r=cc(ra,dec)
    
    plt.scatter(ra_r,dec_r,color='blue',s=5,label='MW' if b==0 else None)

#DR2 area
ravals = []
decvals = []
for r in results:
    ravals.append(r['ra'])
    decvals.append(r['decl'])

ra_r,dec_r=cc(ravals,decvals)
plt.scatter(ra_r,dec_r,marker='o',color='blue',alpha=0.2,zorder=-5,edgecolors='none',s=50,label='DR2')

_,r=plot_select(results,lambda r:r['dr1']==1,label='Complete (DR1)',color='magenta')
_,r=plot_select(r,lambda r:r['lgz']==1,label='Complete (internal)',color='orange')
_,r=plot_select(r,lambda r:r['gz_status']=='Complete',label='Complete',color='green')
_,r=plot_select(r,lambda r:r['gz_status'] in ['In progress'],label='In progress',color='cyan')
_,r=plot_select(r,lambda r:r['gz_status'] in ['Created','Downloading'],label='Ready to upload',color='yellow')
#_,r=plot_select(r,lambda r:r['gz_status'] in ['Failed'],label='Failed',color='red')


ax=plt.gca()

tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
tick_labels = np.remainder(tick_labels+360+org,360)
tick_labels = list(tick_labels)
for i in range(0,len(tick_labels)):
    tick_labels[i] = ''+str(tick_labels[i])+'\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \,'

ax.set_xticklabels(tick_labels,verticalalignment='top',rotation='vertical')

plt.xlabel('R.A.')
plt.ylabel('Decl.')
plt.grid(True)
plt.legend(loc=4)
plt.tight_layout()
plt.title('RGZ status at '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),loc='right')
if len(sys.argv)>1:
    plt.savefig(sys.argv[1],dpi=250)
else:
    plt.show()
