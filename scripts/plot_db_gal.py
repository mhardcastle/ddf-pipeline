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
    sdb.cur.execute('select fields.id as id,gal_l as ra,gal_b as decl,fields.status as status,observations.status as ostatus,observations.location as location,sum(nsb*integration/232) as s,count(observations.id) as c,fields.priority from fields left join observations on (observations.field=fields.id) group by fields.id having ostatus is not null')
    #sdb.cur.execute('select * from fields where status!="Not started"')
    results=sdb.cur.fetchall()

print(len(results),'fields have some observations')
        
fig = plt.figure(figsize=(16, 8))
fig.add_subplot(111, projection='aitoff')

# GP

for b in [-10,0,10]:

    lon=np.linspace(-180,180,1000)
    lat=b*np.ones_like(lon)

    sc=SkyCoord(ra=lon,dec=lat,unit=(u.deg,u.deg),frame='icrs')

    ra=np.array(sc.galactic.l)
    dec=np.array(sc.galactic.b)

    ra_r,dec_r=cc(ra,dec)
    
    plt.scatter(ra_r,dec_r,color='blue',s=5,label='CE' if b==0 else None)

#DR2 area
ravals = []
decvals = []
infile = open(os.environ['DDF_DIR']+'/ddf-pipeline/misc/DR2-pointings.txt','r')
for line in infile:
    line = line[:-1]
    line = line.split(' ')
    while '' in line:
        line.remove('')
    pointing = line[0]
    ravals.append(float(line[1]))
    decvals.append(float(line[2]))
infile.close()
sc=SkyCoord(ra=ravals,dec=decvals,unit=(u.deg,u.deg),frame='icrs')
ravals=np.array(sc.galactic.l)
decvals=np.array(sc.galactic.b)
ra_r,dec_r=cc(ravals,decvals)
plt.scatter(ra_r,dec_r,marker='o',color='blue',alpha=0.2,zorder=-5,edgecolors='none',s=50,label='DR2')
    
_,r=plot_select(results,lambda r:r['status'] in ['Archived','Complete'],label='Complete',color='green')
_,r=plot_select(r,lambda r:r['status'] in ['Running'],label='Running',color='cyan')
_,r=plot_select(r,lambda r:r['status'] in ['Downloaded','Downloading','Unpacking','Averaging','Ready','Queued','Unpacked'],label='In progress',color='yellow')
_,r=plot_select(r,lambda r:r['status'] in ['Failed','List failed','D/L failed'],label='Failed',color='red')

_,r=plot_select(r,lambda r:r['ostatus']=='DI_processed' and r['s']>7,label='Ready',color='black')

_,r=plot_select(r,lambda r:r['location']!="Sara",label='Observed (not Sara)',color='red',alpha=1.0,s=5)
_,r=plot_select(r,lambda r:True,label='Observed',color='black',alpha=0.5,s=5)


ax=plt.gca()

tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
tick_labels = np.remainder(tick_labels+360+org,360)
tick_labels = list(tick_labels)
for i in range(0,len(tick_labels)):
    tick_labels[i] = ''+str(tick_labels[i])+'\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \,'

ax.set_xticklabels(tick_labels,verticalalignment='top',rotation='vertical')

plt.xlabel('$\ell$')
plt.ylabel('$b$')
plt.grid(True)
plt.legend(loc=4)
plt.tight_layout()
plt.title('DR2 processing status at '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),loc='right')

plt.show()
