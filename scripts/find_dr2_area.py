from __future__ import print_function
from builtins import range
import numpy as np
import os
from astropy_healpix import HEALPix
from astropy import units as u
import sys
from surveys_db import SurveysDB

#pos=np.loadtxt(os.environ['DDF_DIR']+'/ddf-pipeline/misc/DR2-pointings.txt',usecols=(1,2))

#if len(sys.argv)==1:
#    pos=pos[(pos[:,0]>137) & (pos[:,0]<250)]

hp = HEALPix(nside=1024)
print(hp.npix,'total healpix pixels on sky')
area=hp.pixel_area.value*3283
print('area of one healpix is',area,'sq. deg')

for archived in [False,True]:

    if archived:
        print('Doing archived area only')
        command='select ra,decl from fields where dr2>0 and status="Archived"'
    else:
        print('Doing full sky area')
        command='select ra,decl from fields where dr2>0'
    with SurveysDB(readonly=True) as sdb:
        sdb.cur.execute(command)
        result=sdb.cur.fetchall()

    print('Number of pointings is',len(result))
    ra=[]
    dec=[]
    for r in result:
        ra.append(r['ra'])
        dec.append(r['decl'])

    #ra=pos[:,0]
    #dec=pos[:,1]

    print('RA range is',np.min(ra),np.max(ra))
    print('Dec range is',np.min(dec),np.max(dec))

    plist=[]

    for i in range(len(ra)):
        pixels=hp.cone_search_lonlat(ra[i]*u.deg,dec[i]*u.deg,1.85*u.deg)
        plist=plist+list(pixels)

    print(len(plist),'total pixels')
    print(len(set(plist)),'total distinct pixels')
    print('Area is',len(set(plist))*area,'sq. deg')


