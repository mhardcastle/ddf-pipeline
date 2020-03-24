#!/usr/bin/env python

from __future__ import print_function
import sys
import os
from astropy.coordinates import SkyCoord,get_icrs_coordinates
import astropy.units as u
from surveys_db import SurveysDB
import numpy as np

def find_pos(ra,dec,offset=4,name=None,verbose=True):
    sc=SkyCoord(ra,dec,unit='deg')
    minoffset=None
    with SurveysDB() as sdb:
        sdb.cur.execute('select * from fields')
        results=sdb.cur.fetchall()
        ras=[r['ra'] for r in results]
        decs=[r['decl'] for r in results]
        fsc=SkyCoord(ras,decs,unit='deg')
        seps=sc.separation(fsc).value
        for i,r in enumerate(results):
            if seps[i]>offset: continue
            sdb.cur.execute('select * from observations where field="%s"' % r['id'])
            count=len(sdb.cur.fetchall())
            sdb.cur.execute('select * from observations where field="%s" and status="DI_processed"' % r['id'])
            proc_count=len(sdb.cur.fetchall())
            sep=seps[i]
            print('%-16s %-16s %2i %2i %8.3f %8.3f %6.3f %s' % (r['id'],r['status'],count,proc_count,r['ra'],r['decl'],sep,r['location']))
            if r['status']=='Archived':
                if minoffset is None or sep<minoffset:
                    minoffset=sep
                    bestfield=r['id']
    if minoffset is None:
        return None
    else:
        return bestfield

if __name__=='__main__':

    retval=None
    if len(sys.argv)==3:
        try:
            ra=float(sys.argv[1])
            dec=float(sys.argv[2])
        except ValueError:
            if sys.argv[1]=='object':
                c=get_icrs_coordinates(sys.argv[2])
            else:
                c = SkyCoord(sys.argv[1],sys.argv[2], frame='icrs',unit=(u.hourangle, u.deg))
            ra=float(c.ra.degree)
            dec=float(c.dec.degree)
            print(ra,dec)
        retval=find_pos(ra,dec)
    elif len(sys.argv)==2:
        s=sys.argv[1][4:]
        coord=s[0:2]+':'+s[2:4]+':'+s[4:9]+' '+s[9:12]+':'+s[12:14]+':'+s[14:]
        sc = SkyCoord(coord,unit=(u.hourangle,u.deg))
        ra=sc.ra.value
        dec=sc.dec.value
        print('Parsed coordinates to ra=%f, dec=%f' % (ra,dec))
        name=sys.argv[1]
        retval=find_pos(ra,dec,name=name)
    else:
        print('Call me with the name of an ILTJ source OR RA, Dec in degrees OR "object objectname"')
    if retval is not None:
        print('Return value was',retval)
