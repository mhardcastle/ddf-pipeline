#!/usr/bin/python

import sys
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
from surveys_db import SurveysDB
from auxcodes import sepn
import numpy as np

factor=180.0/np.pi

def separation(ra1,dec1,ra2,dec2):
    # same as sepn but in degrees
    return factor*sepn(ra1/factor,dec1/factor,ra2/factor,dec2/factor)

def find_pos(ra,dec,name=None,offset=4):
    raoffset=offset/np.cos(dec/factor)
    with SurveysDB() as sdb:
        sdb.cur.execute('select * from fields where ra>%f and ra<%f and decl>%f and decl<%f' % (ra-raoffset,ra+raoffset,dec-offset,dec+offset))
        results=sdb.cur.fetchall()
        for r in results:
            sdb.cur.execute('select * from observations where field="%s"' % r['id'])
            count=len(sdb.cur.fetchall())
            sdb.cur.execute('select * from observations where field="%s" and status="DI_processed"' % r['id'])
            proc_count=len(sdb.cur.fetchall())
            print '%-16s %-16s %2i %2i %8.3f %8.3f %6.3f %s' % (r['id'],r['status'],count,proc_count,r['ra'],r['decl'],separation(ra,dec,r['ra'],r['decl']),r['location'])

if __name__=='__main__':
    
    if len(sys.argv)==3:
        ra=float(sys.argv[1])
        dec=float(sys.argv[2])
        find_pos(ra,dec)
    elif len(sys.argv)==2:
        s=sys.argv[1][4:]
        coord=s[0:2]+':'+s[2:4]+':'+s[4:9]+' '+s[9:12]+':'+s[12:14]+':'+s[14:]
        sc = SkyCoord(coord,unit=(u.hourangle,u.deg))
        ra=sc.ra.value
        dec=sc.dec.value
        print 'Parsed coordinates to ra=%f, dec=%f' % (ra,dec)
        name=sys.argv[1]
        find_pos(ra,dec,name=name)
    else:
        print 'Call me with the name of a source or RA, Dec in degrees'
