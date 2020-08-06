#!/usr/bin/env python

from __future__ import print_function
import sys
import os
from astropy.coordinates import SkyCoord,get_icrs_coordinates
import astropy.units as u
from surveys_db import SurveysDB
import numpy as np
from astropy.table import Table

def table_from_dict_list(l):
    ''' l is a list of dictionaries with identical string keys '''
    keys=l[0].keys()
    coldict={}
    for k in keys:
        klist=[]
        for d in l:
            klist.append(d[k])
        coldict[k]=klist
    return Table(coldict)

class Finder(object):
    def __init__(self):
        with SurveysDB(readonly=True) as sdb:
            sdb.cur.execute('select * from fields left join quality on fields.id=quality.id order by fields.id')
            results=sdb.cur.fetchall()
            self.t=table_from_dict_list(results)
            self.t['sc']=SkyCoord(self.t['ra'],self.t['decl'],unit=u.deg)
            
    def find(self,ra,dec,offset=4,verbose=False,check_obs=False):
        if check_obs:
            sdb=SurveysDB(readonly=True)

        sc=SkyCoord(ra,dec,unit='deg')
        minoffset=None
        t=self.t
        t['sep']=sc.separation(t['sc']).value
        tdet=t[t['sep']<=offset]
        for r in tdet:
            if check_obs:
                sdb.cur.execute('select * from observations where field="%s"' % r['id'])
                count=len(sdb.cur.fetchall())
                sdb.cur.execute('select * from observations where field="%s" and status="DI_processed"' % r['id'])
                proc_count=len(sdb.cur.fetchall())
            else:
                count=-1
                proc_count=-1
            if verbose: print('%-16s %-16s %2i %2i %8.3f %8.3f %6.3f %s' % (r['id'],r['status'],count,proc_count,r['ra'],r['decl'],r['sep'],r['location']))
            if r['status']=='Archived':
                if minoffset is None or r['sep']<minoffset:
                    minoffset=r['sep']
                    bestfield=r
        if check_obs:
            sdb.close()
            
        if minoffset is None:
            return None
        else:
            return bestfield
        
def find_pos(ra,dec,offset=4,verbose=True):
    # standalone wrapper
    f=Finder()
    result=f.find(ra,dec,offset=offset,verbose=True,check_obs=True)
    if result is None:
        return result
    else:
        return result['id']
    
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
        retval=find_pos(ra,dec)
    else:
        print('Call me with the name of an ILTJ source OR RA, Dec in degrees OR "object objectname"')
    if retval is not None:
        print('Return value was',retval)
