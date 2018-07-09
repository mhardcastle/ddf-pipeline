#!/usr/bin/env python
# Make MS lists, checking for heavily flagged data

import glob
import pyrap.tables as pt
import numpy as np
from auxcodes import warn
from surveys_db import SurveysDB,use_database,get_id,tag_idd

def check_flagged(ms):
    t = pt.table(ms, readonly=True)
    tc = t.getcol('FLAG').flatten()
    return float(np.sum(tc))/len(tc)

def get_timerange(ms):
    t = pt.table(ms +'/OBSERVATION', readonly=True, ack=False)
    return t.getcell('TIME_RANGE',0)

def getpos(ms):
    t = pt.table(ms+ '/OBSERVATION', readonly=True, ack=False)
    name=t[0]['LOFAR_TARGET']
    (start,end)=t[0]['TIME_RANGE']
    t = pt.table(ms+ '/FIELD', readonly=True, ack=False)
    direction = t[0]['PHASE_DIR']
    ra, dec = direction[0]
    if (ra<0):
        ra+=2*np.pi
    t = pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False)
    channels=t[0]['NUM_CHAN']
    t=pt.table(ms, readonly=True, ack=False)
    time=t.getcol('TIME')
    tv=np.sort(np.unique(time))
    dt=tv[1]-tv[0]
    
    return name[0],ra*(180/np.pi),dec*(180/np.pi),channels,start,dt

def make_list():
    g=sorted(glob.glob('*.ms'))
    full_mslist=[]
    start_times=[]
    for ms in g:
        ff=check_flagged(ms)
        t0,t1=get_timerange(ms)
        print ms,ff
        if ff<0.8:
            full_mslist.append(ms)
            start_times.append(t0)
    full_mslist = np.array(full_mslist)
            
    # check for multiple observations
    Ustart_times = np.unique(start_times)

    if len(full_mslist)<18:
        warn('Too few MS found for normal running: only %i' % len(full_mslist))
        return False

    if len(full_mslist)<24:
        warn('Warning -- only %i ms found' % len(full_mslist))
        
    # ensure lists contain all ms from all observations and same subset for each observation
    write_full_mslist = np.array(())
    write_mslist = np.array(())
    for start_time in Ustart_times:        
        write_full_mslist = np.hstack((write_full_mslist,full_mslist[start_times==start_time]))
        write_mslist = np.hstack((write_mslist,full_mslist[start_times==start_time][2::4]))

    open('big-mslist.txt','w').writelines(ms+'\n' for ms in write_full_mslist)
    open('mslist.txt','w').writelines(ms+'\n' for ms in write_mslist)
    return True

def list_db_update(success):
    id=get_id()
    sdb=SurveysDB()
    idd=sdb.get_id(id)
    if success:
        idd['status']='Ready'
        tag_idd(sdb,idd)
        msname=open('mslist.txt').readlines()[0].rstrip()
        name,ra,dec,nchan,start,dt=getpos(msname)
        idd['fieldname']=name
        idd['ra']=ra
        idd['decl']=dec
        idd['nchan']=nchan/10
        idd['date']=start
        idd['integration']=dt
    else:
        idd['status']='List failed'
    sdb.set_id(idd)

if __name__=='__main__':
    success=make_list()
    if use_database():
        list_db_update(success)
