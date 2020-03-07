#!/usr/bin/env python
# Make MS lists, checking for heavily flagged data

from __future__ import print_function
import os
import glob
import pyrap.tables as pt
import numpy as np
from auxcodes import warn
from surveys_db import use_database,update_status

def check_flagged(ms):
    t = pt.table(ms, readonly=True)
    tc = t.getcol('FLAG').flatten()
    return float(np.sum(tc))/len(tc)

def get_timerange(ms):
    t = pt.table(ms +'/OBSERVATION', readonly=True, ack=False)
    return t.getcell('TIME_RANGE',0)

def make_list(workdir='.',force=False):
    g=sorted(glob.glob(workdir+'/*.ms'))
    full_mslist=[]
    start_times=[]
    for ms in g:
        ff=check_flagged(ms)
        t0,t1=get_timerange(ms)
        print(ms,ff)
        if ff<0.8:
            full_mslist.append(os.path.basename(ms))
            start_times.append(t0)
    full_mslist = np.array(full_mslist)
            
    # check for multiple observations
    Ustart_times = np.unique(start_times)

    if len(full_mslist)<18:
        warn('Too few MS found for normal running: only %i' % len(full_mslist))
        if not force: return False

    if len(full_mslist)<24:
        warn('Warning -- only %i ms found' % len(full_mslist))
        
    # ensure lists contain all ms from all observations and same subset for each observation
    write_full_mslist = np.array(())
    write_mslist = np.array(())
    for start_time in Ustart_times:        
        write_full_mslist = np.hstack((write_full_mslist,full_mslist[start_times==start_time]))
        write_mslist = np.hstack((write_mslist,full_mslist[start_times==start_time][2::4]))

    open(workdir+'/big-mslist.txt','w').writelines(ms+'\n' for ms in write_full_mslist)
    open(workdir+'/mslist.txt','w').writelines(ms+'\n' for ms in write_mslist)
    return True

def list_db_update(success,workdir=None):
    if success:
        update_status(None,'Ready',workdir=workdir)
    else:
        update_status(None,'List failed',workdir=workdir)

if __name__=='__main__':
    import sys
    if len(sys.argv)>1:
        force=sys.argv[1]=='force'
        print('Force is',force)
    else:
        force=False
    success=make_list(workdir=os.getcwd(),force=force)
    if use_database():
        list_db_update(success)
