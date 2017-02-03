#!/usr/bin/env python
# Make MS lists, checking for heavily flagged data

import glob
import pyrap.tables as pt
import numpy as np
from auxcodes import warn

def check_flagged(ms):
    t = pt.table(ms)
    tc = t.getcol('FLAG').flatten()
    return float(np.sum(tc))/len(tc)

def make_list():
    g=sorted(glob.glob('*.ms'))
    full_mslist=[]
    for ms in g:
        ff=check_flagged(ms)
        print ms,ff
        if ff<0.8:
            full_mslist.append(ms)
    full_mslist = np.array(full_mslist)
            
    # check for multiple observations
    obsids = np.array([ms.split('_')[0] for ms in full_mslist])
    Uobsids = np.unique(obsids)

    if len(full_mslist)<24:
        warn('Warning -- only %i ms found' % len(full_mslist))
        
    # ensure lists contain all ms from all observations and same subset for each observation
    write_full_mslist = np.array(())
    write_mslist = np.array(())
    for obsid in Uobsids:        
        write_full_mslist = np.hstack((write_full_mslist,full_mslist[obsids==obsid]))
        write_mslist = np.hstack((write_mslist,full_mslist[obsids==obsid][2::4]))

    open('big-mslist.txt','w').writelines(ms+'\n' for ms in write_full_mslist)
    open('mslist.txt','w').writelines(ms+'\n' for ms in write_mslist)

if __name__=='__main__':
    make_list()
