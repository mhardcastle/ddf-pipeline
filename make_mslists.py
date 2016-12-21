#!/usr/bin/env python
# Make MS lists, checking for heavily flagged data

import glob
import pyrap.tables as pt
import numpy as np

def check_flagged(ms):
    t = pt.table(ms)
    tc = t.getcol('FLAG').flatten()
    return float(np.sum(tc))/len(tc)

g=glob.glob('*.ms')
full_mslist=[]
for ms in g:
    ff=check_flagged(ms)
    print ms,ff
    if ff<0.8:
        full_mslist.append(ms)

open('big_mslist.txt','w').writelines(ms+'\n' for ms in full_mslist)
open('mslist.txt','w').writelines(ms+'\n' for ms in full_mslist[2::4])
