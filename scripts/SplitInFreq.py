#!/usr/bin/env python
import os
from pyrap.tables import table
import sys
import numpy as np


def split(MSName):
    tf=table("%s/SPECTRAL_WINDOW/"%MSName)
    _,nch=tf.getcol("CHAN_FREQ").shape

    nMS=20

    chStep=nch//nMS
    chEdges=np.linspace(0,nch,nMS+1)
    TemplateName="DP3_split.parset"
    f=open(TemplateName,"w")
    f.write("msin.datacolumn = DATA\n")
    f.write("msin.autoweight = false\n")
    f.write("steps = []\n")
    f.close()
    

    for iMS in range(nMS):
        ch0,ch1=chEdges[iMS],chEdges[iMS+1]
        nchThis=ch1-ch0
        ss="DP3 %s msin=%s msout=%s.chunk%03i.ms msin.startchan=%i msin.nchan=%i"%(TemplateName,MSName,MSName,iMS,ch0,nchThis)
        print("===============================================")
        print(ss)
        os.system(ss)

if __name__=="__main__":
    MSName=sys.argv[1]
    split(MSName)
