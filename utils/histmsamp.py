#!/usr/bin/python

# make and optionally plot a histogram of amplitudes of MSs

from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
import pyrap.tables as pt
import numpy as np
from astropy.io import fits

def sumdico(rootname):

    model=fits.open(rootname+'.int.model.fits')
    fluxim=fits.open(rootname+'.Norm.fits')

    cmodel=model[0].data*np.sqrt(fluxim[0].data)
    model.close()
    fluxim.close()

    return np.sum(cmodel)

def find_uvmin(listname,level,colname='CORRECTED_DATA',plot=False,tstep=30):

    if plot:
        import matplotlib.pyplot as plt

    mss=[l.rstrip() for l in open(listname).readlines()]

    minuv=0 # km
    maxuv=5
    bins=np.linspace(0,5,100)

    bvl=[]
    avl=[]
    for ms in mss:
        print('Doing',ms)
        t = pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False)
        freq = t[0]['REF_FREQUENCY']
        channels = t[0]['NUM_CHAN']
        chlist=t[0]['CHAN_FREQ']
        lamb=old_div(3e8,freq)
        max=100000
        print('Frequency is',freq,'Hz','wavelength is',lamb,'m')
        t.close()
        t=pt.table(ms, readonly=True, ack=False)
        newt=pt.taql('select UVW,TIME,'+colname+' from $t where FLAG[0,0]=False')
        time=newt.getcol('TIME')
        tuniq=np.unique(time)
        print('There are',len(tuniq),'unique times')
        if len(tuniq)<=tstep:
            print('Too much data flagged, skipping')
            continue
        tvals=tuniq[::tstep]
        print('Using',len(tvals),'times')
        uv=newt.getcol('UVW')/1000.0

        data=newt.getcol(colname)[:,:,0::3] # XX, YY
        stokesi=data[:,:,0]+data[:,:,1]
        adata=np.mean(np.absolute(stokesi),axis=1)

        newt.close()
        t.close()
        maxd=0
        for t in tvals:
            uvd=np.sqrt(np.sum(uv[time==t]**2.0,axis=1))
            amp=adata[time==t]
            bvals=np.digitize(uvd,bins)
            amps=amp[bvals<len(bins)]
            fbvals=bvals[bvals<len(bins)]
            bvl.append(fbvals)
            avl.append(amps)

    bvl=np.concatenate(bvl)
    avl=np.concatenate(avl)

    histv=np.zeros_like(bins)

    for i in range(len(bins)-1):
        if np.sum(bvl==i)==0: continue
        histv[i]=np.percentile(avl[bvl==i],100.0)
        if plot: plt.scatter(0.5*(bins[i]+bins[i+1]),histv[i])

    if plot:
        plt.xlabel('UV dist (km)')
        plt.ylabel('Apparent flux (Jy)')

    for i in range(len(bins)-1):
        if np.all(histv[i:]<level):
            print('Cut at',bins[i],'km')
            uvmin=bins[i]
            break

    if plot:
        plt.show()

    return uvmin

if __name__=='__main__':

    import sys
    level=sumdico(sys.argv[2])
    print('Total apparent flux is',level,'Jy')
    listname=sys.argv[1]
    
    print(find_uvmin(listname,level,plot=True))
