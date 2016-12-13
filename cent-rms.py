#!/usr/bin/python

# find the central rms of an image by iteratively removing the outliers

import numpy as np
from astropy.io import fits

def get_rms(hdu,boxsize=1000,niter=20,eps=1e-6,verbose=False):

    data=hdu[0].data
    _,_,ys,xs=data.shape
    subim=data[0,0,ys/2-boxsize/2:ys/2+boxsize/2,xs/2-boxsize/2:xs/2+boxsize/2].flatten()
    oldrms=1
    for i in range(niter):
        rms=np.std(subim)
        if verbose: print len(subim),rms
        if np.abs(oldrms-rms)/rms < eps:
            return rms
        subim=subim[np.abs(subim)<5*rms]
        oldrms=rms
    raise Exception('Failed to converge')

if __name__=='__main__':
    import sys
    for name in sys.argv[1:]:
        hdu=fits.open(name)
        print name,get_rms(hdu)
