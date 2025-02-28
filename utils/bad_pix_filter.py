from astropy.io import fits
import scipy.ndimage as nd
import numpy as np
import os
import sys
from copy import copy

def make_blank_mask(indata,sizecut=20,verbose=True):
    '''
    take as input indata, a 2D array which may have small non-NaN islands in a sea of NaN, and return a mask of the same shape as indata which selects all those islands for NaN-blanking
    '''
    if verbose: print('Making blank mask')
    data=copy(indata)
    data[data==0]==-9999
    data[np.isnan(data)]=0
    l,n=nd.label(data)
    size = np.bincount(l.ravel())
    if verbose: print(size)
    big=np.arange(len(size))
    big=big[1:]
    size=size[1:]
    big=big[size>sizecut]
    if verbose: print(big)
    return ~np.isin(l,big) 

def blanker(infile,verbose=False):

    hdu=fits.open(infile)
    swapdata = hdu[0].data.byteswap().newbyteorder()
    mask=make_blank_mask(swapdata)
    hdu[0].data[mask]=np.nan

    os.rename(infile,infile.replace('.fits','-old.fits'))
    hdu.writeto(infile,overwrite=True)

if __name__=='__main__':
    for f in sys.argv[1:]:
        blanker(f)

