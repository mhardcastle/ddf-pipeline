from astropy.io import fits
import scipy.ndimage as nd
import numpy as np
import os
import sys

def make_blank_mask(indata,sizecut=20):
    '''
    take as input indata, a 2D array which may have small non-NaN islands in a sea of NaN, and return a mask of the same shape as indata which selects all those islands for NaN-blanking
    '''
    data = indata.byteswap().newbyteorder()
    data[data==0]==-9999
    data[np.isnan(data)]=0
    l,n=nd.label(data)
    size = np.bincount(l.ravel())
    #biggest_label = size[1:].argmax() + 1
    big=np.arange(len(size))[size>sizecut]
    big=big[1:]
    return ~np.isin(l,big) 

def blanker(infile,verbose=False):

    hdu=fits.open(infile)
    mask=blank_mask(hdu[0].data)
    hdu[0].data[mask]=np.nan

    os.rename(infile,infile.replace('.fits','-old.fits'))
    hdu.writeto(infile,overwrite=True)

if __name__=='__main__':
    for f in sys.argv[1:]:
        blanker(f)

