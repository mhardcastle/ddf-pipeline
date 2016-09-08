#!/usr/bin/python

from astropy.io import fits
import numpy as np

def make_cube(outname,imlist,freqs=None):
    # imlist must be in the desired frequency order
    # dummy frequency header will be created if not present; this need not be correct but PyBDSM won't work without it.

    if freqs is None:
        freqs=[1400e9,1401e9]

    hdus=[]
    for image in imlist:
        hdus.append(fits.open(image))

    d1,d2,y,x=hdus[0][0].data.shape

    if d1>1 or d2>1:
        raise Exception('Images to form cube from must be single-channel and single-Stokes')

    newdata=np.zeros((len(hdus),1,y,x),dtype=np.float32)

    for i,h in enumerate(hdus):
        newdata[i,0,:,:]=h[0].data[0,0,:,:]

    hdus[0][0].data=newdata
    hdus[0][0].header['NAXIS4']=len(hdus)
    hdus[0][0].header['CTYPE4']='FREQ'
    hdus[0][0].header['CRPIX4']=1.0
    hdus[0][0].header['CRVAL4']=freqs[0]
    hdus[0][0].header['CDELT4']=freqs[1]-freqs[0] # need not be correct
    hdus[0][0].header['CUNIT4']='Hz'
    hdus[0][0].header['NAXIS3']=len(hdus)
    hdus[0][0].header['CTYPE3']='STOKES'
    hdus[0][0].header['CRPIX3']=1.0
    hdus[0][0].header['CRVAL3']=1.0
    hdus[0][0].header['CDELT3']=1.0
    hdus[0][0].header['CUNIT3']=''
    hdus[0].writeto(outname,clobber=True)
    for h in hdus: h.close()

if __name__=='__main__':
    import sys

    make_cube(sys.argv[1],sys.argv[2:])
