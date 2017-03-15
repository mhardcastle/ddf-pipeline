#!/usr/bin/env python

import glob
from astropy.io import fits

def stack_bootstrap():
    g=glob.glob('image_low_*_SSD.app.restored.fits')
    hdus=[]
    for f in g:
        hdus.append(fits.open(f))
    for h in hdus[1:]:
        hdus[0][0].data+=h[0].data
    hdus[0][0].data/=len(hdus)
    hdus[0].writeto('bootstrap_stack.fits',clobber=True)

if __name__=='__main__':
    stack_bootstrap()

