#!/usr/bin/python

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import numpy as np

def modify_mask(infile,outfile,table,radius,fluxlim,save_filtered=None):
    """Take a pre-existing mask file, in infile: find all entries in FITS
    table table that lie in the map region, and add their positions to
    the mask with a fixed radius radius in pixels: write the mask out
    to outfile, which may be the same as infile.

    save_filtered, if not None, should be a filename to save the
    filtered table to.

    """

    t=Table.read(table)
    t=t[t['Peak_flux']>fluxlim]
    hdu=fits.open(infile)
    w=WCS(hdu[0].header)
    map=hdu[0].data
    mask=np.zeros_like(map,dtype=int)
    _,_,ymax,xmax=map.shape
    x,y,_,_=w.wcs_world2pix(t['RA'],t['DEC'],0,0,0)
    filter=(x>=0) & (x<xmax) & (y>=0) & (y<ymax)
    x=x[filter]
    y=y[filter]
    for xv,yv in zip(x,y):
        cxmin=xv-radius-1
        if cxmin<0: cxmin=0
        cxmax=xv+radius+1
        if cxmin>xmax: cxmax=xmax
        cymin=yv-radius-1
        if cymin<0: cymin=0
        cymax=yv+radius+1
        if cymin>ymax: cymax=ymax
        X, Y = np.meshgrid(np.arange(cxmin,cxmax,1.0), np.arange(cymin,cymax,1.0))
        rv=np.sqrt((X+0.5-xv)**2.0+(Y+0.5-yv)**2.0)
        mask[0,0,Y[rv<radius].astype(int),X[rv<radius].astype(int)]=1

    hdu[0].data=(map.astype(int) | mask).astype(np.float32)
    hdu.writeto(outfile,clobber=True)
    if save_filtered is not None:
        newt=t[filter]
        newt.write(save_filtered,overwrite=True)
        

if __name__=='__main__':
    import sys

    infile=sys.argv[1]
    outfile=sys.argv[2]
    table=sys.argv[3]
    radius=int(sys.argv[4])
    try:
        sf=sys.argv[5]
    except:
        sf=None
    modify_mask(infile,outfile,table,radius,500,save_filtered=sf)
