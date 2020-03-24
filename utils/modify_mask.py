#!/usr/bin/env python

from __future__ import division
from builtins import zip
from past.utils import old_div
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import numpy as np

def ellipse(x,y,xp,yp,major,minor,angle):
    #tests if a point[xp,yp] is within
    #boundaries defined by the ellipse
    #of center[x,y], diameter d D, and tilted at angle

    cosa=np.cos(old_div(angle*np.pi,180))
    sina=np.sin(old_div(angle*np.pi,180))
    dd=(minor/2.0)**2.0
    DD=(major/2.0)**2.0

    a=np.power(cosa*(xp-x)+sina*(yp-y),2)
    b=np.power(sina*(xp-x)-cosa*(yp-y),2)
    return (old_div(a,dd))+(old_div(b,DD))

def modify_mask(infile,outfile,table,radius,fluxlim,save_filtered=None,do_extended=False,cellsize=1.5,pointsize=30.0):
    """Take a pre-existing mask file, in infile: find all entries in FITS
    table table that lie in the map region, and add their positions to
    the mask with a fixed radius radius in pixels: write the mask out
    to outfile, which may be the same as infile.

    save_filtered, if not None, should be a filename to save the
    filtered table to.

    """

    t=Table.read(table)
    t=t[t['Peak_flux']>fluxlim]
    if len(t)==0:
        raise Exception('Flux-filtered table is zero-length. Check your table fluxes and/or positions')
    hdu=fits.open(infile)
    w=WCS(hdu[0].header)
    map=hdu[0].data
    mask=np.zeros_like(map,dtype=int)
    _,_,ymax,xmax=map.shape
    x,y,_,_=w.wcs_world2pix(t['RA'],t['DEC'],0,0,0)
    filter=(x>=0) & (x<xmax) & (y>=0) & (y<ymax)
    t=t[filter]
    x=x[filter]
    y=y[filter]
    for i,(xv,yv) in enumerate(zip(x,y)):
        do_ellipse=False
        r=radius
        if do_extended:
            # check whether the major axis exceeds a limit
            major=t[i]['Maj']
            minor=t[i]['Min']
            pa=t[i]['PA']
            if major>pointsize:
                do_ellipse=True
                r=major+radius
                # print 'Doing ellipse with',major,minor,pa

        cxmin=xv-r-1
        if cxmin<0: cxmin=0
        cxmax=xv+r+1
        if cxmax>xmax: cxmax=xmax
        cymin=yv-r-1
        if cymin<0: cymin=0
        cymax=yv+r+1
        if cymax>ymax: cymax=ymax
        X, Y = np.meshgrid(np.arange(cxmin,cxmax,1.0), np.arange(cymin,cymax,1.0))
        if do_ellipse:
            ellv=ellipse(xv,yv,X,Y,old_div(major,cellsize)+radius*2.0,old_div(minor,cellsize)+radius*2.0,pa)
            mask[0,0,Y[ellv<1.0].astype(int),X[ellv<1.0].astype(int)]=1
        else:
            rv=np.sqrt((X+0.5-xv)**2.0+(Y+0.5-yv)**2.0)
            mask[0,0,Y[rv<radius].astype(int),X[rv<radius].astype(int)]=1

    hdu[0].data=(map.astype(int) | mask).astype(np.float32)
    hdu.writeto(outfile,overwrite=True)
    if save_filtered is not None:
        t.write(save_filtered,overwrite=True)
        

if __name__=='__main__':
    import sys

    infile=sys.argv[1]
    outfile=sys.argv[2]
    table=sys.argv[3]
    radius=float(sys.argv[4])
    try:
        sf=sys.argv[5]
    except:
        sf=None
    modify_mask(infile,outfile,table,radius,500,save_filtered=sf,do_extended=True)
