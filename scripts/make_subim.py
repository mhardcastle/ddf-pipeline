#!/usr/bin/env python

# Make cutouts of pipeline products

# from git/lotss-catalogue/utils/subim.py

from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sys
from astropy.coordinates import SkyCoord,get_icrs_coordinates
import astropy.units as u

def cube_extract(f,ra,dec,x,y,size,hduid=0,verbose=True):
    """
    Like flatten, but don't flatten. So just make the appropriate slice
    """
    naxis=f[hduid].header['NAXIS']
    if naxis<=2:
        raise RuntimeError('Can\'t make cube from this')

    if verbose:
        print(f[hduid].data.shape)
    ds=f[hduid].data.shape[-2:]
    by,bx=ds
    xmin=int(x-size)
    if xmin<0:
        xmin=0
    xmax=int(x+size)
    if xmax>bx:
        xmax=bx
    ymin=int(y-size)
    if ymin<0:
        ymin=0
    ymax=int(y+size)
    if ymax>by:
        ymax=by
    
    if ymax<=ymin or xmax<=xmin:
        # this can only happen if the required position is not on the map
        print(xmin,xmax,ymin,ymax)
        raise RuntimeError('Failed to make subimage!')

    w = WCS(f[hduid].header)
    w.wcs.crpix[0]=w.wcs.crpix[0]-xmin
    w.wcs.crpix[1]=w.wcs.crpix[1]-ymin

    header = w.to_header()
    slice=[]
    for i in range(naxis,0,-1):
        if i==1:
            slice.append(np.s_[xmin:xmax])
        elif i==2:
            slice.append(np.s_[ymin:ymax])
        else:
            slice.append(np.s_[:])
    if verbose:
        print(slice)

    hdu=fits.PrimaryHDU(f[hduid].data[slice],header)
    copy=('EQUINOX','EPOCH','BMAJ','BMIN','BPA','RESTFRQ','DATE-OBS','TELESCOP','OBJECT','DATE-MAP','ORIGIN')
    for k in copy:
        r=f[hduid].header.get(k)
        if r:
            hdu.header[k]=r
    if 'TAN' in hdu.header['CTYPE1']:
        hdu.header['LATPOLE']=f[hduid].header['CRVAL2']
    hdulist=fits.HDUList([hdu])
    return hdulist
    
    
def flatten(f,ra,dec,x,y,size,hduid=0,channel=0,freqaxis=3,verbose=True):
    """ 
    Flatten a fits file so that it becomes a 2D image. Return new header and
    data
    This version also makes a sub-image of specified size.
    """

    naxis=f[hduid].header['NAXIS']
    if naxis<2:
        raise RuntimeError('Can\'t make map from this')

    if verbose:
        print(f[hduid].data.shape)
    ds=f[hduid].data.shape[-2:]
    by,bx=ds
    xmin=int(x-size)
    if xmin<0:
        xmin=0
    xmax=int(x+size)
    if xmax>bx:
        xmax=bx
    ymin=int(y-size)
    if ymin<0:
        ymin=0
    ymax=int(y+size)
    if ymax>by:
        ymax=by
    
    if ymax<=ymin or xmax<=xmin:
        # this can only happen if the required position is not on the map
        print(xmin,xmax,ymin,ymax)
        raise RuntimeError('Failed to make subimage!')

    w = WCS(f[hduid].header)
    wn=WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]-xmin
    wn.wcs.crpix[1]=w.wcs.crpix[1]-ymin
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    try:
        wn.wcs.pc=w.wcs.pc[0:2,0:2]
    except AttributeError:
        pass # pc is not present
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2

    slice=[]
    for i in range(naxis,0,-1):
        if i==1:
            slice.append(np.s_[xmin:xmax])
        elif i==2:
            slice.append(np.s_[ymin:ymax])
        elif i==freqaxis:
            slice.append(channel)
        else:
            slice.append(0)
    if verbose:
        print(slice)

    hdu=fits.PrimaryHDU(f[hduid].data[slice],header)
    copy=('EQUINOX','EPOCH','BMAJ','BMIN','BPA','RESTFRQ','DATE-OBS','TELESCOP','OBJECT','DATE-MAP','ORIGIN','BTYPE','BUNIT')
    for k in copy:
        r=f[hduid].header.get(k)
        if r:
            hdu.header[k]=r
    if 'TAN' in hdu.header['CTYPE1']:
        hdu.header['LATPOLE']=f[hduid].header['CRVAL2']
    hdulist=fits.HDUList([hdu])
    return hdulist

def extract_subim(filename,ra,dec,size,hduid=0,verbose=True,cubemode=False):
    if verbose:
        print('Opening',filename)
    orighdu=fits.open(filename)
    psize=int(old_div(size,orighdu[hduid].header['CDELT2']))
    if verbose:
        print('pix size is',psize,size)
    ndims=orighdu[hduid].header['NAXIS']
    pvect=np.zeros((1,ndims))
    lwcs=WCS(orighdu[hduid].header)
    pvect[0][0]=ra
    pvect[0][1]=dec
    imc=lwcs.wcs_world2pix(pvect,0)
    x=imc[0][0]
    y=imc[0][1]
    if cubemode:
        hdu=cube_extract(orighdu,ra,dec,x,y,psize,hduid=hduid,verbose=verbose)
    else:
        hdu=flatten(orighdu,ra,dec,x,y,psize,hduid=hduid,verbose=verbose)
    return hdu

def extract_and_save(filename,ra,dec,size,outname='cutout.fits',cubemode=False,scale=None):
    hdu=extract_subim(filename,ra,dec,size,verbose=False,cubemode=cubemode)
    if scale is not None:
        hdu[0].data*=scale
    hdu.writeto(outname,overwrite=True)


if __name__=='__main__':
    filename=sys.argv[1]
    size=float(sys.argv[2])
    if len(sys.argv)==5:
        try:
            ra=float(sys.argv[3])
            dec=float(sys.argv[4])
        except ValueError:
            if sys.argv[3]=='object':
                c=get_icrs_coordinates(sys.argv[4])
            else:
                c = SkyCoord(sys.argv[3],sys.argv[4], frame='icrs',unit=(u.hourangle, u.deg))
            ra=float(c.ra.degree)
            dec=float(c.dec.degree)
            print(ra,dec)
        extract_and_save(filename,ra,dec,size)
    elif len(sys.argv)==4:
        s=sys.argv[3][4:]
        coord=s[0:2]+':'+s[2:4]+':'+s[4:9]+' '+s[9:12]+':'+s[12:14]+':'+s[14:]
        sc = SkyCoord(coord,unit=(u.hourangle,u.deg))
        ra=sc.ra.value
        dec=sc.dec.value
        print('Parsed coordinates to ra=%f, dec=%f' % (ra,dec))
        name=sys.argv[1]
        extract_and_save(filename,ra,dec,size)
    else:
        print('Call: filename size RA DEC _OR_ filename size ILTname.')

