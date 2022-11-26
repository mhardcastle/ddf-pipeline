from __future__ import print_function
from builtins import range
from reproject.utils import parse_input_data, parse_output_projection
try:
    from reproject.interpolation.core_celestial import _reproject_celestial as reproj_interp
except ModuleNotFoundError:
    from reproject.interpolation.core import _reproject_full as reproj_interp
try:
    from reproject.spherical_intersect.core import _reproject_celestial as reproj_exact
except ModuleNotFoundError:
    from reproject.interpolation.core import _reproject_full as reproj_exact    
from reproject.interpolation.high_level import ORDER

import six

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from multiprocessing import Pool
import sys

def reproject_interp_chunk_2d(input_data, output_projection, shape_out=None, hdu_in=0,
                              order='bilinear', blocks=(1000, 1000), parallel=False):
    """
    For a 2D image, reproject in chunks
    """

    array_in, wcs_in = parse_input_data(input_data, hdu_in=hdu_in)
    wcs_out, shape_out = parse_output_projection(output_projection, shape_out=shape_out)

    if isinstance(order, six.string_types):
        order = ORDER[order]

    # find corners
    xs=[]
    ys=[]
    ymax,xmax=input_data.data.shape
    for x,y in ((0,0),(xmax,0),(xmax,ymax),(0,ymax)):
        sky=wcs_in.wcs_pix2world(x,y,0)
        xt,yt=wcs_out.wcs_world2pix(sky[0],sky[1],0)
        xs.append(xt)
        ys.append(yt)

    xtmin=int(np.min(xs)/1.1)
    if xtmin<0: xtmin=0
    xtmax=int(np.max(xs)*1.1)
    if xtmax>shape_out[1]: xtmax=shape_out[1]
    ytmin=int(np.min(ys)/1.1)
    if ytmin<0: ytmin=0
    ytmax=int(np.max(ys)*1.1)
    if ytmax>shape_out[0]: ytmax=shape_out[0]

    print(xtmin,xtmax,ytmin,ytmax)

    print('There will be',int(round((ytmax-ytmin)/blocks[0])*round((xtmax-xtmin)/blocks[1])),'chunks')
    
    # Create output arrays
    array = np.nan*np.ones(shape_out, dtype='float32')
    footprint = np.zeros(shape_out, dtype='float32')
    for imin in range(ytmin, ytmax, blocks[0]):
        imax = min(imin + blocks[0], array.shape[0])
        for jmin in range(xtmin, xtmax, blocks[1]):
            print('.', end=' ')
            sys.stdout.flush()
            jmax = min(jmin + blocks[1], array.shape[1])
            shape_out_sub = (imax - imin, jmax - jmin)
            wcs_out_sub = wcs_out.deepcopy()
            wcs_out_sub.wcs.crpix[0] -= jmin
            wcs_out_sub.wcs.crpix[1] -= imin
            array_sub, footprint_sub = reproj_interp(array_in, wcs_in, wcs_out_sub,
                                                            shape_out=shape_out_sub,
                                                            order=order)
            array[imin:imax, jmin:jmax] = array_sub
            footprint[imin:imax, jmin:jmax] = footprint_sub

    print() 
    return array, footprint

def reproject_exact_chunk_2d(input_data, output_projection, shape_out=None, hdu_in=0,
                              order='bilinear', blocks=(1000, 1000), parallel=False):
    """
    For a 2D image, reproject in chunks
    """

    array_in, wcs_in = parse_input_data(input_data, hdu_in=hdu_in)
    wcs_out, shape_out = parse_output_projection(output_projection, shape_out=shape_out)

    if isinstance(order, six.string_types):
        order = ORDER[order]

    # Create output arrays
    array = np.zeros(shape_out, dtype=float)
    footprint = np.zeros(shape_out, dtype=float)

    for imin in range(0, array.shape[0], blocks[0]):
        imax = min(imin + blocks[0], array.shape[0])
        for jmin in range(0, array.shape[1], blocks[1]):
            print('.', end=' ')
            jmax = min(jmin + blocks[1], array.shape[1])
            shape_out_sub = (imax - imin, jmax - jmin)
            wcs_out_sub = wcs_out.deepcopy()
            wcs_out_sub.wcs.crpix[0] -= jmin
            wcs_out_sub.wcs.crpix[1] -= imin
            array_sub, footprint_sub = reproj_exact(array_in, wcs_in, wcs_out_sub,
                                                            shape_out=shape_out_sub,
                                                            parallel=parallel)
            array[imin:imax, jmin:jmax] = array_sub
            footprint[imin:imax, jmin:jmax] = footprint_sub
    print()
    return array, footprint


def reproject_worker(a):
    imin,imax,jmin,jmax,array_in,wcs_in,wcs_out_sub,shape_out_sub,order=a
    print('worker:',imin,jmin)
    array_sub, footprint_sub = _reproject_celestial(array_in, wcs_in, wcs_out_sub,
                                                    shape_out=shape_out_sub, order=order)
    return imin,imax,jmin,jmax,array_sub,footprint_sub

def reproject_interp_chunk_2d_multi(input_data, output_projection, shape_out=None, hdu_in=0,
                                    order='bilinear', blocks=(1000, 1000), parallel=True):
    """
    For a 2D image, reproject in chunks
    """
    
    array_in, wcs_in = parse_input_data(input_data, hdu_in=hdu_in)
    wcs_out, shape_out = parse_output_projection(output_projection, shape_out=shape_out)

    if isinstance(order, six.string_types):
        order = ORDER[order]

    # Create output arrays
    array = np.zeros(shape_out, dtype=float)
    footprint = np.zeros(shape_out, dtype=float)
    # Arguments for pool
    args=[] 
    print('chunking for pool')
    for imin in range(0, array.shape[0], blocks[0]):
        imax = min(imin + blocks[0], array.shape[0])
        for jmin in range(0, array.shape[1], blocks[1]):
            jmax = min(jmin + blocks[1], array.shape[1])
            shape_out_sub = (imax - imin, jmax - jmin)
            wcs_out_sub = wcs_out.deepcopy()
            wcs_out_sub.wcs.crpix[0] -= jmin
            wcs_out_sub.wcs.crpix[1] -= imin
            args.append((imin,imax,jmin,jmax,array_in,wcs_in,wcs_out_sub,shape_out_sub,order))

    print('pooling')
    p=Pool(16)
    results=p.map(reproject_worker,args)
    for r in results:
        imin,imax,jmin,jmax,array_sub,footprint_sub=r
        array[imin:imax, jmin:jmax] = array_sub
        footprint[imin:imax, jmin:jmax] = footprint_sub

    return array, footprint

if __name__=='__main__':

    reproj=reproject_exact_chunk_2d

    mra=180.0
    mdec=45.0
    cellsize=1.0/3600
    ctype=('RA---SIN','DEC--SIN')

    # pick sizes here
    rxsize=rysize=20000
    xsize=ysize=1024

    rwcs=WCS(naxis=2)
    rwcs.wcs.ctype=ctype
    rwcs.wcs.cdelt=(-cellsize,cellsize)
    rwcs.wcs.crval=[mra,mdec]
    rwcs.wcs.crpix=[1,1]

    rheader=rwcs.to_header()
    rheader['NAXIS']=2
    rheader['NAXIS1']=rxsize
    rheader['NAXIS2']=rysize

    header=rwcs.to_header()
    header['NAXIS']=2
    header['NAXIS1']=xsize
    header['NAXIS2']=ysize

    data=np.random.rand(ysize,xsize)

    hdu=fits.PrimaryHDU(header=header,data=data)
    r,footprint=reproj(hdu, rheader, hdu_in=0, parallel=True)
    fits.PrimaryHDU(header=rheader,data=r).writeto('output.fits',clobber=True)
