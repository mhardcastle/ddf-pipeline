#!/usr/bin/env python

from astropy.convolution import convolve_fft, Gaussian2DKernel
from radio_beam import Beam
from astropy.io import fits
import astropy.units as u
import scipy
import sys
from auxcodes import report
import argparse
from scipy.signal import medfilt2d
import numpy as np

def do_convolve(filename,resolution,outfile=None,scale=True,restore_image=None,intimage=None,appimage=None,appout=None):
    '''Convolve FITS file filename to resolution resolution and save to
    outfile outfile or else to modified version of the input filename.

    If scale is True, assume there is an old beam, map units are
    Jy/beam, find the correct convolution and scale the map
    appropriately.

    If scale is False, don't assume there is an old beam, map units are Jy.

    '''
    if appimage:
        # our input image is a true flux image and we are going to
        # make an apparent flux convolved image using intimage and appimage as a
        # template. Compute the beam image and smooth it to remove
        # outliers.
        report('Reading files to construct beam')
        hdu_int=fits.open(intimage,memmap=False)
        hdu_app=fits.open(appimage,memmap=False)
        report('Making smooth beam image')
        with np.errstate(divide='ignore'):
            beam=hdu_app[0].data[0,0]/hdu_int[0].data[0,0]
        beam=medfilt2d(beam)
        hdu_int.close()
        hdu_app.close()
    report('Reading files')
    hdu=fits.open(filename,memmap=False)
    pixelscale=3600.0*hdu[0].header['CDELT2']
    
    new_b=Beam(resolution*u.arcsec)
    if scale:
        old_b=Beam.from_fits_header(hdu[0].header)
        rr=new_b.sr/old_b.sr
        print('Beam area ratio is %.2f' % rr)
        dcb=new_b.deconvolve(old_b)
        gauss_kern = dcb.as_kernel(pixelscale*u.arcsec)
    else:
        gauss_kern = new_b.as_kernel(pixelscale*u.arcsec)
        rr=new_b.sr/(pixelscale*u.arcsec)**2 # beam area in pixels
        rr=rr.to(u.dimensionless_unscaled)
        print('Beam area in pixels is %.2f' % rr)
    image=hdu[0].data[0,0]
    report('FFT')
    fft_mp = lambda a: scipy.fft.fftn(a, workers=-1)
    ifft_mp = lambda a: scipy.fft.ifftn(a, workers=-1)

    smoothed_data_gauss = convolve_fft(image, gauss_kern, allow_huge=True, fftn=fft_mp, ifftn=ifft_mp)
    report('Scaling convolved image')
    hdu[0].data[0,0]=smoothed_data_gauss*rr
    hdu[0].header.update(new_b.to_header_keywords())
    if restore_image is not None:
        report('Restoring to residual')
        hdu2=fits.open(restore_image)
        hdu[0].data+=hdu2[0].data
    report('Write to disk')
    if outfile is None:
        outfile=filename.replace('.fits','_convolved.fits')
    hdu.writeto(outfile,overwrite=True)
    if appimage:
        hdu[0].data[0,0]*=beam
        if appout is not None:
            outfile=appout
        else:
            outfile=outfile.replace('.int.','.app.')
        hdu.writeto(outfile,overwrite=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Convolve and optionally restore an image')
    parser.add_argument('--model',action='store_true',help='Convolving a model image (units Jy)')
    parser.add_argument('--restore',help='Image to restore')
    parser.add_argument('--resolution',help='Resolution to convolve to in arcsec')
    parser.add_argument('--outfile',help='Output file')
    parser.add_argument('--appimage',help='Apparent flux image to generate beam')
    parser.add_argument('--intimage',help='True flux image to generate beam')
    parser.add_argument('--appout',help='Apparent image output file')
    
    parser.add_argument('filename')
    args = parser.parse_args()
    filename=args.filename
    resolution=float(args.resolution)
    do_convolve(filename,resolution,scale=not args.model,restore_image=args.restore,outfile=args.outfile,appimage=args.appimage,intimage=args.intimage,appout=args.appout)
    
