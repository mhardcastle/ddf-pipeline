#!/usr/bin/env python
from __future__ import print_function
import sys
from auxcodes import get_rms_map3
import os
import glob
from convolve import do_convolve
from astropy.io import fits

d=sys.argv[1]

if not os.path.isdir(d):
    raise RuntimeError('Directory '+d+' does not exist')

os.chdir(d)

g=glob.glob('*.tessel.reg')
for appname in ['image_full_ampphase_di_m.NS.app.restored.fits','image_full_low_m.app.restored.fits']:
    if not os.path.isfile(appname):
        raise RuntimeError('File is missing: '+appname)
    scalefactor=appname.replace('.fits','.scalefactors.fits')
    if os.path.isfile(scalefactor):
        scalefactor_applied=appname.replace('.fits','.scaled.fits')
        if not os.path.isfile(scalefactor_applied):
            hdu_orig=fits.open(appname)
            hdu_scalefactor=fits.open(scalefactor)
            hdu_orig[0].data*=hdu_scalefactor[0].data
            hdu_orig.writeto(scalefactor_applied,overwrite=True)
            appname=scalefactor_applied
        
    noisename=appname.replace('.fits','_facetnoise.fits')
    if not os.path.isfile(noisename):
        print('Making facet noise map from',appname)
        get_rms_map3(appname,g[0],noisename,database=False)

decl=float(d.split('+')[1])
if decl<=14:
    print('Low-declination field!')
    filename='image_full_ampphase_di_m.NS.int.restored.fits'
    convname=filename.replace('.fits','_convolved.fits')
    if not os.path.isfile(convname):
        print('Convolving to low-dec resolution')
        do_convolve('image_full_ampphase_di_m.NS.int.model.fits',9.0,outfile=convname,scale=False,restore_image='image_full_ampphase_di_m.NS.int.residual.fits',intimage=filename,appimage=filename.replace('.int','.app'))

    appname='image_full_ampphase_di_m.NS.app.restored_convolved.fits'
    scalefactor=appname.replace('.fits','.scalefactors.fits')
    if os.path.isfile(scalefactor):
        scalefactor_applied=appname.replace('.fits','.scaled.fits')
        if not os.path.isfile(scalefactor_applied):
            hdu_orig=fits.open(appname)
            hdu_scalefactor=fits.open(scalefactor)
            hdu_orig[0].data*=hdu_scalefactor[0].data
            hdu_orig.writeto(scalefactor_applied,overwrite=True)
            appname=scalefactor_applied

    noisename=appname.replace('.fits','_facetnoise.fits')
    if not os.path.isfile(noisename):
        print('Making facet noise map from',appname)
        get_rms_map3(appname,g[0],noisename,database=False)
