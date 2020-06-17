#!/usr/bin/env python

# intended as a one-stop shop for mosaicing
# contains some of the same arguments as mosaic.py
from astropy.coordinates import SkyCoord
from astropy.table import Table
import argparse
import os
import glob
from astropy.io import fits
import numpy as np
import pickle
import os,sys
try:
    import bdsf as bdsm
except ImportError:
    import lofar.bdsm as bdsm
from crossmatch_utils import *
import pyregion 
from auxcodes import flatten
from astropy import wcs
from astropy.wcs import WCS

def filter_outside_extract(ds9region,infilename,catalogue):

    hdu=fits.open(infilename)
    hduflat = flatten(hdu)
    map=hdu[0].data
    w = WCS(flatten(hdu).header)

    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)
    inregion = []
    for element in catalogue:
	i,j = w.wcs_world2pix(element['RA'],element['DEC'],0)
	print i,j,manualmask[int(i),int(j)]
	inregion.append(manualmask[int(i),int(j)])
    return catalogue[inregion]

# Run bdsm on the image

parser = argparse.ArgumentParser(description='fitsimage')
parser.add_argument('fitsimage', type=str, help='fitsimage')
parser.add_argument('catalogue', type=str, help='The LoTSS-DR2 catalogue)')
parser.add_argument('regionfile', type=str, help='extractionregion')

args = parser.parse_args()

infile = args.fitsimage
catalogue = args.catalogue
regionfile = args.regionfile

restfrq=143.65e6 # should work this out from the FITS headers eventually

f = fits.open(infile)
ref_ra = f[0].header['CRVAL1']
if ref_ra < 0.0:
	ref_ra  = 360+ref_ra
ref_dec = f[0].header['CRVAL2']

if  not os.path.exists(infile.replace('.fits','cat.srl.fits')):
	img = bdsm.process_image(infile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=True, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None, frequency=restfrq)
	img.write_catalog(outfile=infile.replace('.fits','cat.srl.fits'),catalog_type='srl',format='fits',correct_proj='True')
        img.write_catalog(outfile=infile.replace('.fits','cat.srl.reg'),catalog_type='srl',format='ds9',correct_proj='True')
# LoTSS-DR2 catalogue
# Filter DR2 cat

lotssdr2=Table.read(catalogue)
print 'original length:',len(lotssdr2)
lotssdr2=filter_catalogue(lotssdr2,ref_ra,ref_dec,1.0)
print 'filter around',ref_ra,ref_dec
print 'filter to 1.0 deg:',len(lotssdr2)
lotssdr2=select_isolated_sources(lotssdr2,30)
print 'isolated sources',len(lotssdr2)
lotssdr2=lotssdr2[lotssdr2['Total_flux']/lotssdr2['Isl_rms']>20.0]
print 'snr  more than 20 sources',len(lotssdr2)
lotssdr2=lotssdr2[lotssdr2['S_Code'] == 'S']
print 'S_Code = S sources',len(lotssdr2)
lotssdr2 = lotssdr2[ lotssdr2['Total_flux']/lotssdr2['Peak_flux'] < 1.25 + 3.1*(lotssdr2['Peak_flux']/lotssdr2['Isl_rms'])**-0.53]
print 'Compact sources',len(lotssdr2)

# Filter cutout cat
cutout=Table.read(infile.replace('.fits','cat.srl.fits'))
if cutout['RA'][0] < 0.0:
	cutout['RA'] = 360.0+cutout['RA']
print 'original length:',len(cutout)
cutout=filter_catalogue(cutout,ref_ra,ref_dec,1.0)
print 'filter to 1.0 deg:',len(cutout)
cutout=select_isolated_sources(cutout,60)
print 'isolated sources',len(cutout)
cutout=cutout[cutout['Total_flux']/cutout['Isl_rms']>20.0]
print 'snr  more than 20 sources',len(cutout)
cutout=cutout[cutout['S_Code'] == 'S']
print 'S_Code = S sources',len(cutout)
cutout = cutout[ cutout['Total_flux']/cutout['Peak_flux'] < 1.25 + 3.1*(cutout['Peak_flux']/cutout['Isl_rms'])**-0.53]
print 'Compact sources',len(cutout)


# Simply nearest neighbour match
matched = match_catalogues(lotssdr2,cutout,1,'cutout')
lotssdr2=lotssdr2[~np.isnan(lotssdr2['cutout_separation'])]
print 'After cross match',len(lotssdr2)

# Cut to match extraction region
lotssdr2 = filter_outside_extract(regionfile,infile,lotssdr2)
print len(lotssdr2),'region filtered'

ratios=lotssdr2['Total_flux']/lotssdr2['cutout_Total_flux']/1000.0
lotssdr2.write('matched.fits',overwrite=True)

print 'Median,mean,std',np.median(ratios),np.mean(ratios),np.std(ratios)
print 'Multiply image by ',np.median(ratios)
ds9file = open('matched.ds9.reg','w')
ds9file.write('# Region file format: DS9 version 4.1 \n')
ds9file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
ds9file.write('fk5\n')

for source in lotssdr2:
	ds9file.write('circle(%s,%s,20")\n'%(source['RA'],source['DEC']))
ds9file.close()
