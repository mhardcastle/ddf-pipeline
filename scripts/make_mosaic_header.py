#!/usr/bin/env python

# make a header for a mosaic using an input image as a template.
# then we can run
# mosaic.py --save --use_shifted --load_layout --find_noise --directories  ...
from __future__ import division
from past.utils import old_div
from auxcodes import flatten
import sys
from astropy.io import fits
import pickle

template=sys.argv[1]
hdus=fits.open(template)

hdu=flatten(hdus)
size=2.5
cellsize=1.5/3600.0
himsize=int(old_div(size,cellsize))
hdu.header['NAXIS1']=2*himsize
hdu.header['NAXIS2']=2*himsize
hdu.header['CRPIX1']=himsize
hdu.header['CRPIX2']=himsize
# fix up headers
hdu.header['TELESCOP']='LOFAR'
hdu.header['RESTFRQ']=143.65e6
hdu.header['OBSERVER']='LoTSS'

with open('mosaic-header.pickle','w') as f:
    pickle.dump(hdu.header,f)
