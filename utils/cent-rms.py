#!/usr/bin/python

# find the central rms of an image by iteratively removing the outliers

from __future__ import print_function
from __future__ import absolute_import
import numpy as np
from astropy.io import fits
from auxcodes import get_rms

if __name__=='__main__':
    import sys
    for name in sys.argv[1:]:
        hdu=fits.open(name)
        print(name,get_rms(hdu))
