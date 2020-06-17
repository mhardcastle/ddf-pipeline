#!/usr/bin/env python
# Call with the ILT name or other object name plus a region size in degrees

from __future__ import print_function
from __future__ import absolute_import
from find_pos import find_pos
import sys
import os
from astropy.coordinates import SkyCoord,get_icrs_coordinates
from make_subim import extract_and_save
import astropy.units as u

objname=sys.argv[1]
try:
    imsize=float(sys.argv[2])
except:
    imsize=0.5

if 'ILTJ' in objname:
    s=objname[4:]
    coord=s[0:2]+':'+s[2:4]+':'+s[4:9]+' '+s[9:12]+':'+s[12:14]+':'+s[14:]
    sc = SkyCoord(coord,unit=(u.hourangle,u.deg))
    ra=sc.ra.value
    dec=sc.dec.value
    print('Parsed coordinates to ra=%f, dec=%f' % (ra,dec))
elif 'pos' in objname:
    ra=float(sys.argv[3])
    dec=float(sys.argv[4])
    objname='pos-%f-%f' % (ra,dec)
else:
    c=get_icrs_coordinates(objname)
    ra=float(c.ra.degree)
    dec=float(c.dec.degree)
    print('Coordinate lookup gives ra=%f, dec=%f' % (ra,dec))

field=find_pos(ra,dec,offset=3.5)

if field is None:
    print('Archived image does not exist')
else:
    wd='/data/lofar/DR2/fields/'+field
    if not os.path.isdir(wd):
        print('DR2 directory does not exist')
    else:
        print('Extracting FULL total intensity cutout')
        extract_and_save(wd+'/image_full_ampphase_di_m.NS_shift.int.facetRestored.fits',ra,dec,imsize,outname=objname+'_I.fits')
        print('Extracting LOW total intensity cutout')
        extract_and_save(wd+'/image_full_low_m.int.restored.fits',ra,dec,imsize,outname=objname+'_I_low.fits')
        
