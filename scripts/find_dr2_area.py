import numpy as np
import os
from astropy_healpix import HEALPix
from astropy import units as u

pos=np.loadtxt(os.environ['DDF_DIR']+'/ddf-pipeline/misc/DR2-pointings.txt',usecols=(1,2))

pos=pos[(pos[:,0]>137) & (pos[:,0]<250)]

print 'Number of pointings is',len(pos)
ra=pos[:,0]
dec=pos[:,1]

print 'RA range is',np.min(ra),np.max(ra)
print 'Dec range is',np.min(dec),np.max(dec)

hp = HEALPix(nside=1024)
print hp.npix,'total healpix pixels on sky'
area=hp.pixel_area.value*3283
print 'area of one healpix is',area,'sq. deg'

plist=[]

for i in range(len(ra)):
    pixels=hp.cone_search_lonlat(ra[i]*u.deg,dec[i]*u.deg,1.85*u.deg)
    plist=plist+list(pixels)

print len(plist),'total pixels'
print len(set(plist)),'total distinct pixels'
print 'Area is',len(set(plist))*area,'sq. deg'


