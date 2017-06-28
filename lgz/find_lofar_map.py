from astropy.table import Table
from astropy.io import fits
import glob
import numpy as np

t=Table.read('subset_with_id.fits')

# read the maps

g=glob.glob('/data/lofar/mjh/hetdex_v3/mosaics/P*')
outfile=open('lofar-maps.txt','w')

files=[]
ras=[]
decs=[]
for d in g:
    file=d+'/mosaic.fits'
    hdu=fits.open(file)
    ras.append(hdu[0].header['CRVAL1'])
    decs.append(hdu[0].header['CRVAL2'])
    files.append(file)
#    print file,ras[-1],decs[-1]

ras=np.array(ras)
decs=np.array(decs)
for r in t:
    dist=np.cos(decs*np.pi/180.0)*(ras-r['RA'])**2.0 + (decs-r['DEC'])**2.0
    i=np.argmin(dist)
    print >>outfile, r['Source_id'],files[i]

outfile.close()
