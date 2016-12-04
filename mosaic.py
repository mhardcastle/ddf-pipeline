#!/usr/bin/python

# Mosaic final images

# arguments are directories with final images
# we use the .smooth.int.restored.fits and .fluxscale.fits files

from reproject import reproject_interp
import sys
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Mosaic ddf-pipeline directories')
parser.add_argument('directories', metavar='D', nargs='+',
                    help='directory name')
parser.add_argument('--beamcut', dest='beamcut', default=0.3, help='Beam level to cut at')

args = parser.parse_args()

threshold=float(args.beamcut)
hdus=[]
app=[]
wcs=[]
for d in args.directories:
    hdus.append(fits.open(d+'/image_full_ampphase1m.smooth.int.restored.fits'))
    app.append(fits.open(d+'/image_full_ampphase1m.app.restored.fits'))
    wcs.append(WCS(hdus[-1][0].header))
    print wcs[-1]

for i in range(len(app)):
    app[i][0].data/=hdus[i][0].data
    app[i][0].data[app[i][0].data<threshold]=0

ras=np.array([w.wcs.crval[0] for w in wcs])
decs=np.array([w.wcs.crval[1] for w in wcs])

mra=np.mean(ras)
mdec=np.mean(decs)
print 'Will make mosaic at',mra,mdec

# we make a reference WCS and use it to find the extent in pixels
# needed for the combined image

rwcs=WCS(naxis=4)
rwcs.wcs.ctype=wcs[0].wcs.ctype
rwcs.wcs.cdelt=wcs[0].wcs.cdelt
rwcs.wcs.crval=[mra,mdec,wcs[0].wcs.crval[2],wcs[0].wcs.crval[3]]
rwcs.wcs.crpix=[1,1,1,1]

xmin=0
xmax=0
ymin=0
ymax=0
for a,w in zip(app,wcs):
    _,_,ys,xs=np.where(a[0].data)
    axmin=xs.min()
    aymin=ys.min()
    axmax=xs.max()
    aymax=ys.max()
    print 'non-zero',axmin,aymin,axmax,aymax
    for x,y in ((axmin,aymin),(axmax,aymin),(axmin,aymax),(axmax,aymax)):
        ra,dec,_,_ =[float(f) for f in w.wcs_pix2world(x,y,0,0,0)]
        #print ra,dec
        nx,ny,_,_=[float (f) for f in rwcs.wcs_world2pix(ra,dec,0,0,0)]
        print nx,ny
        if nx<xmin: xmin=nx
        if nx>xmax: xmax=nx
        if ny<ymin: ymin=ny
        if ny>ymax: ymax=ny

print 'co-ord range:', xmin, xmax, ymin, ymax

xsize=int(xmax-xmin)
ysize=int(ymax-ymin)

rwcs.wcs.crpix=[-int(xmin)+1,-int(ymin)+1,1,1]
print 'checking:', rwcs.wcs_world2pix(mra,mdec,0,0,0)
print rwcs

header=rwcs.to_header()
header['NAXIS']=4
header['NAXIS1']=xsize
header['NAXIS2']=ysize
header['NAXIS3']=header['NAXIS4']=1

isum=np.zeros([1,1,ysize,xsize])
wsum=np.zeros_like(isum)
mask=np.zeros_like(isum,dtype=np.bool)
print 'now reprojecting'
for i in range(len(hdus)):
    print 'reproject image',i
    r, footprint = reproject_interp(hdus[i], header)
    r[np.isnan(r)]=0
    #hdu = fits.PrimaryHDU(header=header,data=r)
    #hdu.writeto('reproject-%i.fits' % i,clobber=True)
    print '...'
    w, footprint = reproject_interp(app[i], header)
    mask|=~np.isnan(w)
    w[np.isnan(w)]=0
    isum+=r*w
    wsum+=w
    print

# mask now contains True where a non-nan region was present in either map
isum/=wsum
isum[~mask]=np.nan
for ch in ('BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER'):
    header[ch]=hdus[0][0].header[ch]
header['ORIGIN']='ddf-pipeline-mosaic'
header['UNITS']='Jy/beam'

hdu = fits.PrimaryHDU(header=header,data=isum)
hdu.writeto('mosaic.fits',clobber=True)
