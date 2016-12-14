#!/usr/bin/python

# Mosaic final images

# arguments are directories with final images
# we use the .smooth.int.restored.fits and .fluxscale.fits files

from reproject import reproject_interp,reproject_exact
import sys
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import argparse

def flatten(f):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        return f[0].header,f[0].data

    w = WCS(f[0].header)
    wn=WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2
    copy=('EQUINOX','EPOCH','BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
    for k in copy:
        r=f[0].header.get(k)
        if r:
            header[k]=r

    slice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            slice.append(np.s_[:],)
        else:
            slice.append(0)
        
    hdu = fits.PrimaryHDU(header=header,data=f[0].data[slice])
    return hdu

parser = argparse.ArgumentParser(description='Mosaic ddf-pipeline directories')
parser.add_argument('directories', metavar='D', nargs='+',
                    help='directory name')
parser.add_argument('--beamcut', dest='beamcut', default=0.3, help='Beam level to cut at')
parser.add_argument('--exact', dest='exact', action='store_true', help='Do exact reproection (slow)')

args = parser.parse_args()

if args.exact:
    reproj=reproject_exact
else:
    reproj=reproject_interp

threshold=float(args.beamcut)
hdus=[]
app=[]
wcs=[]
for d in args.directories:
    hdus.append(flatten(fits.open(d+'/image_full_ampphase1m.smooth.int.restored.fits')))
    app.append(flatten(fits.open(d+'/image_full_ampphase1m.app.restored.fits')))
    wcs.append(WCS(hdus[-1].header))
    print wcs[-1]

for i in range(len(app)):
    app[i].data/=hdus[i].data
    app[i].data[app[i].data<threshold]=0

ras=np.array([w.wcs.crval[0] for w in wcs])
decs=np.array([w.wcs.crval[1] for w in wcs])

mra=np.mean(ras)
mdec=np.mean(decs)
print 'Will make mosaic at',mra,mdec

# we make a reference WCS and use it to find the extent in pixels
# needed for the combined image

rwcs=WCS(naxis=2)
rwcs.wcs.ctype=wcs[0].wcs.ctype
rwcs.wcs.cdelt=wcs[0].wcs.cdelt
rwcs.wcs.crval=[mra,mdec]
rwcs.wcs.crpix=[1,1]

xmin=0
xmax=0
ymin=0
ymax=0
for a,w in zip(app,wcs):
    ys,xs=np.where(a.data)
    axmin=xs.min()
    aymin=ys.min()
    axmax=xs.max()
    aymax=ys.max()
    del(xs)
    del(ys)
    print 'non-zero',axmin,aymin,axmax,aymax
    for x,y in ((axmin,aymin),(axmax,aymin),(axmin,aymax),(axmax,aymax)):
        ra,dec=[float(f) for f in w.wcs_pix2world(x,y,0)]
        #print ra,dec
        nx,ny=[float (f) for f in rwcs.wcs_world2pix(ra,dec,0)]
        print nx,ny
        if nx<xmin: xmin=nx
        if nx>xmax: xmax=nx
        if ny<ymin: ymin=ny
        if ny>ymax: ymax=ny

print 'co-ord range:', xmin, xmax, ymin, ymax

xsize=int(xmax-xmin)
ysize=int(ymax-ymin)

rwcs.wcs.crpix=[-int(xmin)+1,-int(ymin)+1]
print 'checking:', rwcs.wcs_world2pix(mra,mdec,0)
print rwcs

header=rwcs.to_header()
header['NAXIS']=2
header['NAXIS1']=xsize
header['NAXIS2']=ysize

isum=np.zeros([ysize,xsize])
wsum=np.zeros_like(isum)
mask=np.zeros_like(isum,dtype=np.bool)
print 'now reprojecting'
for i in range(len(hdus)):
    print 'reproject image',i
    r, footprint = reproj(hdus[i], header, hdu_in=0)
    r[np.isnan(r)]=0
    hdu = fits.PrimaryHDU(header=header,data=r)
    hdu.writeto('reproject-%i.fits' % i,clobber=True)
    print '...'
    w, footprint = reproj(app[i], header, hdu_in=0)
    mask|=~np.isnan(w)
    w[np.isnan(w)]=0
    hdu = fits.PrimaryHDU(header=header,data=w)
    hdu.writeto('weight-%i.fits' % i,clobber=True)
    isum+=r*w
    wsum+=w
    print

# mask now contains True where a non-nan region was present in either map
isum/=wsum
isum[~mask]=np.nan
for ch in ('BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER'):
    header[ch]=hdus[0].header[ch]
header['ORIGIN']='ddf-pipeline-mosaic'
header['UNITS']='Jy/beam'

hdu = fits.PrimaryHDU(header=header,data=isum)
hdu.writeto('mosaic.fits',clobber=True)
