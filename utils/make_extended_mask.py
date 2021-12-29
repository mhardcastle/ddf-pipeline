from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from astropy.io import fits
from astropy.wcs import WCS
from auxcodes import get_rms,get_rms_map2,flatten
import scipy.ndimage as nd
import numpy as np
import pyregion
from scipy.signal import convolve2d

def add_manual_mask(infile,ds9region,outfile):
    hdu=fits.open(infile)
    hduflat = flatten(hdu)

    map=hdu[0].data
    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)

    hdu[0].data=(map.astype(int) | manualmask).astype(np.float32)
    hdu.writeto(outfile,overwrite=True)

def merge_mask(in1,in2,outfile):
    hdu1=fits.open(in1)
    hdu2=fits.open(in2)

    map1=hdu1[0].data.astype(np.int32)
    map2=hdu2[0].data.astype(np.int32)

    hdu1[0].data = (map1 | map2).astype(np.float32)
    hdu1.writeto(outfile,overwrite=True)

def make_extended_mask(infile,fullresfile,rmsthresh=3.0,sizethresh=2500,maxsize=25000,rootname=None,verbose=False,rmsfacet=False,ds9region='image_dirin_SSD_m_c.tessel.reg'):
    ''' infile is the input low-res image, fullresfile is the full-resolution template image, sizethresh the minimum island size in pixels '''

    if rootname is None:
        prefix=''
    else:
        prefix=rootname+'-'

    hdu=fits.open(infile)
    if rmsfacet == False:
        rms=get_rms(hdu,ignore_error=True)
    if rmsfacet == True:
        get_rms_map2(infile,ds9region,prefix+'rms-low.fits')
        hdu2=fits.open(prefix+'rms-low.fits')
        rms=hdu2[0].data[0,0,:]

    det=hdu[0].data[0,0,:]>rmsthresh*rms
    labels, count = nd.label(det)

    print('found',count,'islands')
    #label, counts = np.unique(labels, return_counts=True)
    label=np.unique(labels)
    counts=np.bincount(labels.flatten())

    big=(counts>sizethresh) & (counts<maxsize)
    big_regions=label[big]

    print('Found',len(big_regions),'large islands')
    if verbose: print(counts[big])

    mask=np.zeros_like(det,dtype=int)
    for l in big_regions:
        if l: mask+=l*(labels==l)

    slices=nd.find_objects(mask)
    big_slices=[slices[i-1] for i in big_regions if i]
    kernel = np.ones((3,3))
    mask = convolve2d(mask, kernel, mode='same', fillvalue=0)
    mask = (mask>1)
    w=WCS(hdu[0].header)
    hdu[0].data[0,0]=mask.astype(np.float32)
    hdu.writeto(prefix+'mask-low.fits',clobber=True)

    if fullresfile is not None:

        # regrid all the objects onto the full-res image, if a template is supplied
        hduf=fits.open(fullresfile)
        maskf=np.zeros_like(hduf[0].data[0,0,:,:])
        wf=WCS(hduf[0].header)


        for slice in big_slices:
            yslice=slice[0]
            xslice=slice[1]
            ymin=yslice.start
            ymax=yslice.stop
            xmin=xslice.start
            xmax=xslice.stop
            ppos=[]
            for x in [xmax,xmin]:
                for y in [ymax,ymin]:
                   ppos.append([x,y,0,0])
            worldlim=w.wcs_pix2world(ppos,0)
            wpos=[]
            for ra in [worldlim[:,0].min(),worldlim[:,0].max()]:
                for dec in [worldlim[:,1].min(),worldlim[:,1].max()]:
                    wpos.append([ra,dec,0,0])
            pixlim=wf.wcs_world2pix(wpos,0)
            xminf=int(pixlim[:,0].min())
            xmaxf=int(pixlim[:,0].max())
            yminf=int(pixlim[:,1].min())
            ymaxf=int(pixlim[:,1].max())
            xs=np.arange(xminf,xmaxf)
            ys=np.arange(yminf,ymaxf)
            x,y=np.meshgrid(xs,ys)
            x=x.flatten()
            y=y.flatten()
            pix=np.array([x,y,np.zeros_like(x),np.zeros_like(x)]).T
            world=wf.wcs_pix2world(pix,0)
            opix=w.wcs_world2pix(world,0)
            for xv,yv,op in zip(x,y,opix):
                try:
                    if (xv>=0 and yv>=0) and mask[int(op[1]),int(op[0])]>0:
                        maskf[yv,xv]=1
                except IndexError:
                    # catch wcs mismatches or similar
                    pass

        hduf[0].data[0,0]=maskf.astype(np.float32)
        hduf.writeto(prefix+'mask-high.fits',clobber=True)

if __name__=='__main__':
    import sys
    make_extended_mask(sys.argv[1],sys.argv[2],sizethresh=int(sys.argv[3]),rmsthresh=float(sys.argv[4]),rmsfacet=eval(sys.argv[5]),ds9region=sys.argv[6],rootname='test',verbose=True)
