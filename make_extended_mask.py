from astropy.io import fits
from astropy.wcs import WCS
from auxcodes import get_rms
from auxcodes import flatten
import scipy.ndimage as nd
import numpy as np
import pyregion

def add_manual_mask(infile,ds9region,outfile):
    hdu=fits.open(infile)
    hduflat = flatten(hdu)

    map=hdu[0].data
    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)

    hdu[0].data=(map.astype(int) | manualmask).astype(np.float32)
    hdu.writeto(outfile,clobber=True)


def make_extended_mask(infile,fullresfile,rmsthresh=3.0,sizethresh=2500):

    hdu=fits.open(infile)
    rms=get_rms(hdu)

    det=hdu[0].data[0,0,:]>rmsthresh*rms
    labels, count = nd.label(det)

    print 'found',count,'islands'
    label, counts = np.unique(labels, return_counts=True)

    big=(counts>sizethresh)
    big_regions=label[big]

    print 'Found',len(big_regions)-1,'large islands'

    mask=np.zeros_like(det,dtype=int)
    for l in big_regions:
        if l: mask+=l*(labels==l)

    slices=nd.find_objects(mask)
    big_slices=[slices[i-1] for i in big_regions if i]

    w=WCS(hdu[0].header)
    hdu[0].data=mask
    hdu.writeto('mask-low.fits',overwrite=True)

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
                if mask[int(op[1]),int(op[0])]>0:
                    maskf[yv,xv]=1

        hduf[0].data=maskf
        hduf.writeto('mask-high.fits',overwrite=True)

if __name__=='__main__':
    make_extended_mask('image_low_initial_MSMF.app.restored.fits','image_dirin_MSMF.app.restored.fits')
