#!/usr/bin/env python

# Mosaic final images

# arguments are directories with final images

from __future__ import print_function
from __future__ import division
from builtins import zip
from builtins import range
from pipeline_version import version
from reproject import reproject_interp,reproject_exact
from reproj_test import reproject_interp_chunk_2d
from auxcodes import die, get_rms, flatten
import sys
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import numpy as np
import argparse
import pickle
import os.path

def make_mosaic(args):
    if args.scale is not None:
        if len(args.scale) != len(args.directories):
            die('Scales provided must match directories',database=False)

    if args.noise is not None:
        if len(args.noise) != len(args.directories):
            die('Noises provided must match directories',database=False)

    if args.band is not None:
        band=int(args.band)
        print('Doing band %i'%band)
    else:
        band=None
            
    if args.rootname:
        rootname=args.rootname+'-'
    else:
        rootname=''

    if band is not None:
        rootname=('band%i-' % band) + rootname
                
    if args.exact:
        reproj=reproject_exact
    else:
        reproj=reproject_interp_chunk_2d

    if args.do_vlow:
        intname='image_full_vlow_nocut_m.int.restored.fits'
        appname='image_full_vlow_nocut_m.app.restored.fits'
    elif args.do_wsclean:
        intname='WSCLEAN_low-MFS-image-int.fits'
        appname='WSCLEAN_low-MFS-image.fits'
    elif args.do_lowres:
        intname='image_full_low_m.int.restored.fits'
        appname='image_full_low_m.app.restored.fits'
    elif args.do_stokesV:
        intname='image_full_high_stokesV.dirty.corr.fits'
        appname='image_full_high_stokesV.dirty.fits'
    elif band is not None:
        intname='image_full_ampphase_di_m.NS_Band%i_shift.int.facetRestored.fits' % band
        appname='image_full_ampphase_di_m.NS_Band%i_shift.app.facetRestored.fits' % band
    elif args.use_shifted:
        intname='image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'
        appname='image_full_ampphase_di_m.NS_shift.app.facetRestored.fits'
    else:
        intname='image_full_ampphase_di_m.NS.int.restored.fits'
        appname='image_full_ampphase_di_m.NS.app.restored.fits'

    # astromap blanking if required
    bth=None
    try:
        bth=float(args.astromap_blank)
    except:
        pass

    threshold=float(args.beamcut)
    hdus=[]
    app=[]
    astromaps=[]
    wcs=[]
    print('Reading files...')
    noise=[]
    name=[]
    for d in args.directories:
        name.append(d.split('/')[-1])
        hdu=fits.open(d+'/'+intname)

        if args.do_stokesV:
            hdu[0].data[0][0] = hdu[0].data[0][1]

        if args.find_noise:
            print('Estimating noise for', d+'/' + intname)
            if args.do_vlow or args.do_wsclean:
                noise.append(get_rms(hdu,boxsize=500,niter=50))
            elif args.do_lowres:
                noise.append(get_rms(hdu,boxsize=1500))
            else:
                noise.append(get_rms(hdu))
        hdus.append(flatten(hdu))
        if args.do_stokesV:
                tmp = fits.open(d+'/'+appname)
                tmp[0].data[0][0] = tmp[0].data[0][1]
                app.append(flatten(tmp))
        else:
                app.append(flatten(fits.open(d+'/'+appname)))
        if bth:
            astromaps.append(flatten(fits.open(d+'/astromap.fits')))

    if args.find_noise:
        args.noise=noise
        print('Noise values are:')
        for t,n in zip(name,noise):
            print(t,n)

    print('Computing noise/beam factors...')
    for i in range(len(app)):
        np.seterr(divide='ignore')
        app[i].data=np.divide(app[i].data,hdus[i].data)
        app[i].data[app[i].data<threshold]=0
        # at this point this is the beam factor: we want 1/sigma**2.0, so divide by central noise and square
        if args.noise is not None:
                app[i].data/=args.noise[i]

        app[i].data=app[i].data**2.0



    if args.shift:
        print('Finding shifts (NOTE THIS CODE IS OBSOLETE)...')
        # shift according to the FIRST delta ra/dec from quality pipeline
        dras=[]
        ddecs=[]
        for d in args.directories:
            t=Table.read(d+'/image_full_ampphase1m.cat.fits_FIRST_match_filtered.fits')
            dras.append(np.mean(t['FIRST_dRA']))
            ddecs.append(np.mean(t['FIRST_dDEC']))
        print('Applying shifts:',dras,ddecs)
        for i in range(len(app)):
            for hdu in [hdus[i],app[i]]:
                ra=hdu.header['CRVAL1']
                dec=hdu.header['CRVAL2']
                hdu.header['CRVAL1']-=dras[i]/(3600.0*np.cos(np.pi*dec/180.0))
                hdu.header['CRVAL2']-=ddecs[i]/3600.0

    for i in range(len(app)):
        wcs.append(WCS(hdus[i].header))

    # astromap blanking
    if bth:
        print('Blanking using astrometry quality maps with threshold',bth,'arcsec')
        for i in range(len(app)):
            outname=rootname+'astroblank-'+name[i]+'.fits'
            if args.load and os.path.isfile(outname):
                print('Loading previously blanked image')
                hdu=fits.open(outname)
                hdus[i].data=hdu[0].data
            else:
                print('Blanking image',i)
                dmaxy,dmaxx=hdus[i].data.shape
                count=0
                am=astromaps[i]
                awcs=WCS(am.header)
                maxy,maxx=am.data.shape
                for y in range(maxy):
                    for x in range(maxx):
                        value=am.data[y,x]
                        if np.isnan(value):
                            if y<maxy-1:
                                value=am.data[y+1,x]
                        if value>bth:
                            ra,dec=[float(f) for f in awcs.wcs_pix2world(x,y,0)]
                            rx,ry=[int(p) for p in wcs[i].wcs_world2pix(ra,dec,0)]
                            rxp=rx+21 # astromap pix size, with margin
                            ryp=ry+21
                            if rx<0: rx=0
                            if ry<0: ry=0
                            if rxp>dmaxx: rxp=dmaxx
                            if ryp>dmaxy: ryp=dmaxy
                            hdus[i].data[ry:ryp,rx:rxp]=np.nan
                            count+=1
                print('... blanked',count*900.0/3600,'square arcmin')
                if args.save: hdus[i].writeto(outname,overwrite=True)
            app[i].data[np.isnan(hdus[i].data)]=np.nan

    # If the header is directly passed in, use it
    try:
        header=args.header
        xsize=header['NAXIS1']
        ysize=header['NAXIS2']
        print('Mosaic using header passed from calling program')
    except:
        header=None
    if header is None:
        if args.load_layout:
            with open(rootname+'mosaic-header.pickle') as f:
                header=pickle.load(f)
            xsize=header['NAXIS1']
            ysize=header['NAXIS2']
            print('Mosaic using loaded header')
        else:
            print('Creating the mosaic header')
            ras=np.array([w.wcs.crval[0] for w in wcs])
            decs=np.array([w.wcs.crval[1] for w in wcs])

            mra=np.mean(ras)
            mdec=np.mean(decs)
            print('Will make mosaic at',mra,mdec)

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
                print('non-zero',axmin,aymin,axmax,aymax)
                for x,y in ((axmin,aymin),(axmax,aymin),(axmin,aymax),(axmax,aymax)):
                    ra,dec=[float(f) for f in w.wcs_pix2world(x,y,0)]
                    #print ra,dec
                    nx,ny=[float (f) for f in rwcs.wcs_world2pix(ra,dec,0)]
                    print(nx,ny)
                    if nx<xmin: xmin=nx
                    if nx>xmax: xmax=nx
                    if ny<ymin: ymin=ny
                    if ny>ymax: ymax=ny

            print('co-ord range:', xmin, xmax, ymin, ymax)

            xsize=int(xmax-xmin)
            ysize=int(ymax-ymin)

            rwcs.wcs.crpix=[-int(xmin)+1,-int(ymin)+1]
            print('checking:', rwcs.wcs_world2pix(mra,mdec,0))
            print(rwcs)

            header=rwcs.to_header()
            header['NAXIS']=2
            header['NAXIS1']=xsize
            header['NAXIS2']=ysize

            with open(rootname+'mosaic-header.pickle','wb') as f:
                pickle.dump(header,f)

    isum=np.zeros([ysize,xsize])
    wsum=np.zeros_like(isum)
    mask=np.zeros_like(isum,dtype=np.bool)
    print('now making the mosaic')
    for i in range(len(hdus)):
        print('image',i,'(',name[i],')')
        outname=rootname+'reproject-'+name[i]+'.fits'
        if args.load and os.path.exists(outname):
            print('loading...')
            hdu=fits.open(outname)
            r=hdu[0].data
        else:
            print('reprojecting...')
            r, footprint = reproj(hdus[i], header, hdu_in=0, parallel=False)
            r[np.isnan(r)]=0
            hdu = fits.PrimaryHDU(header=header,data=r)
            if args.save: hdu.writeto(outname,overwrite=True)
        print('weights',i,'(',name[i],')')
        outname=rootname+'weight-'+name[i]+'.fits'
        if args.load and os.path.exists(outname):
            print('loading...')
            hdu=fits.open(outname)
            w=hdu[0].data
            mask|=(w>0)
        else:
            print('reprojecting...')
            w, footprint = reproj(app[i], header, hdu_in=0, parallel=False)
            mask|=~np.isnan(w)
            w[np.isnan(w)]=0
            hdu = fits.PrimaryHDU(header=header,data=w)
            if args.save: hdu.writeto(outname,overwrite=True)
        print('add to mosaic...')
        if args.scale is not None:
            print('Applying scale %s to %s'%(args.scale[i],name[i]))
            r = r*args.scale[i]
            w /= args.scale[i]**2.0
        isum+=r*w
        wsum+=w

    if not(args.no_write):
        isum/=wsum
        # mask now contains True where a non-nan region was present in either map
        isum[~mask]=np.nan
        for ch in ('BMAJ', 'BMIN', 'BPA'):
            try:
                header[ch]=hdus[0].header[ch]
            # Exception for Stokes V images where dong have a BMAJ
            except KeyError:
                print('No entry in header for %s and not creating one'%ch)

        header['ORIGIN']='ddf-pipeline '+version()

        hdu = fits.PrimaryHDU(header=header,data=isum)
        if args.do_lowres:
            mosname='low-mosaic.fits'
        else:
            mosname='mosaic.fits'
        hdu.writeto(rootname+mosname,overwrite=True)

        hdu = fits.PrimaryHDU(header=header,data=wsum)
        hdu.writeto(rootname+mosname.replace('.fits','-weights.fits'),overwrite=True)

    else:
        mosname=None

    return rootname+mosname
    
            
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Mosaic ddf-pipeline directories')
    parser.add_argument('--directories', metavar='D', nargs='+',
                        help='directory name')
    parser.add_argument('--rootname', dest='rootname', default='', help='Root name for output files, default uses no prefix')
    parser.add_argument('--beamcut', dest='beamcut', default=0.3, help='Beam level to cut at')
    parser.add_argument('--band', dest='band', default=None, help='Band number to mosaic, leave unset for full-bw image')
    parser.add_argument('--exact', dest='exact', action='store_true', help='Do exact reprojection (slow)')
    parser.add_argument('--save', dest='save', action='store_true', help='Save intermediate images')
    parser.add_argument('--load', dest='load', action='store_true', help='Load existing intermediate images')
    parser.add_argument('--noise', dest='noise', type=float, nargs='+', help='UNSCALED Central noise level for weighting: must match numbers of maps')
    parser.add_argument('--scale', dest='scale', type=float, nargs='+', help='Scale factors by which maps should be multiplied: must match numbers of maps')
    parser.add_argument('--use_shifted', dest='use_shifted', action='store_true', help='Use the shifted images from the pipeline')
    parser.add_argument('--shift', dest='shift', action='store_true', help='Shift images before mosaicing')
    parser.add_argument('--no_write', dest='no_write', action='store_true', help='Do not write final mosaic')
    parser.add_argument('--find_noise', dest='find_noise', action='store_true', help='Find noise from image')
    parser.add_argument('--do_lowres',dest='do_lowres', action='store_true', help='Mosaic low-res images instead of high-res')
    parser.add_argument('--do_vlow',dest='do_vlow', action='store_true', help='Mosaic vlow images instead of high-res')
    parser.add_argument('--do_wsclean',dest='do_wsclean', action='store_true', help='Mosaic subtracted vlow images instead of high-res')
    parser.add_argument('--do_stokesV',dest='do_stokesV', action='store_true', help='Mosaic stokes V images instead of high-res')
    parser.add_argument('--astromap_blank',dest='astromap_blank', help='')
    parser.add_argument('--load_layout', dest='load_layout', action='store_true', help='Load a previously defined mosaic layout rather than determining from the images.')

    args = parser.parse_args()
    make_mosaic(args)
