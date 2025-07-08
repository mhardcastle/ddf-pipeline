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
from auxcodes import die, get_rms, flatten, convert_regionfile_to_poly, get_rms_map3
import sys
from bad_pix_filter import make_blank_mask
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import numpy as np
import argparse
import pickle
import os.path
import glob
import pyregion
import scipy.ndimage as nd
from copy import deepcopy
import multiprocessing as mp
from queue import Empty
from convolve import do_convolve
from time import sleep

def reproj_inner(q,reproj,hdu,header,shift,direction,ds9region,guard=20):
    print('Direction',direction,'starting')
    r = pyregion.parse(ds9region)
    manualmask = r.get_mask(hdu=hdu)
    # Find bounding box
    yv = np.any(manualmask, axis=1)
    xv = np.any(manualmask, axis=0)
    xmin, xmax = np.where(xv)[0][[0, -1]]
    ymin, ymax = np.where(yv)[0][[0, -1]]
    # Add guard
    xmin-=guard
    if xmin<0: xmin=0
    ymin-=guard
    if ymin<0: ymin=0
    xmax+=guard
    if xmax>hdu.data.shape[1]:
        xmax=hdu.data.shape[1]
    ymax+=guard
    if ymax>hdu.data.shape[0]:
        ymax=hdu.data.shape[0]
    print('Bounding box is',xmin,xmax,ymin,ymax)
    newdata=hdu.data[ymin:ymax,xmin:xmax]
    newheader=deepcopy(hdu.header)
    # adjust the header both to shift and to take account of the subregion
    cellsize=3600*newheader['CDELT2'] # arcsec per pixel
    #newheader['CRPIX1']-=xmin
    #newheader['CRPIX2']-=ymin
    #newheader['CRVAL1']-=shift['RA_offset']/3600.0
    #newheader['CRVAL2']-=shift['DEC_offset']/3600.0
    newheader['CRPIX1']-=xmin+shift['RA_offset']/cellsize
    newheader['CRPIX2']-=ymin-shift['DEC_offset']/cellsize
    shhdu=fits.PrimaryHDU(data=newdata,header=newheader)
    rpm,_=reproj(shhdu,header,hdu_in=0,parallel=False)
    rphdu=fits.PrimaryHDU(header=header,data=rpm)
    newmask = r.get_mask(hdu=rphdu)
    rpm[~newmask]=0
    print('Direction',direction,'returning result to queue')
    q.put(rpm)


def do_reproj_mp(reproj,hdu,header,shift=None,polylist=None,badfacet=None):
    # Wrapper around reproj which handles per-facet reprojection if required
    if shift is None:
        return reproj(hdu,header,hdu_in=0,parallel=False)
    else:
        rpm=None
        q=mp.Queue()
        for direction,ds9region in enumerate(polylist):
            if badfacet and direction in badfacet: continue
            p=mp.Process(target=reproj_inner,args=(q,reproj,hdu,header,shift[direction],direction,ds9region))
            p.start()
        while True:
            n=len(mp.active_children())
            print('In main loop,',n,'active children')
            if n==0: break
            try:
                result=q.get(block=False)
            except Empty:
                print('Tick...')
                sleep(1)
            else:
                if rpm is None:
                    rpm=result
                else:
                    rpm+=result
        return rpm,None  # footprint is not used

def do_reproj(reproj,hdu,header,shift=None,polylist=None,debug=True):
    # Wrapper around reproj which handles per-facet reprojection if required
    if shift is None:
        return reproj(hdu,header,hdu_in=0,parallel=False)
    else:
        rpm=None
        for direction,ds9region in enumerate(polylist):
            print(direction,ds9region)
            r = pyregion.parse(ds9region)
            manualmask = r.get_mask(hdu=hdu)
            print('Convolving to get the new size...')
            manualmask = nd.gaussian_filter(manualmask.astype(float),sigma=3)
            #manualmask=np.binary_dilation(manualmask,structure=np.ones(shape=(5,5)))
            print('Done')
            manualmask = manualmask>0.0 # now a bool array
            yv = np.any(manualmask, axis=1)
            xv = np.any(manualmask, axis=0)
            xmin, xmax = np.where(xv)[0][[0, -1]]
            ymin, ymax = np.where(yv)[0][[0, -1]]
            print('Bounding box is',xmin,xmax,ymin,ymax)
            newdata=hdu.data[ymin:ymax,xmin:xmax]
            newheader=deepcopy(hdu.header)
            # adjust the header both to shift and to take account of the subregion
            newheader['CRPIX1']-=xmin
            newheader['CRPIX2']-=ymin
            #newheader['CRVAL1']+=shift['RA_offset'][direction]/3600.0
            #newheader['CRVAL2']-=shift['DEC_offset'][direction]/3600.0
            newheader['CRVAL1']-=shift['RA_offset'][direction]/3600.0
            newheader['CRVAL2']+=shift['DEC_offset'][direction]/3600.0
            shhdu=fits.PrimaryHDU(data=newdata,header=newheader)
            if debug:
                shhdu.writeto('facet-%i.fits' % direction, overwrite=True)
            if rpm is None:
                rpm,_=reproj(shhdu,header,hdu_in=0,parallel=False)
                rphdu=fits.PrimaryHDU(header=header,data=rpm)
                newmask = r.get_mask(hdu=rphdu)
                rpm[~newmask]=0
            else:
                newmask = r.get_mask(hdu=rphdu)
                rpm+=np.where(newmask,reproj(shhdu,header,hdu_in=0,parallel=False)[0],0)
            if debug:
                rphdu=fits.PrimaryHDU(header=header,data=rpm)
                rphdu.writeto('direction-%i.fits' % direction, overwrite=True)
        return rpm,None  # footprint is not used so safe to return none

def make_mosaic(args):
    if args.find_noise and args.read_noise:
        raise RuntimeError('Cannot both find noise and read it')
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
        reproj=reproject_exact # may not work
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
    elif band is not None and args.use_shifted:
        intname='image_full_ampphase_di_m.NS_Band%i_shift.int.facetRestored.fits' % band
        appname='image_full_ampphase_di_m.NS_Band%i_shift.app.facetRestored.fits' % band
    elif band is not None:
        intname='image_full_ampphase_di_m.NS_Band%i.int.facetRestored.fits' % band
        appname='image_full_ampphase_di_m.NS_Band%i.app.facetRestored.fits' % band
    elif args.use_shifted:
        intname='image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'
        appname='image_full_ampphase_di_m.NS_shift.app.facetRestored.fits'
    else:
        intname='image_full_ampphase_di_m.NS.int.restored.fits'
        appname='image_full_ampphase_di_m.NS.app.restored.fits'

    if args.use_scalefactor:
        sfname=appname.replace('.fits','.scalefactors.fits')
        
    if args.convolve:
        orig_intname=intname
        orig_appname=appname
        intname=intname.replace('.fits','_convolved.fits')
        appname=appname.replace('.fits','_convolved.fits')
        
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
    noisefiles=[]
    name=[]
    shifts=[]
    polylists=[]
    badfacets=[]
    if args.directories is None:
        raise RuntimeError("At least one directory name must be supplied")
    for d in args.directories:
        name.append(d.split('/')[-1])
        infile=d+'/'+intname
        if not os.path.isfile(infile):
            if args.convolve:
                print('Convolved file',infile,'does not exist, making it')
                do_convolve(d+'/'+orig_intname,float(args.convolve),d+'/'+intname)
                do_convolve(d+'/'+orig_appname,float(args.convolve),d+'/'+appname)
            else:
                raise RuntimeError('Expected file',infile,'does not exist')
        hdu=fits.open(infile)

        if args.convolve:
            if hdu[0].header['BMAJ']*3600.0!=float(args.convolve):
                raise RuntimeError('Resolution of convolved file on disk '+infile+' does not match required')

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

        imagefilename=d+'/'+appname
        print('Reading image file',imagefilename)
        if args.do_stokesV:
                tmp = fits.open(imagefilename)
                tmp[0].data[0][0] = tmp[0].data[0][1]
                app.append(flatten(tmp))
        else:
                app.append(flatten(fits.open(imagefilename)))

        if args.use_scalefactor:
            print('Applying scale factor',d+'/'+sfname)
            sfimg=flatten(fits.open(d+'/'+sfname))
            hdus[-1].data*=sfimg.data
            app[-1].data*=sfimg.data
                
        if bth:
            astromaps.append(flatten(fits.open(d+'/astromap.fits')))

        if args.read_noise:
            noisename=d+'/'+appname.replace('.fits','_facetnoise.fits')
            if not os.path.isfile(noisename):
                g=glob.glob(d+'/*.tessel.reg')
                get_rms_map3(d+'/'+appname,g[0],noisename,database=False)
            noisefiles.append(flatten(fits.open(noisename)))
        if args.apply_shift or args.facet_only:
            print('Reading the tessel file')
            g=glob.glob(d+'/*.tessel.reg')
            if len(g)==0:
                raise RuntimeError('apply_shift specified but no tessel file present in '+d)
            else:
                polylists.append(convert_regionfile_to_poly(g[0]))
            
            if args.apply_shift:
                print('Reading the shift file and tessel file')
                t=Table.read(d+'/pslocal-facet_offsets.fits')
                bad=(t['RA_peak']/t['RA_peak_error']<2) | (t['DEC_peak']/t['DEC_peak_error']<2)
                print('Found',np.sum(bad),'bad fits')
                if np.all(bad):
                    print('All bad, zeroing out offsets')
                    t[bad]['RA_offset']=0
                    t[bad]['DEC_offset']=0
                else:
                    print('Replacing bad fits with median shift')
                    ra_median=np.median(t[~bad]['RA_offset'])
                    dec_median=np.median(t[~bad]['DEC_offset'])
                    t[bad]['RA_offset']=ra_median
                    t[bad]['DEC_offset']=dec_median
                shifts.append(t)
            else:
                print('Generating a dummy shift file')
                t=Table([np.zeros(len(polylists[-1])),np.zeros(len(polylists[-1]))],names=('RA_offset','DEC_offset'))
                shifts.append(t)
        else:
            shifts.append(None)
            polylists.append(None)

        if args.use_badfacet:
            badfacetfile=d+'/Badfacets.txt'
            if os.path.isfile(badfacetfile):
                print('Reading the bad facet file')
                lines=open(badfacetfile).readlines()
                bflist=eval(','.join(lines[1].rstrip().split(',')[1:]))
                badfacets.append(bflist)
                print('Adding',len(bflist),'bad facets')
            else:
                badfacets.append([])
        else:
            badfacets.append(None)

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
 # at this point this is the beam factor: we want 1/sigma**2.0, so divide by noise and square
        if args.noise is not None:
            app[i].data/=args.noise[i]
        elif noisefiles:
            app[i].data/=noisefiles[i].data

        app[i].data=app[i].data**2.0


    '''
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
    '''

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

    # ----------------------------------------
    # mosaic main loop
    # ----------------------------------------
    isum=np.zeros([ysize,xsize],dtype=np.float32)
    wsum=np.zeros_like(isum)
    mask=np.zeros_like(isum,dtype=bool)

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
            r, footprint = do_reproj_mp(reproj, hdus[i], header, shift=shifts[i],polylist=polylists[i],badfacet=badfacets[i])
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
            w, footprint = do_reproj_mp(reproj, app[i], header, shift=shifts[i],polylist=polylists[i],badfacet=badfacets[i])
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

        if np.sum(np.isnan(isum))>100:
            # blank small islands -- due to beam errors
            island_mask=make_blank_mask(isum,verbose=True)
            isum[island_mask]=np.nan
        for ch in ('BMAJ', 'BMIN', 'BPA'):
            try:
                header[ch]=hdus[0].header[ch]
            # Exception for Stokes V images which don't have a BMAJ
            except KeyError:
                print('No entry in header for %s and not creating one'%ch)

        header['ORIGIN']='ddf-pipeline '+version()

        hdu = fits.PrimaryHDU(header=header,data=isum)
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
    parser.add_argument('--convolve', default=None, help='Resolution in arcsec to convolve to')
    parser.add_argument('--band', dest='band', default=None, help='Band number to mosaic, leave unset for full-bw image')
    parser.add_argument('--exact', dest='exact', action='store_true', help='Do exact reprojection (slow)')
    parser.add_argument('--save', dest='save', action='store_true', help='Save intermediate images')
    parser.add_argument('--load', dest='load', action='store_true', help='Load existing intermediate images')
    parser.add_argument('--noise', dest='noise', type=float, nargs='+', help='UNSCALED Central noise level for weighting: must match numbers of maps')
    parser.add_argument('--scale', dest='scale', type=float, nargs='+', help='Scale factors by which maps should be multiplied: must match numbers of maps')
    parser.add_argument('--use_shifted', dest='use_shifted', action='store_true', help='Use the shifted images from the pipeline')
    #parser.add_argument('--shift', dest='shift', action='store_true', help='Shift images before mosaicing')
    parser.add_argument('--apply_shift', action='store_true', help='Apply per-facet shift from an offset file')
    parser.add_argument('--facet_only', action='store_true', help='Do not do per-facet shift, but mosaic per facet')
    parser.add_argument('--no_write', dest='no_write', action='store_true', help='Do not write final mosaic')
    parser.add_argument('--find_noise', dest='find_noise', action='store_true', help='Find noise from image')
    parser.add_argument('--read_noise', action='store_true', help='Read noise from a pre-existing per-facet noise file')
    parser.add_argument('--use_badfacet', action='store_true', help='Read a bad facet file')
    parser.add_argument('--use_scalefactor', action='store_true', help='Read a scale factor image')
    parser.add_argument('--do_lowres',dest='do_lowres', action='store_true', help='Mosaic low-res images instead of high-res')
    parser.add_argument('--do_vlow',dest='do_vlow', action='store_true', help='Mosaic vlow images instead of high-res')
    parser.add_argument('--do_wsclean',dest='do_wsclean', action='store_true', help='Mosaic subtracted vlow images instead of high-res')
    parser.add_argument('--do_stokesV',dest='do_stokesV', action='store_true', help='Mosaic stokes V images instead of high-res')
    parser.add_argument('--astromap_blank',dest='astromap_blank', help='')
    parser.add_argument('--load_layout', dest='load_layout', action='store_true', help='Load a previously defined mosaic layout rather than determining from the images.')

    args = parser.parse_args()
    make_mosaic(args)
