#!/usr/bin/env python

# intended as a one-stop shop for mosaicing
# contains some of the same arguments as mosaic.py

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import range
import argparse
from find_mosaic_pointings import read_pointingfile, find_pointings_to_mosaic
import os
from auxcodes import getpos,getposim,dotdict
import glob
from mosaic import make_mosaic
from astropy.io import fits
import numpy as np
from surveys_db import SurveysDB
import pickle
import os,sys
from sourcefind import run_tiered_bdsf as run_bdsf

def make_header(maxsep,name,ra,dec,cellsize,resolution,history=None):
    # construct template FITS header
    header=fits.Header()
    size=(maxsep/2.0)*1.15
    cellsize/=3600.0
    himsize=int(size/cellsize)
    header['SIMPLE']=True
    header['BITPIX']=-32
    header['NAXIS']=2
    header['WCSAXES']=2
    header['NAXIS1']=2*himsize
    header['NAXIS2']=2*himsize
    header['CTYPE1']='RA---SIN'
    header['CTYPE2']='DEC--SIN'
    header['CUNIT1']='deg'
    header['CUNIT2']='deg'
    header['CRPIX1']=himsize
    header['CRPIX2']=himsize
    header['CRVAL1']=ra
    header['CRVAL2']=dec
    header['CDELT1']=-cellsize
    header['CDELT2']=cellsize
    header['RADESYS']='ICRS'
    header['EQUINOX']=2000.0
    header['LONPOLE']=180.0
    header['LATPOLE']=header['CRVAL2']
    header['BMAJ']=resolution/3600.0
    header['BMIN']=resolution/3600.0
    header['BPA']=0
    header['TELESCOP']='LOFAR'
    header['RESTFRQ']=144e6
    header['OBSERVER']='LoTSS'
    header['BUNIT']='JY/BEAM'
    header['BSCALE']=1.0
    header['BZERO']=0
    header['BTYPE']='Intensity'
    header['OBJECT']=name
    if history is not None:
        for h in history:
            header.add_history(h)
    return header,himsize

def blank_mosaic(imname,himsize):
    hdu=fits.open(imname)
    x=np.array(list(range(0,2*himsize)))
    xv, yv = np.meshgrid(x, x)
    xv-=himsize
    yv-=himsize
    hdu[0].data[np.sqrt(xv**2.0+yv**2.0)>himsize]=np.nan
    hdu.writeto(imname.replace('.fits','-blanked.fits'), overwrite=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Mosaic LoTSS pointings')
    parser.add_argument('--directories', metavar='D', nargs='+',
                        help='directories to search for pipeline output')
    parser.add_argument('--astromap_blank', dest='astromap_blank', default=0.5, help='Acceptable astrometry error in arcsec')
    parser.add_argument('--beamcut', dest='beamcut', default=0.3, help='Beam level to cut at')
    parser.add_argument('--band', dest='band', default=None, help='Band to mosaic for multi-band data')
    parser.add_argument('--no-check',dest='no_check', action='store_true', help='Do not check for missing images')
    parser.add_argument('--no-highres',dest='no_highres', action='store_true', help='Skip the high-resolution images')
    parser.add_argument('--no-bdsf',dest='no_bdsf', action='store_true', help='Skip the source extraction')
    parser.add_argument('--do-lowres',dest='do_lowres', action='store_true', help='Mosaic low-res images as well')
    parser.add_argument('--do-wsclean',dest='do_wsclean', action='store_true', help='Mosaic subtracted WSCLEAN images')
    parser.add_argument('--do-vlow',dest='do_vlow', action='store_true', help='Mosaic vlow images as well')
    parser.add_argument('--do-stokesV',dest='do_stokesV', action='store_true', help='Mosaic stokes V images as well')
    parser.add_argument('--record-checksum',dest='record_checksum', action='store_true', help='Store checksum values in history')
    parser.add_argument('--do_scaling',dest='do_scaling',action='store_true',help='Apply scale factor from quality database')
    parser.add_argument('--save-header',dest='save_header',action='store_true',help='Save the mosaic header')
    parser.add_argument('--apply-shift',dest='apply_shift',action='store_true',help='Apply per-facet shifts from an offset file')
    parser.add_argument('mospointingname', type=str, help='Mosaic central pointing name')
    parser.add_argument('--ignorepointings', type=str, default='', help='Pointings to ignore')
    parser.add_argument('--ignore_field', type=str, default='', help='Ignore pointings without this DB field set positive')
    
    args = parser.parse_args()
    mospointingname = args.mospointingname
    pointingdict = read_pointingfile()
    ignorepointings = args.ignorepointings
    check_convolve=False # check resolution if True
    use_badfacet=False
    with SurveysDB(readonly=True) as sdb:
        field_dict=sdb.get_field(mospointingname)

    if args.do_wsclean:
        fname='WSCLEAN_low-MFS-image-int.fits'
        args.no_highres=True
    elif args.do_vlow:
        fname='image_full_vlow_nocut_m.int.restored.fits'
        args.no_highres=True
    elif args.do_stokesV:
        fname='image_full_high_stokesV.dirty.corr.fits'
        check_convolve=True
        args.no_highres=True
        use_badfacet=True
    elif args.apply_shift:
        fname='image_full_ampphase_di_m.NS.app.restored.fits'
        check_convolve=True
        use_badfacet=True
    else:
        fname='image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'
        check_convolve=True
        
    print('Now searching for results directories')
    cwd=os.getcwd()

    # find what we need to put in the mosaic
    mosaicpointings,mosseps = find_pointings_to_mosaic(pointingdict,mospointingname)
    if ignorepointings != '':
        ignorepointings = ignorepointings.split(',')

    maxsep=np.max(mosseps)
    # now find whether we have got these pointings somewhere!
    mosaicdirs=[]
    missingpointing = False
    scales = []
    resolutions = []
    if args.record_checksum:
        checksums=[]
    else:
        checksums=None
    sdb = SurveysDB(readonly=True)
    for p in mosaicpointings:
        if p in ignorepointings:
            continue
        currentdict = sdb.get_field(p)
        if args.ignore_field and not currentdict[args.ignore_field]:
            continue
        print('Wanting to put pointing %s in mosaic'%p)
        for d in args.directories:
            rd=d+'/'+p
            print(rd)
            if os.path.isfile(rd+'/'+fname):
                print(rd+'/'+fname,'exists!')
                mosaicdirs.append(rd)
                if args.record_checksum:
                    with open(rd+'/checksums.txt') as cs:
                        csl=[l.rstrip() for l in cs.readlines()]
                    for l in csl:
                        checksums.append(p+','+l)
                # check resolution 
                if check_convolve:
                    hdu=fits.open(rd+'/'+fname)
                    resolutions.append((hdu[0].header['BMAJ'],hdu[0].header['BMIN']))
                    print('For',fname,'appending resolution',np.array(resolutions[-1])*3600)
                    hdu.close()
                # check quality
                try:
                    qualitydict = sdb.get_quality(p)
                    print(qualitydict)
                    #scale=qualitydict['scale']
                    scale= 1.0/(qualitydict['nvss_scale']/5.9124)
                    if scale is None:
                        print('Missing scaling factor for',p)
                        missingpointing=True
                        scale=1.0
                    scales.append(scale)

                except TypeError:
                    missingpointing = True
                    print('No scaling factor for ',p)
                    scales.append(1.0)
                break
        else:
            print('Pointing',p,'not found')
            missingpointing = True
            if not missingpointing and (currentdict['status'] != 'Archived' or currentdict['archive_version'] != 4):
                print('Pointing',p,'not archived with archive_version 4')
                missingpointing = True
    sdb.close()
    if not(args.no_check) and missingpointing == True:
        raise RuntimeError('Failed to find a required pointing')

    print('Mosaicing using directories', mosaicdirs)

    # now construct the inputs for make_mosaic
    high_resolution=6.0

    if check_convolve:
        # We should convolve to a fixed resolution of 9 x 9 arcsec if:
        # All images do not have the same resolution OR
        # Any image has a non-circular beam.
        # This should catch all cases
        resolutions=np.array(resolutions)
        print('Resolution checking')
        print(resolutions[:,0]-np.mean(resolutions[:,0]))
        different=np.any(np.abs(resolutions[:,0]-np.mean(resolutions[:,0]))>1e-5)
        non_circ=not np.all(np.equal(resolutions[:,0],resolutions[:,1]))
        print('Resolutions are different:',different)
        print('Some beams are non-circular:',non_circ)

    mos_args=dotdict({'save':True, 'load':True,'exact':False})
    if args.apply_shift:
        if np.abs(field_dict['gal_b'])<=10:
            mos_args.facet_only=True
            mos_args.apply_shift=False
        else:
            mos_args.facet_only=False
            mos_args.apply_shift=True
        mos_args.read_noise=True
        mos_args.astromap_blank=None
    else:
        mos_args.use_shifted=True
        mos_args.find_noise=True
        mos_args.astromap_blank=args.astromap_blank
    if check_convolve and (different or non_circ):
        mos_args.convolve=9.0
        high_resolution=9.0
    mos_args.beamcut=args.beamcut
    mos_args.directories=mosaicdirs
    mos_args.band=args.band
    if use_badfacet:
        mos_args.use_badfacet=True
    if args.do_scaling:
        print('Applying scales',scales)
        mos_args.scale=scales

    if not(args.no_highres):
        print('Making the high-resolution mosaic')
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],1.5,high_resolution,history=checksums)

        mos_args.header=header
        print('Calling make_mosaic')
        if args.save_header:
            with open('mosaic-header.pickle','wb') as f:
                pickle.dump(header,f)

        mosname=make_mosaic(mos_args)

        print('Blanking the mosaic...')

        blank_mosaic(mosname,himsize)

    if args.do_lowres:
        print('Making the low-resolution mosaic...')
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],4.5,20.0,history=checksums)
        mos_args.header=header
        mos_args.rootname='low'
        mos_args.do_lowres=True
        mos_args.astromap_blank=False # don't bother with low-res map
        mos_args.convolve=None

        if args.save_header:
            with open('low-mosaic-header.pickle','wb') as f:
                pickle.dump(header,f)

        make_mosaic(mos_args)

        print('Blanking the mosaic...')

        blank_mosaic('low-mosaic.fits',himsize)

    if args.do_wsclean:
        print('Making the WSCLEAN subtracted mosaic...')
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],15,60.0,history=checksums)
        mos_args.header=header
        mos_args.rootname='vlow'
        mos_args.do_wsclean=True
        mos_args.astromap_blank=False # don't bother with low-res map

        if args.save_header:
            with open('vlow-mosaic-header.pickle','w') as f:
                pickle.dump(header,f)

        make_mosaic(mos_args)

        print('Blanking the mosaic...')

        blank_mosaic('vlow-mosaic.fits',himsize)

    if args.do_vlow:
        print('Making the very low-resolution mosaic...')
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],15,60.0,history=checksums)
        mos_args.header=header
        mos_args.rootname='vlow'
        mos_args.do_vlow=True
        mos_args.astromap_blank=False # don't bother with low-res map

        if args.save_header:
            with open('vlow-mosaic-header.pickle','w') as f:
                pickle.dump(header,f)

        make_mosaic(mos_args)

        print('Blanking the mosaic...')

        blank_mosaic('vlow-mosaic.fits',himsize)


    if args.do_stokesV:
        print('Making the stokes V mosaic...')
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],4.5,20.0,history=checksums)
        mos_args.header=header
        mos_args.rootname='stokesV'
        mos_args.do_stokesV=True
        mos_args.astromap_blank=False # don't bother with low-res map

        if args.save_header:
            with open('stokesV-mosaic-header.pickle','w') as f:
                pickle.dump(header,f)

        make_mosaic(mos_args)

        print('Blanking the mosaic...')

        blank_mosaic('stokesV-mosaic.fits',himsize)


    if args.band is None and not args.no_bdsf:
        print('Now running PyBDSF to extract sources')

        catprefix='mosaic'
        infile='mosaic-blanked.fits'
        if args.no_highres:
            catprefix='low-mosaic'
            infile='low-mosaic-blanked.fits'
	    
        run_bdsf(infile,infile,catprefix=catprefix)
