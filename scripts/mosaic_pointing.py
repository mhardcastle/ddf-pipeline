#!/usr/bin/env python

# intended as a one-stop shop for mosaicing
# contains some of the same arguments as mosaic.py
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
try:
    import bdsf as bdsm
except ImportError:
    import lofar.bdsm as bdsm

restfrq=143.65e6 # should work this out from the FITS headers eventually

def make_header(maxsep,name,ra,dec,cellsize,resolution):
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
    header['RESTFRQ']=restfrq
    header['OBSERVER']='LoTSS'
    header['BUNIT']='JY/BEAM'
    header['BSCALE']=1.0
    header['BZERO']=0
    header['BTYPE']='Intensity'
    header['OBJECT']=name
    return header,himsize

def blank_mosaic(imname,himsize):
    hdu=fits.open(imname)
    x=np.array(range(0,2*himsize))
    xv, yv = np.meshgrid(x, x)
    xv-=himsize
    yv-=himsize
    hdu[0].data[np.sqrt(xv**2.0+yv**2.0)>himsize]=np.nan
    hdu.writeto(imname.replace('.fits','-blanked.fits'), clobber=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Mosaic LoTSS pointings')
    parser.add_argument('--directories', metavar='D', nargs='+',
                        help='directories to search for pipeline output')
    parser.add_argument('--astromap_blank', dest='astromap_blank', default=0.5, help='Acceptable astrometry error in arcsec')
    parser.add_argument('--beamcut', dest='beamcut', default=0.3, help='Beam level to cut at')
    parser.add_argument('--no-check',dest='no_check', action='store_true', help='Do not check for missing images')
    parser.add_argument('--no-highres',dest='no_highres', action='store_true', help='Mosaic low-res images as well')
    parser.add_argument('--do-lowres',dest='do_lowres', action='store_true', help='Mosaic low-res images as well')
    parser.add_argument('--do_scaling',dest='do_scaling',action='store_true',help='Apply scale factor from quality database')
    parser.add_argument('mospointingname', type=str, help='Mosaic central pointing name')
    
    args = parser.parse_args()
    mospointingname = args.mospointingname
    pointingdict = read_pointingfile()

    print 'Now searching for results directories'
    cwd=os.getcwd()

    # find what we need to put in the mosaic
    mosaicpointings,mosseps = find_pointings_to_mosaic(pointingdict,mospointingname)

    maxsep=np.max(mosseps)
    # now find whether we have got these pointings somewhere!
    mosaicdirs=[]
    missingpointing = False
    scales = []
    sdb = SurveysDB()
    for p in mosaicpointings:
        print 'Wanting to put pointing %s in mosaic'%p
        for d in args.directories:
            rd=d+'/'+p
            print rd
            if os.path.isfile(rd+'/image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'):
                mosaicdirs.append(rd)
                try:
                    qualitydict = sdb.get_quality(p)
                    currentdict = sdb.get_field(p)
                    print qualitydict
                    scales.append(qualitydict['scale'])
                except TypeError:
                    missingpointing = True
                    print 'No scaling factor for ',p             
                break
        else:
            print 'Pointing',p,'not found'
            missingpointing = True
        if not missingpointing and (currentdict['status'] != 'Archived' or currentdict['archive_version'] != 4):
            print 'Pointing',p,'not archived with archive_version 4'
            missingpointing = True
    if not(args.no_check) and missingpointing == True:
        sdb.close()
        raise RuntimeError('Failed to find a required pointing')
    sdb.close()
    print 'Mosaicing using directories', mosaicdirs

    # now construct the inputs for make_mosaic

    mos_args=dotdict({'save':True, 'load':True,'exact':False,'use_shifted':True,'find_noise':True})
    mos_args.astromap_blank=args.astromap_blank
    mos_args.beamcut=args.beamcut
    mos_args.directories=mosaicdirs
    if args.do_scaling:
        print 'Applying scales',scales
        mos_args.scale=scales

    if not(args.no_highres):
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],1.5,6.0)

        mos_args.header=header
        print 'Calling make_mosaic'
        #with open('mosaic-header.pickle','w') as f:
        #    pickle.dump(header,f)

        make_mosaic(mos_args)

        print 'Blanking the mosaic...'

        blank_mosaic('mosaic.fits',himsize)

    if args.do_lowres:
        print 'Making the low-resolution mosaic...'
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],4.5,20.0)
        mos_args.header=header
        mos_args.rootname='low'
        mos_args.do_lowres=True
        mos_args.astromap_blank=False # don't bother with low-res map

        #with open('low-mosaic-header.pickle','w') as f:
        #    pickle.dump(header,f)

        make_mosaic(mos_args)

        print 'Blanking the mosaic...'

        blank_mosaic('low-mosaic.fits',himsize)

    print 'Now running PyBDSF to extract sources'
    
    catprefix='mosaic'
    infile='mosaic-blanked.fits'
    if args.no_highres:
        catprefix='low-mosaic'
        infile='low-mosaic-blanked.fits'
    img = bdsm.process_image(infile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(160,50), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0,output_opts=True, output_all=True, atrous_do=True,atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None,frequency=restfrq)
    img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True')
    img.export_image(outfile=catprefix +'.rms.fits',img_type='rms',img_format='fits',clobber=True)
    img.export_image(outfile=catprefix +'.resid.fits',img_type='gaus_resid',img_format='fits',clobber=True)
    img.export_image(outfile=catprefix +'.pybdsmmask.fits',img_type='island_mask',img_format='fits',clobber=True)
    img.write_catalog(outfile=catprefix +'.cat.reg',catalog_type='srl',format='ds9',correct_proj='True')
