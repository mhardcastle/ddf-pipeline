#!/usr/bin/env python

# intended as a one-stop shop for mosaicing
# contains some of the same arguments as mosaic.py
# Tiered pybds approach develoepd by C. Hale whose script was adapted slightly here to run on mosaics rather than app and int images.

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
from shutil import copyfile
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from photutils.aperture import SkyEllipticalAperture
from astropy import units as u
try:
    import bdsf as bdsf
except ImportError:
    import lofar.bdsm as bdsf

restfrq=143.65e6 # should work this out from the FITS headers eventually

def run_tiered_bdsf(imf):
    # Tiered source finding code adapted from Catherine Hale code.

    # --- Create Mean 0 map ---
    folder = os.getcwd()+'/'

    imagef=folder+imf
    meanf=imf[:-5]+'_Mean0.fits'
    label=''

    copyfile(imagef, meanf)

    d=fits.open(meanf, mode='update')
    d[0].data*=0
    d.flush()
    d.close()


    # --- Run initial PyBDSF ---

    restfrq=144000000.0
 
    img = bdsf.process_image(imf, thresh_isl=3.0, thresh_pix=4.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=False, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq)

    img.export_image(outfile=imagef[:-5]+'-Default-'+label+'.rms.fits',img_format='fits', img_type='rms', clobber=True)
    img.export_image(outfile=imagef[:-5]+'-Default-'+label+'.resid.fits',img_format='fits', img_type='gaus_resid', clobber=True)
    img.write_catalog(outfile=imagef[:-5]+'-Default-'+label+'.srl.fits',format='fits', catalog_type='srl', clobber=True)
    img.write_catalog(outfile=imagef[:-5]+'-Default-'+label+'.gaul.fits',format='fits', catalog_type='gaul', clobber=True)
    img.export_image(outfile=imagef[:-5]+'-Default-'+label+'.model.fits',img_format='fits', img_type='gaus_model', clobber=True)

    # --- Add Bright Sources back into field ---

    snr_thresh=150
    source_dat=Table.read(imagef[:-5]+'-Default-'+label+'.srl.fits')
    gaus_dat=Table.read(imagef[:-5]+'-Default-'+label+'.gaul.fits')

    bright_sources=Table(source_dat[source_dat['Peak_flux']/source_dat['Isl_rms']>=snr_thresh], copy=True)
    bright_sources.write(folder+'Bright_Sources_SNR_%i-%s.fits'%(snr_thresh, label), overwrite=True)

    print("Number of bright sources: ", len(bright_sources))

    bright_gaus=Table()
    for i in range(len(bright_sources)):
        bright_gaus=vstack([bright_gaus, gaus_dat[gaus_dat['Source_id']==bright_sources['Source_id'][i]]])

    bright_gaus.write(folder+'Bright_Gaussians_SNR_%i-%s.fits'%(snr_thresh, label), overwrite=True)

    model_im_new_f=imagef[:-5]+'-Default-'+label+'.model-Mask-Bright.fits'
    copyfile(imagef[:-5]+'-Default-'+label+'.model.fits', model_im_new_f)

    wcs_im=WCS(imagef[:-5]+'-Default-'+label+'.model.fits', naxis=2)
    model_im_new=fits.open(model_im_new_f, mode='update')
    model_im_new[0].header = header_N4_to_N2(model_im_new[0].header)
    model_im_new[0].header['WCSAXES']=2
    model_im_new[0].data=model_im_new[0].data[0,0]

    # Mask regions 

    mask_final=np.zeros((np.shape(model_im_new[0].data)[0], np.shape(model_im_new[0].data)[1]))
    for i in range(len(bright_gaus)):
        position = SkyCoord(ra=bright_gaus['RA'][i], dec=bright_gaus['DEC'][i], unit='deg')
        aper = SkyEllipticalAperture(position, 1.5*bright_gaus['Maj'][i]*u.degree, 1.5*bright_gaus['Min'][i]*u.degree, bright_gaus['PA'][i]*u.degree)
        pix_aper=aper.to_pixel(wcs_im)
        m=pix_aper.to_mask()
        out=m.to_image((np.shape(model_im_new[0].data)[0], np.shape(model_im_new[0].data)[1]))
        mask_final+=out

    mask_final[mask_final!=0]=1

    model_im_new[0].data=1-mask_final
    model_im_new.flush()
    model_im_new.close()

    # Make final residual map with bright sources in

    model_im_orig=fits.getdata(imagef[:-5]+'-Default-'+label+'.model.fits')
    model_im_no_brightf= imagef[:-5]+'-Default-'+label+'.model-Bright-Removed.fits'

    copyfile(imagef[:-5]+'-Default-'+label+'.model.fits', model_im_no_brightf)

    model_im_no_bright=fits.open(model_im_no_brightf, mode='update')
    model_im_no_bright[0].data=model_im_new[0].data*model_im_orig
    model_im_no_bright.flush()
    model_im_no_bright.close()

    resid_im_no_brightf=imf[:-5]+'-Default-'+label+'.resid-Bright-Remain.fits'
    copyfile(imagef[:-5]+'-Default-'+label+'.resid.fits', resid_im_no_brightf)

    full_im_orig=fits.getdata(folder+imf)
    resid_im_no_bright=fits.open(resid_im_no_brightf, mode='update')
    resid_im_no_bright[0].data=full_im_orig-model_im_no_bright[0].data
    resid_im_no_bright.flush()
    resid_im_no_bright.close()

    # --- Run PyBDSF on Residual ---

    img_resid = bdsf.process_image(resid_im_no_brightf, thresh_isl=3.0, thresh_pix=4.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=False, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq)

    rms_image_withBright_f=resid_im_no_brightf[:-5]+'.rms.fits'
    img_resid.export_image(outfile=rms_image_withBright_f, img_format='fits', img_type='rms', clobber=True)


    # --- Run PyBDSF - supplying RMS image --- 

    img_supply = bdsf.process_image(imf, thresh_isl=3.0, thresh_pix=4.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=False, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq, rmsmean_map_filename=[meanf, rms_image_withBright_f])

    img_supply.export_image(outfile=imagef[:-5]+'-Default-'+label+'-SupplyMaps.rms.fits',img_format='fits', img_type='rms', clobber=True)
    img_supply.export_image(outfile=imagef[:-5]+'-Default-'+label+'-SupplyMaps.resid.fits',img_format='fits', img_type='gaus_resid', clobber=True)
    img_supply.write_catalog(outfile=imagef[:-5]+'-Default-'+label+'-SupplyMaps.srl.fits',format='fits', catalog_type='srl', clobber=True)
    img_supply.write_catalog(outfile=imagef[:-5]+'-Default-'+label+'-SupplyMaps.gaul.fits',format='fits', catalog_type='gaul', clobber=True)
    img_supply.export_image(outfile=imagef[:-5]+'-Default-'+label+'-SupplyMaps.model.fits',img_format='fits', img_type='gaus_model', clobber=True)
    img_supply.export_image(outfile=imagef[:-5]+'-Default-'+label+'-SupplyMaps.mask.fits',img_format='fits', img_type='island_mask', clobber=True)


    # --- Extract Bright sources from Residual ---

    int_residf=imf[:-5]+'-Default-'+label+'-SupplyMaps.resid.fits'

    img_supply_large = bdsf.process_image(int_residf,thresh_isl=3.0, thresh_pix=10.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=False, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq, rmsmean_map_filename=[meanf, rms_image_withBright_f], flag_maxsize_bm=100)

    img_supply_large.write_catalog(outfile=int_residf+'-Default-'+label+'-SupplyMaps-FlagBeam.gaul.fits',format='fits', catalog_type='gaul', clobber=True)
    img_supply_large.write_catalog(outfile=int_residf+'-Default-'+label+'-SupplyMaps-FlagBeam.srl.fits',format='fits', catalog_type='srl', clobber=True)
    img_supply_large.export_image(outfile=int_residf+'-Default-'+label+'-SupplyMaps-FlagBeam.model.fits',img_format='fits', img_type='gaus_model', clobber=True)
    img_supply_large.export_image(outfile=int_residf+'-Default-'+label+'-SupplyMaps-FlagBeam.mask.fits',img_format='fits', img_type='island_mask', clobber=True)



    # --- Update residual with large sources ---

    resid_ext_emf=imagef+'-Default-'+label+'-SupplyMaps-withFlagBeam.resid.fits'
    model_final=imagef+'-Default-'+label+'-SupplyMaps-withFlagBeam.model.fits'

    model_flagged=int_residf+'-Default-'+label+'-SupplyMaps-FlagBeam.model.fits'
    model_original=imagef[:-5]+'-Default-'+label+'-SupplyMaps.model.fits'

    model_flag_dat=fits.open(model_flagged)
    model_orig_dat=fits.open(model_original)

    copyfile(model_original, model_final)
    copyfile(imagef[:-5]+'-Default-'+label+'-SupplyMaps.resid.fits', resid_ext_emf)

    model_withflag_dat=fits.open(model_final, mode='update')
    resid_withflag_dat=fits.open(resid_ext_emf, mode='update')

    model_withflag_dat[0].data +=model_flag_dat[0].data
    model_withflag_dat.flush()

    resid_withflag_dat[0].data -=model_flag_dat[0].data
    resid_withflag_dat.flush()

    model_withflag_dat.close()
    resid_withflag_dat.close()

    # --- Update Catalogue ---

    dat_srl_flagbeam=Table.read(int_residf+'-Default-'+label+'-SupplyMaps-FlagBeam.srl.fits')
    dat_gaul_flagbeam=Table.read(int_residf+'-Default-'+label+'-SupplyMaps-FlagBeam.gaul.fits')

    dat_srl_orig=Table.read(imagef[:-5]+'-Default-'+label+'-SupplyMaps.srl.fits')
    dat_gaul_orig=Table.read(imagef[:-5]+'-Default-'+label+'-SupplyMaps.gaul.fits')

    max_src_id = np.max(dat_gaul_orig['Source_id'])+1
    max_isl_id = np.max(dat_gaul_orig['Isl_id'])+1
    max_gaus_id = np.max(dat_gaul_orig['Gaus_id'])+1

    dat_srl_flagbeam['Source_id']+=max_src_id+1
    dat_srl_flagbeam['Isl_id']+=max_isl_id+1

    dat_gaul_flagbeam['Source_id']+=max_src_id+1
    dat_gaul_flagbeam['Isl_id']+=max_isl_id+1
    dat_gaul_flagbeam['Gaus_id']+=max_gaus_id+1

    #

    dat_srl_orig['Flag_beam']=np.zeros(len(dat_srl_orig), dtype=int)
    dat_gaul_orig['Flag_beam']=np.zeros(len(dat_gaul_orig), dtype=int)

    dat_srl_flagbeam['Flag_beam']=np.ones(len(dat_srl_flagbeam), dtype=int)
    dat_gaul_flagbeam['Flag_beam']=np.ones(len(dat_gaul_flagbeam), dtype=int)

    #

    dat_srl_final=vstack([dat_srl_orig, dat_srl_flagbeam])
    dat_gaul_final=vstack([dat_gaul_orig, dat_gaul_flagbeam])

    dat_srl_final.write(imagef[:-5]+'-Default-'+label+'-SupplyMaps-withFlagBeam.srl.fits', overwrite=True)
    dat_gaul_final.write(imagef[:-5]+'-Default-'+label+'-SupplyMaps-withFlagBeam.gaul.fits', overwrite=True)


    mask1 = fits.open(imagef[:-5]+'-Default-'+label+'-SupplyMaps.mask.fits')
    mask2 = fits.open(int_residf+'-Default-'+label+'-SupplyMaps-FlagBeam.mask.fits')

    mask1[0].data = mask1[0].data+mask2[0].data
    mask1[0].data[mask1!=0]= 1

    mask1.writeto(imagef[:-5]+'-Default-'+label+'-SupplyMaps-withFlagBeam.mask.fits', overwrite=True)

    # Final products are:
    # mosaic-blanked-Default-P23Hetdex20-SupplyMaps.rms.fits - final rms map (this is the rms map that is supplied to the last source finding step and the rms is not recalcuated in the last step)
    # mosaic-blanked.fits-Default-P23Hetdex20-SupplyMaps-withFlagBeam.resid.fits - final resid map
    # mosaic-blanked.fits-Default-P23Hetdex20-SupplyMaps-withFlagBeam.mask.fits - final mask map
    # mosaic-blanked-Default-P23Hetdex20-SupplyMaps-withFlagBeam.srl.fits - final source catalogue
    # mosaic-blanked-Default-P23Hetdex20-SupplyMaps-withFlagBeam.gaul.fits final guassian catalogue
    return

def header_N4_to_N2(header):
    id=0
    for i in range(len(header.cards)):
        c=header.cards[id][0]
        if c!='' and c!='HISTORY':
            #print(c, c[-1], c[:4])
            if c[-1]=='3' or c[-1]=='4':
                del header[c]
            elif c[:4]=='PC03' or  c[:4]=='PC04':
                del header[c]
            else:
                id+=1
        else:
            id+=1
    return header


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
    parser.add_argument('--do_scaling',dest='do_scaling',action='store_true',help='Apply scale factor from quality database')
    parser.add_argument('--save-header',dest='save_header',action='store_true',help='Save the mosaic header')
    parser.add_argument('mospointingname', type=str, help='Mosaic central pointing name')
    parser.add_argument('--ignorepointings', type=str, default='', help='Pointings to ignore')
    
    args = parser.parse_args()
    mospointingname = args.mospointingname
    pointingdict = read_pointingfile()
    ignorepointings = args.ignorepointings

    if args.do_wsclean:
        fname='WSCLEAN_low-MFS-image-int.fits'
        args.no_highres=True
    elif args.do_vlow:
        fname='image_full_vlow_nocut_m.int.restored.fits'
        args.no_highres=True
    elif args.do_stokesV:
        fname='image_full_high_stokesV.dirty.corr.fits'
        args.no_highres=True
    else:
        fname='image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'
        
    
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
    sdb = SurveysDB()
    for p in mosaicpointings:
        if p in ignorepointings:
            continue
        print('Wanting to put pointing %s in mosaic'%p)
        for d in args.directories:
            rd=d+'/'+p
            print(rd)
            if os.path.isfile(rd+'/'+fname):
                print(rd+'/'+fname,'exists!')
                mosaicdirs.append(rd)
                try:
                    qualitydict = sdb.get_quality(p)
                    currentdict = sdb.get_field(p)
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
    if not(args.no_check) and missingpointing == True:
        sdb.close()
        raise RuntimeError('Failed to find a required pointing')
    sdb.close()
    print('Mosaicing using directories', mosaicdirs)

    # now construct the inputs for make_mosaic

    mos_args=dotdict({'save':True, 'load':True,'exact':False,'use_shifted':True,'find_noise':True})
    mos_args.astromap_blank=args.astromap_blank
    mos_args.beamcut=args.beamcut
    mos_args.directories=mosaicdirs
    mos_args.band=args.band
    if args.do_scaling:
        print('Applying scales',scales)
        mos_args.scale=scales

    if not(args.no_highres):
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],1.5,6.0)

        mos_args.header=header
        print('Calling make_mosaic')
        if args.save_header:
            with open('mosaic-header.pickle','w') as f:
                pickle.dump(header,f)

        mosname=make_mosaic(mos_args)

        print('Blanking the mosaic...')

        blank_mosaic(mosname,himsize)

    if args.do_lowres:
        print('Making the low-resolution mosaic...')
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],4.5,20.0)
        mos_args.header=header
        mos_args.rootname='low'
        mos_args.do_lowres=True
        mos_args.astromap_blank=False # don't bother with low-res map

        if args.save_header:
            with open('low-mosaic-header.pickle','w') as f:
                pickle.dump(header,f)

        make_mosaic(mos_args)

        print('Blanking the mosaic...')

        blank_mosaic('low-mosaic.fits',himsize)

    if args.do_wsclean:
        print('Making the WSCLEAN subtracted mosaic...')
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],15,60.0)
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
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],15,60.0)
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
        header,himsize=make_header(maxsep,mospointingname,pointingdict[mospointingname][1],pointingdict[mospointingname][2],4.5,20.0)
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
	    


        run_tiered_bdsf(infile)

        # OLD (pre Nov 2023) SOURCE FINDING APPROACH
        #img = bdsf.process_image(infile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=True, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None, frequency=restfrq)    
      	#img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True')
        #img.export_image(outfile=catprefix +'.rms.fits',img_type='rms',img_format='fits',clobber=True)
        #img.export_image(outfile=catprefix +'.resid.fits',img_type='gaus_resid',img_format='fits',clobber=True)
        #img.export_image(outfile=catprefix +'.pybdsfmask.fits',img_type='island_mask',img_format='fits',clobber=True)
        #img.write_catalog(outfile=catprefix +'.cat.reg',catalog_type='srl',format='ds9',correct_proj='True')
