#!/usr/bin/env python

# Routine to check quality of LOFAR images
import matplotlib
matplotlib.use('Agg')

import os,sys
import os.path
from quality_parset import option_list
from options import options,print_options
from astropy.io import fits
from astropy.table import Table
try:
    import bdsf as bdsm
except ImportError:
    import lofar.bdsm as bdsm
from auxcodes import report,run,get_rms,warn,die,sepn
import numpy as np
from crossmatch_utils import match_catalogues,filter_catalogue,select_isolated_sources,bootstrap
from quality_make_plots import plot_flux_ratios,plot_flux_errors,plot_position_offset
from facet_offsets import do_plot_facet_offsets,label_table

#Define various angle conversion factors
arcsec2deg=1.0/3600
arcmin2deg=1.0/60
deg2rad=np.pi/180
deg2arcsec = 1.0/arcsec2deg
rad2deg=180.0/np.pi
arcmin2rad=arcmin2deg*deg2rad
arcsec2rad=arcsec2deg*deg2rad
rad2arcmin=1.0/arcmin2rad
rad2arcsec=1.0/arcsec2rad
steradians2degsquared = (180.0/np.pi)**2.0
degsquared2steradians = 1.0/steradians2degsquared

def logfilename(s,options=None):
    if options is None:
        options=o
    if options['logging'] is not None:
        return options['logging']+'/'+s 
    else:
        return None

def filter_catalog(singlecat,matchedcat,fitsimage,outname,auxcatname,options=None):
    if options is None:
        options = o

    if options['restart'] and os.path.isfile(outname):
        warn('File ' + outname +' already exists, skipping source filtering step')
    else:

        matchedcat = Table.read(matchedcat)
        singlecat = Table.read(singlecat)

        fitsimage = fits.open(fitsimage)

        fieldra = fitsimage[0].header['CRVAL1']
        fielddec = fitsimage[0].header['CRVAL2']
        fitsimage.close()

        print 'Originally',len(matchedcat),'sources'
        matchedcat=filter_catalogue(matchedcat,fieldra,fielddec,3.0)

        print '%i sources after filtering for 3.0 deg from centre' % len(matchedcat)

        matchedcat=matchedcat[matchedcat['DC_Maj']<10.0]

        print '%i sources after filtering for sources over 10arcsec in LOFAR' % len(matchedcat)

        # not implemented yet!
        #tooextendedsources_aux = np.array(np.where(matchedcat[1].data[options['%s_match_majkey2'%auxcatname]] > options['%s_filtersize'%auxcatname])).flatten()
        #print '%s out of %s sources filtered out as over %sarcsec in %s'%(np.size(tooextendedsources_aux),len(allsources),options['%s_filtersize'%auxcatname],auxcatname)

        matchedcat=select_isolated_sources(matchedcat,30.0)
        print '%i sources after filtering for isolated sources in LOFAR' % len(matchedcat)

        matchedcat.write(outname)

def sfind_image(catprefix,pbimage,nonpbimage,sfind_pixel_fraction,options=None):

    if options is None:
        options = o
    f = fits.open(nonpbimage)
    imsizex = f[0].header['NAXIS1']
    imsizey = f[0].header['NAXIS2']
    f.close()
    kwargs={}
    if o['sfind_pixel_fraction']<1.0:
        lowerx,upperx = int(((1.0-sfind_pixel_fraction)/2.0)*imsizex),int(((1.0-sfind_pixel_fraction)/2.0)*imsizex + sfind_pixel_fraction*imsizex)
        lowery,uppery = int(((1.0-sfind_pixel_fraction)/2.0)*imsizey),int(((1.0-sfind_pixel_fraction)/2.0)*imsizey + sfind_pixel_fraction*imsizey)
        kwargs['trim_box']=(lowerx,upperx,lowery,uppery)

    if options['restart'] and os.path.isfile(catprefix +'.cat.fits'):
        warn('File ' + catprefix +'.cat.fits already exists, skipping source finding step')
    else:
        img = bdsm.process_image(pbimage, detection_image=nonpbimage, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0,output_opts=True, output_all=True, atrous_do=True,atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None,**kwargs)
        img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True')
        img.export_image(outfile=catprefix +'.rms.fits',img_type='rms',img_format='fits',clobber=True)
        img.export_image(outfile=catprefix +'.resid.fits',img_type='gaus_resid',img_format='fits',clobber=True)
        img.export_image(outfile=catprefix +'.pybdsmmask.fits',img_type='island_mask',img_format='fits',clobber=True)
        img.write_catalog(outfile=catprefix +'.cat.reg',catalog_type='srl',format='ds9',correct_proj='True')

def crossmatch_image(lofarcat,auxcatname,options=None):

    if options is None:
        options = o
    auxcat = options[auxcatname]
    if options['restart'] and os.path.isfile(lofarcat + '_' + auxcatname + '_match.fits'):
        warn('File ' + lofarcat + '_' + auxcatname + '_match.fits already exists, skipping source matching step')
    else:
        t=Table.read(lofarcat)
        tab=Table.read(auxcat)
        match_catalogues(t,tab,o[auxcatname+'_matchrad'],auxcatname)
        t=t[~np.isnan(t[auxcatname+'_separation'])]
        t.write(lofarcat+'_'+auxcatname+'_match.fits')
        
if __name__=='__main__':
    # Main loop
    if len(sys.argv)<2:
        warn('quality_pipeline.py must be called with at least one parameter file\nor a command-line option list.\nE.g "pipeline.py example.cfg second_example.cfg --solutions-robust=0.1"\nSee below for a complete list of possible options with their default values.')
        print_options(option_list)
        sys.exit(1)

    o=options(sys.argv[1:],option_list)
    if o['pbimage'] is None:
        die('pbimage must be specified')
    if o['nonpbimage'] is None:
        die('nonpbimage must be specified')

    # fix up the new list-type options
    for i,cat in enumerate(o['list']):
        try:
            o[cat]=o['filenames'][i]
        except:
            pass
        try:
            o[cat+'_matchrad']=o['radii'][i]
        except:
            pass
        try:
            o[cat+'_fluxfactor']=o['fluxfactor'][i]
        except:
            pass
        
    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])
        
    # pybdsm source finding
    sfind_image(o['catprefix'],o['pbimage'],o['nonpbimage'],o['sfind_pixel_fraction'])

    # facet labels -- do this now for generality
    t=Table.read(o['catprefix'] + '.cat.fits')
    if 'Facet' not in t.columns:
        t=label_table(t,'image_full_ampphase1m.tessel.reg')
        t.write(o['catprefix'] + '.cat.fits',overwrite=True)

    # matching with catalogs
    for cat in o['list']:
        print 'Doing catalogue',cat
        crossmatch_image(o['catprefix'] + '.cat.fits',cat)
        filter_catalog(o['catprefix'] + '.cat.fits',o['catprefix']+'.cat.fits_'+cat+'_match.fits',o['pbimage'],o['catprefix']+'.cat.fits_'+cat+'_match_filtered.fits',cat,options=o)

    # Filter catalogs (only keep isolated compact sources within 3deg of pointing centre)

    # Astrometric plots
    if 'FIRST' in o['list']:
        report('Plotting position offsets')
        plot_position_offset('%s.cat.fits_FIRST_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_FIRST_match_filtered_positions.png'%o['catprefix'],'FIRST',options=o)

        t=Table.read(o['catprefix']+'.cat.fits_FIRST_match_filtered.fits')
        bsra=np.percentile(bootstrap(t['FIRST_dRA'],np.mean,10000),(16,84))
        bsdec=np.percentile(bootstrap(t['FIRST_dDEC'],np.mean,10000),(16,84))
        mdra=np.mean(t['FIRST_dRA'])
        mddec=np.mean(t['FIRST_dDEC'])
        print 'Mean delta RA is %.3f arcsec (1-sigma %.3f -- %.3f arcsec)' % (mdra,bsra[0],bsra[1])
        print 'Mean delta DEC is %.3f arcsec (1-sigma %.3f -- %.3f arcsec)' % (mddec,bsdec[0],bsdec[1])

        report('Plotting per-facet position offsets')
        do_plot_facet_offsets(t,'image_full_ampphase1m.tessel.reg',o['catprefix']+'.cat.fits_FIRST_match_filtered_offsets.png')
        t['FIRST_dRA']-=mdra
        t['FIRST_dDEC']-=mddec
        do_plot_facet_offsets(t,'image_full_ampphase1m.tessel.reg',o['catprefix']+'.cat.fits_FIRST_match_filtered_offsets_registered.png')

        report('Plotting flux ratios')
        # Flux ratio plots (only compact sources)
        plot_flux_ratios('%s.cat.fits_FIRST_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_FIRST_match_filtered_fluxerrors.png'%o['catprefix'],options=o)
    
    report('Plotting flux scale comparison')
    # Flux scale comparison plots
    if 'TGSS' in o['list']:
        plot_flux_errors('%s.cat.fits_TGSS_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_TGSS_match_filtered_fluxratio.png'%o['catprefix'],'TGSS',options=o)
        t=Table.read(o['catprefix']+'.cat.fits_TGSS_match_filtered.fits')
        ratios=t['Total_flux']/(t['TGSS_Total_flux']/o['TGSS_fluxfactor'])
        bsratio=np.percentile(bootstrap(ratios,np.median,10000),(16,84))
        print 'Median LOFAR/TGSS ratio is %.3f (1-sigma %.3f -- %.3f)' % (np.median(ratios),bsratio[0],bsratio[1])
    if 'NVSS' in o['list']:
        t=Table.read(o['catprefix']+'.cat.fits_NVSS_match_filtered.fits')
        t=t[t['Total_flux']>10e-3]
        ratios=t['Total_flux']/t['NVSS_Total_flux']
        bsratio=np.percentile(bootstrap(ratios,np.median,10000),(16,84))
        print 'Median LOFAR/NVSS ratio is %.3f (1-sigma %.3f -- %.3f)' % (np.median(ratios),bsratio[0],bsratio[1])
    # Noise estimate
    hdu=fits.open(o['pbimage'])
    imagenoise = get_rms(hdu)
    print 'An estimate of the image noise is %.3f muJy/beam' % (imagenoise*1E6)
