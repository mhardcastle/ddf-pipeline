#!/usr/bin/env python

# Routine to check quality of LOFAR images
from __future__ import print_function
from __future__ import division
from past.utils import old_div
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
from auxcodes import report,get_rms,warn,die,sepn,get_centpos
import numpy as np
from crossmatch_utils import match_catalogues,filter_catalogue,select_isolated_sources,bootstrap
from quality_make_plots import plot_flux_ratios,plot_flux_errors,plot_position_offset
from facet_offsets import label_table,RegPoly,plot_offsets
from dr_checker import do_dr_checker
from surveys_db import SurveysDB,get_id,use_database

#Define various angle conversion factors
arcsec2deg=1.0/3600
arcmin2deg=1.0/60
deg2rad=old_div(np.pi,180)
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

        print('Originally',len(matchedcat),'sources')
        matchedcat=filter_catalogue(matchedcat,fieldra,fielddec,3.0)

        print('%i sources after filtering for 3.0 deg from centre' % len(matchedcat))

        matchedcat=matchedcat[matchedcat['DC_Maj']<10.0] # ERROR!

        print('%i sources after filtering for sources over 10arcsec in LOFAR' % len(matchedcat))

        # not implemented yet!
        #tooextendedsources_aux = np.array(np.where(matchedcat[1].data[options['%s_match_majkey2'%auxcatname]] > options['%s_filtersize'%auxcatname])).flatten()
        #print '%s out of %s sources filtered out as over %sarcsec in %s'%(np.size(tooextendedsources_aux),len(allsources),options['%s_filtersize'%auxcatname],auxcatname)

        matchedcat=select_isolated_sources(matchedcat,30.0)
        print('%i sources after filtering for isolated sources in LOFAR' % len(matchedcat))

        matchedcat.write(outname)

def sfind_image(catprefix,pbimage,nonpbimage,sfind_pixel_fraction,options=None):

    if options is None:
        options = o
    f = fits.open(nonpbimage)
    imsizex = f[0].header['NAXIS1']
    imsizey = f[0].header['NAXIS2']
    f.close()
    kwargs={}
    if options['sfind_pixel_fraction']<1.0:
        lowerx,upperx = int(((1.0-sfind_pixel_fraction)/2.0)*imsizex),int(((1.0-sfind_pixel_fraction)/2.0)*imsizex + sfind_pixel_fraction*imsizex)
        lowery,uppery = int(((1.0-sfind_pixel_fraction)/2.0)*imsizey),int(((1.0-sfind_pixel_fraction)/2.0)*imsizey + sfind_pixel_fraction*imsizey)
        kwargs['trim_box']=(lowerx,upperx,lowery,uppery)

    if options['restart'] and os.path.isfile(catprefix +'.cat.fits'):
        warn('File ' + catprefix +'.cat.fits already exists, skipping source finding step')
    else:
        img = bdsm.process_image(pbimage, detection_image=nonpbimage, thresh_isl=4.0, thresh_pix=5.0, rms_box=(160,50), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0,output_opts=True, output_all=True, atrous_do=True,atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, ncores=options['NCPU'], blank_limit=None,**kwargs)
        img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True')
        img.export_image(outfile=catprefix +'.rms.fits',img_type='rms',img_format='fits',clobber=True)
        img.export_image(outfile=catprefix +'.resid.fits',img_type='gaus_resid',img_format='fits',clobber=True)
        img.export_image(outfile=catprefix +'.pybdsmmask.fits',img_type='island_mask',img_format='fits',clobber=True)
        img.write_catalog(outfile=catprefix +'.cat.reg',catalog_type='srl',format='ds9',correct_proj='True')

def crossmatch_image(lofarcat,auxcatname,options=None,catdir='.'):

    if options is None:
        options = o
    auxcat = options[auxcatname]
    crossmatchname=lofarcat + '_' + auxcatname + '_match.fits'
    if options['restart'] and os.path.isfile(crossmatchname):
        warn('File ' + crossmatchname+ ' already exists, skipping source matching step')
        t=Table.read(crossmatchname)
        matches=len(t)
        del(t)
    else:
        t=Table.read(lofarcat)
        tab=Table.read(catdir+'/'+auxcat)
        matches=match_catalogues(t,tab,o[auxcatname+'_matchrad'],auxcatname)
        t=t[~np.isnan(t[auxcatname+'_separation'])]
        t.write(lofarcat+'_'+auxcatname+'_match.fits')
    return matches
        
def do_plot_facet_offsets(t,regfile,savefig=None):
    ''' convenience function to plot offsets '''
    if savefig is not None and os.path.isfile(savefig):
        warn('Figure file %s exists, not re-making it' % savefig)
    else:
        cra,cdec=get_centpos()
        r=RegPoly(regfile,cra,cdec)
        if isinstance(t,str):
            t=Table.read(t)
        if 'Facet' not in t.columns:
            r.add_facet_labels(t)
        plot_offsets(t,r.clist,'red')
        if savefig is not None:
            plt.savefig(savefig)

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
    if o['list'] is not None:
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
        
    if "DDF_PIPELINE_CATALOGS" in list(os.environ.keys()):
        o['catdir']=os.environ["DDF_PIPELINE_CATALOGS"]

    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])
        
    # pybdsm source finding
    sfind_image(o['catprefix'],o['pbimage'],o['nonpbimage'],o['sfind_pixel_fraction'],options=o)

    # facet labels -- do this now for generality
    cra,cdec=get_centpos()
    t=Table.read(o['catprefix'] + '.cat.fits')
    tesselfile=o['catprefix']+'.tessel.reg'
    if 'Facet' not in t.columns:
        t=label_table(t,tesselfile,cra,cdec)
        t.write(o['catprefix'] + '.cat.fits',overwrite=True)

    catsources=len(t)
    
    # matching with catalogs
    removelist=[]
    for cat in o['list']:
        print('Doing catalogue',cat)
        if crossmatch_image(o['catprefix'] + '.cat.fits',cat,catdir=o['catdir'])>10:
            filter_catalog(o['catprefix'] + '.cat.fits',o['catprefix']+'.cat.fits_'+cat+'_match.fits',o['pbimage'],o['catprefix']+'.cat.fits_'+cat+'_match_filtered.fits',cat,options=o)
        else:
            print('Insufficient matches, abandoning catalogue')
            removelist.append(cat)
    for cat in removelist:
        o['list'].remove(cat)

    # Astrometric plots
    if 'FIRST' in o['list']:
        report('Plotting position offsets')
        plot_position_offset('%s.cat.fits_FIRST_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_FIRST_match_filtered_positions.png'%o['catprefix'],'FIRST',options=o)

        t=Table.read(o['catprefix']+'.cat.fits_FIRST_match_filtered.fits')
        bsra=np.percentile(bootstrap(t['FIRST_dRA'],np.mean,10000),(16,84))
        bsdec=np.percentile(bootstrap(t['FIRST_dDEC'],np.mean,10000),(16,84))
        mdra=np.mean(t['FIRST_dRA'])
        mddec=np.mean(t['FIRST_dDEC'])
        print('Mean delta RA is %.3f arcsec (1-sigma %.3f -- %.3f arcsec)' % (mdra,bsra[0],bsra[1]))
        print('Mean delta DEC is %.3f arcsec (1-sigma %.3f -- %.3f arcsec)' % (mddec,bsdec[0],bsdec[1]))
        first_ra=mdra
        first_dec=mddec
        
        report('Plotting per-facet position offsets')
        do_plot_facet_offsets(t,tesselfile,o['catprefix']+'.cat.fits_FIRST_match_filtered_offsets.png')
        t['FIRST_dRA']-=mdra
        t['FIRST_dDEC']-=mddec
        do_plot_facet_offsets(t,tesselfile,o['catprefix']+'.cat.fits_FIRST_match_filtered_offsets_registered.png')

        report('Plotting flux ratios')
        # Flux ratio plots (only compact sources)
        plot_flux_ratios('%s.cat.fits_FIRST_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_FIRST_match_filtered_fluxerrors.png'%o['catprefix'],options=o)
    else:
        first_ra=None
        first_dec=None
    
    report('Plotting flux scale comparison')
    # Flux scale comparison plots
    if 'TGSS' in o['list']:
        plot_flux_errors('%s.cat.fits_TGSS_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_TGSS_match_filtered_fluxratio.png'%o['catprefix'],'TGSS',options=o)
        t=Table.read(o['catprefix']+'.cat.fits_TGSS_match_filtered.fits')
        ratios=old_div(t['Total_flux'],(old_div(t['TGSS_Total_flux'],o['TGSS_fluxfactor'])))
        bsratio=np.percentile(bootstrap(ratios,np.median,10000),(16,84))
        print('Median LOFAR/TGSS ratio is %.3f (1-sigma %.3f -- %.3f)' % (np.median(ratios),bsratio[0],bsratio[1]))
        tgss_scale=np.median(ratios)
    else:
        tgss_scale=None
    if 'NVSS' in o['list']:
        t=Table.read(o['catprefix']+'.cat.fits_NVSS_match_filtered.fits')
        t=t[t['Total_flux']>30e-3]
        ratios=old_div(t['Total_flux'],t['NVSS_Total_flux'])
        bsratio=np.percentile(bootstrap(ratios,np.median,10000),(16,84))
        print('Median LOFAR/NVSS ratio is %.3f (1-sigma %.3f -- %.3f)' % (np.median(ratios),bsratio[0],bsratio[1]))
        nvss_scale=np.median(ratios)
    else:
        nvss_scale=None
    # Noise estimate
    hdu=fits.open(o['pbimage'])
    imagenoise = get_rms(hdu)
    rms=imagenoise*1e6
    print('An estimate of the image noise is %.3f muJy/beam' % rms)
    drs=do_dr_checker(o['catprefix']+'.cat.fits',o['pbimage'],verbose=False,peak=0.4)
    dr=np.median(drs)
    print('Median dynamic range is',dr)

    # fit source counts
    if o['fit_sourcecounts']:
        from fit_sourcecounts import do_fit_sourcecounts
        sc_norm,sc_index,scale=do_fit_sourcecounts(rms=imagenoise)
    else:
        sc_norm=sc_index=scale=None
    
    print(rms,dr,catsources,first_ra,first_dec,tgss_scale,nvss_scale,sc_norm,sc_index,scale)

    if use_database():
        id=get_id()
        with SurveysDB() as sdb:
            result=sdb.create_quality(id)
            result['rms']=float(rms)
            result['dr']=float(dr)
            result['catsources']=int(catsources)
            if first_ra is not None:
                result['first_ra']=float(first_ra)
                result['first_dec']=float(first_dec)
            if tgss_scale is not None:
                result['tgss_scale']=float(tgss_scale)
            if nvss_scale is not None:
                result['nvss_scale']=float(nvss_scale)
            if sc_norm is not None:
                result['sc_norm']=float(sc_norm)
                result['sc_index']=float(sc_index)
                result['sc_scale']=float(scale)
            
            sdb.set_quality(result)
