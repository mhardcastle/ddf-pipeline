#!/usr/bin/python

# Routine to check quality of LOFAR images
import os,sys
import os.path
from quality_options import options,print_options
import pyfits
import lofar.bdsm as bdsm
from auxcodes import report,run,find_imagenoise,warn,die
import numpy as np
from quality_make_plots import plot_flux_ratios,plot_flux_errors,plot_position_offset

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

def sepn(r1,d1,r2,d2):
    """
    Calculate the separation between 2 sources, RA and Dec must be
    given in radians. Returns the separation in radians
    """
    # NB slalib sla_dsep does this
    # www.starlink.rl.ac.uk/star/docs/sun67.htx/node72.html
    cos_sepn=np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
    sepn = np.arccos(cos_sepn)
    return sepn

def filter_catalog(singlecat,matchedcat,fitsimage,outname,auxcatname,options=None):
    if options is None:
        options = o
    matchedcat = pyfits.open(matchedcat)
    singlecat = pyfits.open(singlecat)
    
    fitsimage = pyfits.open(fitsimage)
    
    fieldra = fitsimage[0].header['CRVAL1']
    fielddec = fitsimage[0].header['CRVAL2']
    fitsimage.close()
    allsources = matchedcat[1].data['Source_id']
    toofarsources = np.array(np.where(sepn(matchedcat[1].data['RA_1']*deg2rad,matchedcat[1].data['DEC_1']*deg2rad,fieldra*deg2rad,fielddec*deg2rad)*rad2deg > 3.0)).flatten()
    print '%s out of %s sources filtered out as over 3.0 deg from centre'%(np.size(toofarsources),len(allsources))
    tooextendedsources_lofar = np.array(np.where(matchedcat[1].data[options['%s_match_majkey1'%auxcatname]] > 10.0)).flatten()
    print '%s out of %s sources filtered out as over 10arcsec in LOFAR'%(np.size(tooextendedsources_lofar),len(allsources))
    tooextendedsources_aux = np.array(np.where(matchedcat[1].data[options['%s_match_majkey2'%auxcatname]] > options['%s_filtersize'%auxcatname])).flatten()
    print '%s out of %s sources filtered out as over %sarcsec in %s'%(np.size(tooextendedsources_aux),len(allsources),options['%s_filtersize'%auxcatname],auxcatname)
    
    groupsize = np.array(np.where(matchedcat[1].data['Groupsize'] > 1.0)).flatten()
    print '%s out of %s sources filtered out as multiple crossmatches'%(np.size(tooextendedsources_aux),len(allsources))
    
    notsingle = np.array([])
    for i in range(0,len(allsources)):
        allseps = sepn(matchedcat[1].data['RA_1'][i]*deg2rad,matchedcat[1].data['DEC_1'][i]*deg2rad,singlecat[1].data['RA']*deg2rad,singlecat[1].data['DEC']*deg2rad)*rad2arcsec
        print min(allseps)
        numclose = np.size(np.where(allseps < 30.0))
        if numclose > 1.0:
            notsingle = np.append(notsingle,i)
    print '%s out of %s sources filtered out as within 20arcsec of another source'%(len(notsingle),len(allsources))
    sourcestoremove = np.concatenate((toofarsources,tooextendedsources_lofar,tooextendedsources_aux,notsingle,groupsize))
    sourcestoremove = np.unique(sourcestoremove)
    print 'In total removing %s out of %s sources'%(np.size(sourcestoremove),np.size(allsources))
    filtereddata = np.delete(matchedcat[1].data,sourcestoremove)
    matchedcat[1].data = filtereddata
    if options['restart'] and os.path.isfile(outname):
        warn('File ' + outname +' already exists, skipping source filtering step')
    else:
        matchedcat.writeto(outname)

def sfind_image(catprefix,pbimage,nonpbimage,sfind_pixel_fraction,options=None):

    if options is None:
        options = o
    f = pyfits.open(nonpbimage)
    imsizex = f[0].header['NAXIS1']
    imsizey = f[0].header['NAXIS2']
    f.close()
    lowerx,upperx = int(((1.0-sfind_pixel_fraction)/2.0)*imsizex),int(((1.0-sfind_pixel_fraction)/2.0)*imsizex + sfind_pixel_fraction*imsizex)
    lowery,uppery = int(((1.0-sfind_pixel_fraction)/2.0)*imsizey),int(((1.0-sfind_pixel_fraction)/2.0)*imsizey + sfind_pixel_fraction*imsizey)

    if options['restart'] and os.path.isfile(catprefix +'.cat.fits'):
        warn('File ' + catprefix +'.cat.fits already exists, skipping source finding step')
    else:
        img = bdsm.process_image(pbimage,adaptive_rms_box='True',advanced_opts='True',detection_image=nonpbimage,thresh_isl=3,thresh_pix=5,blank_limit=1E-7,adaptive_thresh=100,rms_box_bright=[30,10],trim_box=(lowerx,upperx,lowery,uppery))
        img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True')
        img.export_image(outfile=catprefix +'.rms.fits',img_type='rms',img_format='fits',clobber=True)
        img.export_image(outfile=catprefix +'.resid.fits',img_type='gaus_resid',img_format='fits',clobber=True)
        img.export_image(outfile=catprefix +'.pybdsmmask.fits',img_type='island_mask',img_format='fits',clobber=True)
        img.write_catalog(outfile=catprefix +'.cat.reg',catalog_type='srl',format='ds9',correct_proj='True')

def crossmatch_image(lofarcat,auxcatname,options=None):
    if options is None:
        options = o
    print auxcatname
    auxcat = options[auxcatname]
    if options['restart'] and os.path.isfile(lofarcat + '_' + auxcatname + '_match.fits'):
        warn('File ' + lofarcat + '_' + auxcatname + '_match.fits already exists, skipping source matching step')
    else:
        tmpfile = open('tmp_cmatch.sh','w')
        tmpfile.write('java -jar %s %s <<EOF\n'%(o['jstilts'],os.path.realpath(__file__).replace('quality_pipeline.py','stilts-position-match.py')))
        tmpfile.write('%s\n'%lofarcat)
        tmpfile.write('%s\n'%auxcat)
        tmpfile.write('%s\n'%o['%s_matchrad'%auxcatname])
        tmpfile.write('%s_%s_match.fits\n'%(lofarcat,auxcatname))
        tmpfile.write('<<EOF')
        tmpfile.close()
        os.system('chmod +x tmp_cmatch.sh')
        run('./tmp_cmatch.sh',dryrun=options['dryrun'],log=logfilename('%s_%s_match.log'%(lofarcat,auxcatname),options=options),quiet=options['quiet'])
        os.system('rm tmp_cmatch.sh')
        
if __name__=='__main__':
    # Main loop
    if len(sys.argv)<2:
        warn('quality-pipeline.py must be called with a parameter file.\n')
        print_options()
        sys.exit(1)

    o=options(sys.argv[1])
    if o['pbimage'] is None:
        die('pbimage must be specified')
    if o['nonpbimage'] is None:
        die('nonpbimage must be specified')
        
    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])
        
    # pybdsm source finding
    sfind_image(o['catprefix'],o['pbimage'],o['nonpbimage'],o['sfind_pixel_fraction'])

    # matching with catalogs
    crossmatch_image(o['catprefix'] + '.cat.fits','TGSS')
    crossmatch_image(o['catprefix'] + '.cat.fits','FIRST')

    # Filter catalogs (only keep isolated compact sources within 3deg of pointing centre)
    filter_catalog(o['catprefix'] + '.cat.fits','%s.cat.fits_FIRST_match.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_FIRST_match_filtered.fits'%o['catprefix'],'FIRST',options=o)
    filter_catalog(o['catprefix'] + '.cat.fits','%s.cat.fits_TGSS_match.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_TGSS_match_filtered.fits'%o['catprefix'],'TGSS',options=o)

    # Astrometric plots
    plot_position_offset('%s.cat.fits_FIRST_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_FIRST_match_filtered_positions.png'%o['catprefix'],'FIRST',options=o)
    
    # Flux ratio plots (only compact sources)
    plot_flux_ratios('%s.cat.fits_FIRST_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_FIRST_match_filtered_fluxerrors.png'%o['catprefix'],options=o)
    
    # Flux scale comparison plots
    plot_flux_errors('%s.cat.fits_TGSS_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_TGSS_match_filtered_fluxratio.png'%o['catprefix'],'TGSS',options=o)
    
    # Noise estimate
    imagenoise = find_imagenoise(o['nonpbimage'],0.2E-3)
    print 'An estimate of the image noise is %s'%(round(imagenoise*1E3,3))
