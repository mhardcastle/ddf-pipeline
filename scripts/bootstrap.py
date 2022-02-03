#!/usr/bin/env python

# Code to bootstrap the flux density scale using killMS/DDF
# Should run standalone or as part of the pipeline

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div
import os,sys
import os.path
from auxcodes import run,warn,die,MSList
try:
    import bdsf as bdsm
except ImportError:
    from lofar import bdsm
import pyrap.tables as pt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from pipeline import ddf_image,make_external_mask
import shutil
from astropy.io import fits

def logfilename(s):
    if o['logging'] is not None:
        return o['logging']+'/'+s 
    else:
        return None

def run_bootstrap(o):

    # guess colname. This is necesssary because skip_di means there is
    # no DI corrected column. Will break if e.g. sub_square and
    # skip_di are used together.
    
    if o['skip_di']:
        colname=o['colname']
    else:
        colname='DATA_DI_CORRECTED'
        
    if o['mslist'] is None:
        die('MS list must be specified')

    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])

    # check the data supplied
    if o['frequencies'] is None or o['catalogues'] is None:
        die('Frequencies and catalogues options must be specified')

    if "DDF_PIPELINE_CATALOGS" not in list(os.environ.keys()):
        warn("You need to define the environment variable DDF_PIPELINE_CATALOGS where your catalogs are located")
        sys.exit(2)

    o["tgss"]=o["tgss"].replace("$$",os.environ["DDF_PIPELINE_CATALOGS"])
    o["catalogues"]=[l.replace("$$",os.environ["DDF_PIPELINE_CATALOGS"]) for l in o["catalogues"]]
    lCat=o["catalogues"]+[o["tgss"]]
    for fCat in lCat:
        if not os.path.isfile(fCat):
            warn("Catalog %s does not exist"%fCat)
            sys.exit(2)

    cl=len(o['catalogues'])
    if o['names'] is None:
        o['names']=[os.path.basename(x).replace('.fits','') for x in o['catalogues']]
    if o['radii'] is None:
        o['radii']=[10]*cl
    if o['groups'] is None:
        o['groups']=list(range(cl))
    if (len(o['frequencies'])!=cl or len(o['radii'])!=cl or
        len(o['names'])!=cl or len(o['groups'])!=cl):
        die('Names, groups, radii and frequencies entries must be the same length as the catalogue list')

    low_uvrange=[o['image_uvmin'],2.5*206.0/o['low_psf_arcsec']]
    if o['low_imsize'] is not None:
        low_imsize=o['low_imsize'] # allow over-ride
    else:
        low_imsize=old_div(o['imsize']*o['cellsize'],o['low_cell'])

    low_robust=o['low_robust']

    # Clear the shared memory
    run('CleanSHM.py',dryrun=o['dryrun'])

    # We use the individual ms in mslist.
    m=MSList(o['mslist'])
    Uobsid = set(m.obsids)
    
    for obsid in Uobsid:
        
        warn('Running bootstrap for obsid %s' % obsid)

        freqs=[]
        omslist=[]
        for ms,ob,f in zip(m.mss,m.obsids,m.freqs):
            if ob==obsid:
                omslist.append(ms)
                freqs.append(f)

        if len(freqs)<4:
            die('Not enough frequencies to bootstrap. Check your mslist or MS naming scheme')

        # sort to work in frequency order

        freqs,omslist = (list(x) for x in zip(*sorted(zip(freqs, omslist), key=lambda pair: pair[0])))

        for f,ms in zip(freqs,omslist):
            print(ms,f)

        # generate the sorted input mslist
        with open('temp_mslist.txt','w') as f:
            for line in omslist:
                f.write(line+'\n')

        # Clean in cube mode
        # As for the main pipeline, first make a dirty map
        ddf_image('image_bootstrap_'+obsid+'_init','temp_mslist.txt',
                  cleanmask=None,cleanmode='SSD',ddsols='DDS0',
                  applysols=o['apply_sols'][6],majorcycles=0,robust=low_robust,
                  uvrange=low_uvrange,beamsize=o['low_psf_arcsec'],
                  imsize=low_imsize,cellsize=o['low_cell'],
                  options=o,colname=colname,automask=True,
                  automask_threshold=15,smooth=True,cubemode=True,
                  conditional_clearcache=True)
        external_mask='bootstrap_external_mask.fits'
        make_external_mask(external_mask,'image_bootstrap_'+obsid+'_init.dirty.fits',use_tgss=True,clobber=False,cellsize='low_cell',options=o)
        # Deep SSD clean with this external mask and automasking
        ddf_image('image_bootstrap_'+obsid,'temp_mslist.txt',
                  cleanmask=external_mask,reuse_psf=True,reuse_dirty=True,
                  cleanmode='SSD',ddsols='DDS0',applysols=o['apply_sols'][6],
                  majorcycles=5,robust=low_robust,uvrange=low_uvrange,
                  beamsize=o['low_psf_arcsec'],imsize=low_imsize,
                  cellsize=o['low_cell'],options=o,
                  colname=colname,automask=True,
                  automask_threshold=15,smooth=True,cubemode=True,
                  conditional_clearcache=False)

        if os.path.isfile('image_bootstrap_'+obsid+'.cube.int.restored.pybdsm.srl'):
            warn('Source list exists, skipping source extraction')
        else:
            warn('Running PyBDSM, please wait...')
            img=bdsm.process_image('image_bootstrap_'+obsid+'.cube.int.restored.fits',thresh_pix=5,rms_map=True,atrous_do=True,atrous_jmax=2,group_by_isl=True,rms_box=(80,20), adaptive_rms_box=True, adaptive_thresh=80, rms_box_bright=(35,7),mean_map='zero',spectralindex_do=True,specind_maxchan=1,debug=True,kappa_clip=3,flagchan_rms=False,flagchan_snr=False,incl_chan=True,spline_rank=1)
            # Write out in ASCII to work round bug in pybdsm
            img.write_catalog(catalog_type='srl',format='ascii',incl_chan='true')
            img.export_image(img_type='rms',img_format='fits')

        from make_fitting_product import make_catalogue
        import fitting_factors
        import find_outliers

        # generate the fitting product
        if os.path.isfile(obsid+'crossmatch-1.fits'):
            warn('Crossmatch table exists, skipping crossmatch')
        else:
            t = pt.table(omslist[0]+ '/FIELD', readonly=True, ack=False)
            direction = t[0]['PHASE_DIR']
            ra, dec = direction[0]

            if (ra<0):
                ra+=2*np.pi
            ra*=180.0/np.pi
            dec*=180.0/np.pi

            cats=list(zip(o['catalogues'],o['names'],o['groups'],o['radii']))
            make_catalogue('image_bootstrap_'+obsid+'.cube.int.restored.pybdsm.srl',ra,dec,2.5,cats,outnameprefix=obsid)
    
        freqlist=open(obsid+'frequencies.txt','w')
        for n,f in zip(o['names'],o['frequencies']):
            freqlist.write('%f %s_Total_flux %s_E_Total_flux False\n' % (f,n,n))
        for i,f in enumerate(freqs):
            freqlist.write('%f Total_flux_ch%i E_Total_flux_ch%i True\n' % (f,i+1,i+1))
        freqlist.close()

        # Now call the fitting code

        if os.path.isfile(obsid+'crossmatch-results-1.npy'):
            warn('Results 1 exists, skipping first fit')
        else:
            fitting_factors.run_all(1, name=obsid)

        nreject=-1 # avoid error if we fail somewhere
        if os.path.isfile(obsid+'crossmatch-2.fits'):
            warn('Second crossmatch exists, skipping outlier rejection')
        else:
            nreject=find_outliers.run_all(1, name=obsid)
    
        if os.path.isfile(obsid+'crossmatch-results-2.npy'):
            warn('Results 2 exists, skipping second fit')
        else:
          if nreject==0:
              shutil.copyfile(obsid+'crossmatch-results-1.npy',obsid+'crossmatch-results-2.npy')
        if os.path.isfile(obsid+'crossmatch-results-2.npy'):
            warn('Results 2 exists, skipping first fit')
        else:
            fitting_factors.run_all(2, name=obsid)

        # Now apply corrections

        if o['full_mslist'] is None:
            die('Need big mslist to apply corrections')
        if not(o['dryrun']):
            warn('Applying corrections to MS list')
            scale=np.load(obsid+'crossmatch-results-2.npy')[:,0]
            # InterpolatedUS gives us linear interpolation between points
            # and extrapolation outside it
            spl = InterpolatedUnivariateSpline(freqs, scale, k=1)
            
            bigmslist=[s.strip() for s in open(o['full_mslist']).readlines()]
            obigmslist = [ms for ms in bigmslist if obsid in ms]
            
            for ms in obigmslist:
                t = pt.table(ms)
                try:
                    dummy=t.getcoldesc('SCALED_DATA')
                except RuntimeError:
                    dummy=None
                t.close()
                if dummy is not None:
                    warn('Table '+ms+' has already been corrected, skipping')
                else:
                    # in this version we need to scale both the original data and the data in colname
                    t = pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False)
                    frq=t[0]['REF_FREQUENCY']
                    factor=spl(frq)
                    print(frq,factor)
                    t=pt.table(ms,readonly=False)
                    if o['do_wide']:
                        warn('Using DATA_SUB column in bootstrap')
                        desc=t.getcoldesc('DATA_SUB')
                    else:
                        desc=t.getcoldesc(o['colname'])
                    desc['name']='SCALED_DATA'
                    t.addcols(desc)
                    if o['do_wide']:
                        warn('Using DATA_SUB column in bootstrap')
                        d=t.getcol('DATA_SUB')
                    else:
                        d=t.getcol(o['colname'])
                    d*=factor
                    t.putcol('SCALED_DATA',d)
                    try:
                        dummy=t.getcoldesc(colname)
                    except RuntimeError:
                        dummy=None
                    if dummy is not None:
                        desc=t.getcoldesc(colname)
                        newname=colname+'_SCALED'
                        desc['name']=newname
                        t.addcols(desc)
                        d=t.getcol(colname)
                        d*=factor
                        t.putcol(newname,d)

                    t.close()
    if os.path.isfile('image_bootstrap.app.mean.fits'):
        warn('Mean bootstrap image exists, not creating it')
    else:
        warn('Creating mean bootstrap image')
        hdus=[]
        for obsid in Uobsid:
            hdus.append(fits.open('image_bootstrap_'+obsid+'.app.restored.fits'))
        for i in range(1,len(Uobsid)):
            hdus[0][0].data+=hdus[i][0].data
        hdus[0][0].data/=len(Uobsid)
        hdus[0].writeto('image_bootstrap.app.mean.fits')
                
if __name__=='__main__':
    from parset import option_list
    from options import options
    o=options(sys.argv[1:],option_list)
    run_bootstrap(o)
