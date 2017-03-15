#!/usr/bin/env python

# Code to bootstrap the flux density scale using killMS/DDF
# Should run standalone or as part of the pipeline

import os,sys
import os.path
from auxcodes import run,warn,die
from lofar import bdsm
import pyrap.tables as pt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from pipeline import ddf_image
import shutil

def logfilename(s):
    if o['logging'] is not None:
        return o['logging']+'/'+s 
    else:
        return None

def run_bootstrap(o):
    
    if o['mslist'] is None:
        die('MS list must be specified')

    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])

    # check the data supplied
    if o['frequencies'] is None or o['catalogues'] is None:
        die('Frequencies and catalogues options must be specified')

    cl=len(o['catalogues'])
    if o['names'] is None:
        o['names']=[os.path.basename(x).replace('.fits','') for x in o['catalogues']]
    if o['radii'] is None:
        o['radii']=[10]*cl
    if o['groups'] is None:
        o['groups']=range(cl)
    if (len(o['frequencies'])!=cl or len(o['radii'])!=cl or
        len(o['names'])!=cl or len(o['groups'])!=cl):
        die('Names, groups, radii and frequencies entries must be the same length as the catalogue list')

    low_robust=-0.25
    low_uvrange=[0.1,25.0]

    # Clear the shared memory
    run('CleanSHM.py',dryrun=o['dryrun'])

    # Get the frequencies -- need to take this from the MSs

    mslist=[s.strip() for s in open(o['mslist']).readlines()]
    freqs=[]
    for ms in mslist:
        t = pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False)
        freqs.append(t[0]['REF_FREQUENCY'])

    # sort to work in frequency order

    freqs,mslist = (list(x) for x in zip(*sorted(zip(freqs, mslist), key=lambda pair: pair[0])))

    for f,m in zip(freqs,mslist):
        print m,f

    # Clean in cube mode

    ddf_image('image_bootstrap',o['mslist'],cleanmode='SSD',ddsols='killms_p1',applysols='P',majorcycles=4,robust=low_robust,uvrange=low_uvrange,beamsize=20,imsize=o['bsimsize'],cellsize=o['bscell'],options=o,colname=o['colname'],automask=True,automask_threshold=15,smooth=True,cubemode=True)

    if os.path.isfile('image_bootstrap.cube.int.restored.pybdsm.srl'):
        warn('Source list exists, skipping source extraction')
    else:
        warn('Running PyBDSM, please wait...')
        img=bdsm.process_image('image_bootstrap.cube.int.restored.fits',thresh_pix=5,rms_map=True,atrous_do=True,atrous_jmax=2,group_by_isl=True,rms_box=(80,20), adaptive_rms_box=True, adaptive_thresh=80, rms_box_bright=(35,7),mean_map='zero',spectralindex_do=True,specind_maxchan=1,debug=True,kappa_clip=3,flagchan_rms=False,flagchan_snr=False,incl_chan=True,spline_rank=1)
        # Write out in ASCII to work round bug in pybdsm
        img.write_catalog(catalog_type='srl',format='ascii',incl_chan='true')
        img.export_image(img_type='rms',img_format='fits')

    from make_fitting_product import make_catalogue
    import fitting_factors
    import find_outliers

    # generate the fitting product
    if os.path.isfile('crossmatch-1.fits'):
        warn('Crossmatch table exists, skipping crossmatch')
    else:
 
        t = pt.table(mslist[0]+ '/FIELD', readonly=True, ack=False)
        direction = t[0]['PHASE_DIR']
        ra, dec = direction[0]

        if (ra<0):
            ra+=2*np.pi
        ra*=180.0/np.pi
        dec*=180.0/np.pi

        cats=zip(o['catalogues'],o['names'],o['groups'],o['radii'])
        make_catalogue('image_bootstrap.cube.int.restored.pybdsm.srl',ra,dec,2.5,cats)
    
    freqlist=open('frequencies.txt','w')
    for n,f in zip(o['names'],o['frequencies']):
        freqlist.write('%f %s_Total_flux %s_E_Total_flux False\n' % (f,n,n))
    for i,f in enumerate(freqs):
        freqlist.write('%f Total_flux_ch%i E_Total_flux_ch%i True\n' % (f,i+1,i+1))
    freqlist.close()

    # Now call the fitting code

    if os.path.isfile('crossmatch-results-1.npy'):
        warn('Results 1 exists, skipping first fit')
    else:
        fitting_factors.run_all(1)

    nreject=-1 # avoid error if we fail somewhere
    if os.path.isfile('crossmatch-2.fits'):
        warn('Second crossmatch exists, skipping outlier rejection')
    else:
        nreject=find_outliers.run_all(1)
    
    if os.path.isfile('crossmatch-results-2.npy'):
        warn('Results 1 exists, skipping second fit')
    else:
        if nreject==0:
            shutil.copyfile('crossmatch-results-1.npy','crossmatch-results-2.npy')
        else:
            fitting_factors.run_all(2)

    # Now apply corrections

    if o['full_mslist'] is None:
        die('Need big mslist to apply corrections')
    if not(o['dryrun']):
        warn('Applying corrections to MS list')
        scale=np.load('crossmatch-results-2.npy')[:,0]
        # InterpolatedUS gives us linear interpolation between points
        # and extrapolation outside it
        spl = InterpolatedUnivariateSpline(freqs, scale, k=1)
        bigmslist=[s.strip() for s in open(o['full_mslist']).readlines()]
        for ms in bigmslist:
            t = pt.table(ms)
            try:
                dummy=t.getcoldesc('SCALED_DATA')
            except RuntimeError:
                dummy=None
            t.close()
            if dummy is not None:
                warn('Table '+ms+' has already been corrected, skipping')
            else:
                t = pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False)
                frq=t[0]['REF_FREQUENCY']
                factor=spl(frq)
                print frq,factor
                t=pt.table(ms,readonly=False)
                desc=t.getcoldesc(o['colname'])
                desc['name']='SCALED_DATA'
                t.addcols(desc)
                d=t.getcol(o['colname'])
                d*=factor
                t.putcol('SCALED_DATA',d)
                t.close()

if __name__=='__main__':
    from options import options
    o=options(sys.argv[1])
    run_bootstrap(o)
