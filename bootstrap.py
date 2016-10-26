#!/usr/bin/python

# Code to bootstrap the flux density scale using killMS/DDF
# Should run standalone or as part of the pipeline

import os,sys
import os.path
from auxcodes import run,find_imagenoise,warn,die
import pyrap.tables as pt
import numpy as np
from lofar import bdsm
from make_cube import make_cube
from make_fitting_product import make_catalogue
import fitting_factors
import find_outliers
from scipy.interpolate import InterpolatedUnivariateSpline
from pipeline import ddf_image, make_mask

def logfilename(s):
    if o['logging'] is not None:
        return o['logging']+'/'+s 
    else:
        return None

def run_bootstrap(oa):
    global o
    o=oa
    
    if o['mslist'] is None:
        die('MS list must be specified')

    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])

    low_robust=-0.25
    low_uvrange=[0.1,25.0]

    # Clear the shared memory
    run('CleanSHM.py',dryrun=o['dryrun'])

    # We use the individual ms in mslist.
    mslist=[s.strip() for s in open(o['mslist']).readlines()]

    # Get the frequencies -- need to take this from the MSs

    freqs=[]
    for ms in mslist:
        t = pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False)
        freqs.append(t[0]['REF_FREQUENCY'])

    # sort to work in frequency order

    freqs,mslist = (list(x) for x in zip(*sorted(zip(freqs, mslist), key=lambda pair: pair[0])))

    for f,m in zip(freqs,mslist):
        print m,f

    # First we need to do a MSMF clean to make an initial mask; then
    # we can use this for each of the bands. We use the
    # lowest-frequency dataset.

    ddf_image('image_low_initial_MSMF',mslist[0],cleanmode='MSMF',ddsols='killms_p1',applysols='P',threshold=5e-3,majorcyles=3,robust=low_robust,uvrange=low_uvrange,beamsize=20)
    make_mask('image_low_initial_MSMF.app.restored.fits',20)

    # now loop over the MSs to make the images
    for i,ms in enumerate(mslist):
        imroot='image_low_%i_GA' % i
        ddf_image(imroot,ms,cleanmask='image_low_initial_MSMF.restored.fits.mask.fits',cleanmode='GA',ddsols='killms_p1',applysols='P',majorcyles=3,robust=low_robust,uvrange=low_uvrange,beamsize=20)
        make_mask(imroot+'.app.restored.fits',15)
        ddf_image(imroot+'m',ms,cleanmask=imroot+'.app.restored.fits.mask.fits',previous_image=imroot,reuse_psf=True,use_dicomodel=True,majorcycles=2,cleanmode='GA',ddsols='killms_p1',applysols='P',majorcyles=3,robust=low_robust,uvrange=low_uvrange,beamsize=20,saveimages='H')

    #make the cube
    if os.path.isfile('cube.fits'):
        warn('Cube file exists, skipping cube assembly')
    else:
        warn('Making the cube')
        make_cube('cube.fits',['image_low_%i_GAm.int.restored.fits' % i for i in range(len(mslist))],freqs)
    if os.path.isfile('cube.pybdsm.srl'):
        warn('Source list exists, skipping source extraction')
    else:
        warn('Running PyBDSM, please wait...')
        img=bdsm.process_image('cube.fits',thresh_pix=5,rms_map=True,atrous_do=True,atrous_jmax=2,group_by_isl=True,rms_box=(80,20), adaptive_rms_box=True, adaptive_thresh=80, rms_box_bright=(35,7),mean_map='zero',spectralindex_do=True,specind_maxchan=1,debug=True,kappa_clip=3,flagchan_rms=False,flagchan_snr=False,incl_chan=True,spline_rank=1)
        # Write out in ASCII to work round bug in pybdsm
        img.write_catalog(catalog_type='srl',format='ascii',incl_chan='true')
        img.export_image(img_type='rms',img_format='fits')

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

        # currently VLSS and WENSS are hard-wired in: this should change...

        cats=[['/stri-data/mjh/bootstrap/VLSS.fits','VLSS',40.0],
              ['/stri-data/mjh/bootstrap/wenss.fits','WENSS',10.0]]
        make_catalogue('cube.pybdsm.srl',ra,dec,2.5,cats)
    
    freqlist=open('frequencies.txt','w')
    freqlist.write('74e6 VLSS_flux VLSS_e_flux False\n')
    for i,f in enumerate(freqs):
        freqlist.write('%f Total_flux_ch%i E_Total_flux_ch%i True\n' % (f,i+1,i+1))
    freqlist.write('326e6 WENSS_flux WENSS_e_flux False\n')
    freqlist.close()

    # Now call the fitting code

    if os.path.isfile('crossmatch-results-1.npy'):
        warn('Results 1 exists, skipping first fit')
    else:
        if o['use_mpi']:
            run('mpiexec -np 24 fitting_factors.py 1',dryrun=o['dryrun'],log=None,quiet=o['quiet'])
        else:
            fitting_factors.run_all(1)

    if os.path.isfile('crossmatch-2.fits'):
        warn('Second crossmatch exists, skipping outlier rejection')
    else:
        find_outliers.run_all(1)
    
    if os.path.isfile('crossmatch-results-2.npy'):
        warn('Results 1 exists, skipping second fit')
    else:
        if o['use_mpi']:
            run('mpiexec -np 24 fitting_factors.py 2',dryrun=o['dryrun'],log=None,quiet=o['quiet'])
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
                desc=t.getcoldesc('CORRECTED_DATA')
                desc['name']='SCALED_DATA'
                t.addcols(desc)
                d=t.getcol('CORRECTED_DATA')
                d*=factor
                t.putcol('SCALED_DATA',d)
                t.close()

if __name__=='__main__':
    from options import options
    o=options(sys.argv[1])
    run_bootstrap(o)
