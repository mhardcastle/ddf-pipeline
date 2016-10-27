#!/usr/bin/python

# Code to bootstrap the flux density scale using killMS/DDF
# Currently stand-alone, to incorporate into the pipeline later

# DOES NOT WORK properly until we get beam correction

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

def logfilename(s):
    if o['logging'] is not None:
        return o['logging']+'/'+s 
    else:
        return None

def ddf_image_low(imagename,msname,cleanmask,cleanmode,ddsols,applysols,threshold,majorcycles,dicomodel,robust):
    fname=imagename+'.restored.fits'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF step')
    else:
        runcommand = "DDF.py --ImageName=%s --MSName=%s --NFreqBands=1 --ColName CORRECTED_DATA --NCPU=%i --Mode=Clean --CycleFactor=1.5 --MaxMinorIter=1000000 --MaxMajorIter=%s --MinorCycleMode %s --BeamMode=LOFAR --LOFARBeamMode=A --SaveIms [Residual_i] --Robust %f --Npix=%i --wmax 50000 --Cell %f --UVRangeKm=[0.1,25.0] "%(imagename,msname,o['NCPU_DDF'],majorcycles,cleanmode,robust,o['bsimsize'],o['bscell'])
        if cleanmask != '':
            runcommand += ' --CleanMaskImage=%s'%cleanmask
        if applysols != '':
            runcommand += ' --DDModeGrid=%s --DDModeDeGrid=%s --DDSols=%s'%(applysols,applysols,ddsols)
        if dicomodel != '':
            runcommand += ' --InitDicoModel=%s'%dicomodel
        if threshold != '':
            runcommand += ' --FluxThreshold=%s'%threshold
        run(runcommand,dryrun=o['dryrun'],log=logfilename('DDF-low-'+imagename+'.log'),quiet=o['quiet'])

def make_mask(imagename,thresh):
    fname=imagename+'.mask.fits'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeMask step')
    else:
        runcommand = "MakeMask.py --RestoredIm=%s --Th=%s --Box=50,2"%(imagename,thresh)
        run(runcommand,dryrun=o['dryrun'],log=logfilename('MM-'+imagename+'.log'),quiet=o['quiet'])

def restore(basename,beam):
    fname=basename+'.restoredNew.fits'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping Restore step')
    else:
        runcommand = "Restore.py --BaseImageName=%s --ResidualImage=%s --BeamPix=%f" % (basename, basename+'.residual.fits', beam)
        run(runcommand,dryrun=o['dryrun'],log=logfilename('Restore-'+basename+'.log'),quiet=o['quiet'])

def run_bootstrap(oa):
    global o
    o=oa
    
    if o['mslist'] is None:
        die('MS list must be specified')

    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])

    low_robust=-0.25

    # Clear the shared memory
    #run('CleanSHM.py',dryrun=o['dryrun'])

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

    ddf_image_low('image_low_initial_MSMF',mslist[0],'','MSMF','killms_p1','P',5e-3,3,'',low_robust)
    make_mask('image_low_initial_MSMF.restored.fits',20)

    # now loop over the MSs to make the images
    for i,ms in enumerate(mslist):
        imroot='image_low_%i_GA' % i
        ddf_image_low(imroot,ms,'image_low_initial_MSMF.restored.fits.mask.fits','GA','killms_p1','P','',3,'',low_robust)
        make_mask(imroot+'.restored.fits',15)
        ddf_image_low(imroot+'m',ms,imroot+'.restored.fits.mask.fits','GA','killms_p1','P','',2,imroot+'.DicoModel',low_robust)
        restore(imroot+'m',7.0)

    #make the cube
    if os.path.isfile('cube.fits'):
        warn('Cube file exists, skipping cube assembly')
    else:
        warn('Making the cube')
        make_cube('cube.fits',['image_low_%i_GAm.restoredNew.corr.fits' % i for i in range(len(mslist))],freqs)
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
