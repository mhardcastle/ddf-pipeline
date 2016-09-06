#!/usr/bin/python

# Code to bootstrap the flux density scale using killMS/DDF
# Currently stand-alone, to incorporate into the pipeline later

# DOES NOT WORK properly until we get beam correction

import os,sys
import os.path
from auxcodes import run,find_imagenoise,warn,die
from options import options

def ddf_image_low(imagename,msname,cleanmask,cleanmode,ddsols,applysols,threshold,majorcycles,dicomodel,robust):
    fname=imagename+'.restored.fits'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF step')
    else:
        runcommand = "DDF.py --ImageName=%s --MSName=%s --NFreqBands=1 --ColName CORRECTED_DATA --NCPU=%i --Mode=Clean --CycleFactor=1.5 --MaxMinorIter=1000000 --MaxMajorIter=%s --MinorCycleMode %s --BeamMode=LOFAR --LOFARBeamMode=A --SaveIms [Residual_i] --Robust %f --Npix=%i --wmax 50000 --Cell 5 --UVRangeKm=[0.1,25.0] "%(imagename,msname,o['NCPU_DDF'],majorcycles,cleanmode,robust,6000)
        if cleanmask != '':
            runcommand += ' --CleanMaskImage=%s'%cleanmask
        if applysols != '':
            runcommand += ' --DDModeGrid=%s --DDModeDeGrid=%s --DDSols=%s'%(applysols,applysols,ddsols)
        if dicomodel != '':
            runcommand += ' --InitDicoModel=%s'%dicomodel
        if threshold != '':
            runcommand += ' --FluxThreshold=%s'%threshold
        run(runcommand,dryrun=o['dryrun'],log=o['logging']+'/DDF-low-'+imagename+'.log',quiet=o['quiet'])

def make_mask(imagename,thresh):
    fname=imagename+'.mask.fits'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeMask step')
    else:
        runcommand = "MakeMask.py --RestoredIm=%s --Th=%s --Box=50,2"%(imagename,thresh)
        run(runcommand,dryrun=o['dryrun'],log=o['logging']+'/MM-'+imagename+'.log',quiet=o['quiet'])

if __name__=='__main__':

    # Main loop

    o=options(sys.argv[1])
    if o['mslist'] is None:
        die('MS list must be specified')

    if not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])

    # Clear the shared memory
    #run('CleanSHM.py',dryrun=o['dryrun'])

    # We use the individual ms in mslist.
    mslist=[s.strip() for s in open(o['mslist']).readlines()]
    print 'MS list is', mslist
    low_robust=-0.25

    # First we need to do a MSMF clean to make an initial mask; then
    # we can use this for each of the bands. We use the
    # lowest-frequency dataset, assumed to be first in the mslist for
    # the moment (later we can check this as we will need the MS
    # frequencies anyway)

    ddf_image_low('image_low_initial_MSMF',mslist[0],'','MSMF','killms_p1','P',5e-3,3,'',low_robust)
    make_mask('image_low_initial_MSMF.restored.fits',20)
    
    # now loop over the MSs
    for i,ms in enumerate(mslist):
        imroot='image_low_%i_GA' % i
        ddf_image_low(imroot,ms,'image_low_initial_MSMF.restored.fits.mask.fits','GA','killms_p1','P','',3,'',low_robust)
        make_mask(imroot+'.restored.fits',15)
        ddf_image_low(imroot+'m',ms,imroot+'.restored.fits.mask.fits','GA','killms_p1','P','',2,imroot+'.DicoModel',low_robust)

    # Now restore with matched beam        
    # Restore.py --BaseImageName=test-low-GA --ResidualImage=test-low-GA.residual04.fits --BeamPix=7
