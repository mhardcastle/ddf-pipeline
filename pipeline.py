#!/usr/bin/python

# Routine is to use killms/ddf to selfcalibrate the data using 30 directions.
import os,sys
import os.path
from auxcodes import run,find_imagenoise,warn,die
from options import options

def logfilename(s):
    if o['logging'] is not None:
        return o['logging']+'/'+s 
    else:
        return None

def ddf_image(imagename,mslist,cleanmask,cleanmode,ddsols,applysols,threshold,majorcycles,dicomodel,robust):
    fname=imagename+'.restored.fits'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF step')
    else:
        runcommand = "DDF.py --ImageName=%s --MSName=%s --NFreqBands=2 --ColName %s --NCPU=%i --Mode=Clean --CycleFactor=1.5 --MaxMinorIter=1000000 --MaxMajorIter=%s --MinorCycleMode %s --BeamMode=LOFAR --LOFARBeamMode=A --SaveIms [Residual_i] --Robust %f --Npix=%i --wmax 50000 --Cell 1.5 --NFacets=11 "%(imagename,mslist,colname,o['NCPU_DDF'],majorcycles,cleanmode,robust,o['imsize'])
        if cleanmask != '':
            runcommand += ' --CleanMaskImage=%s'%cleanmask
        if applysols != '':
            runcommand += ' --GlobalNorm=MeanAbs--DDModeGrid=%s --DDModeDeGrid=%s --DDSols=%s'%(applysols,applysols,ddsols)
        if dicomodel != '':
            runcommand += ' --InitDicoModel=%s'%dicomodel
        if threshold != '':
            runcommand += ' --FluxThreshold=%s'%threshold
        run(runcommand,dryrun=o['dryrun'],log=logfilename('DDF-'+imagename+'.log'),quiet=o['quiet'])

def make_mask(imagename,thresh):
    fname=imagename+'.mask.fits'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeMask step')
    else:
        runcommand = "MakeMask.py --RestoredIm=%s --Th=%s --Box=50,2"%(imagename,thresh)
        run(runcommand,dryrun=o['dryrun'],log=logfilename('MM-'+imagename+'.log'),quiet=o['quiet'])

def killms_data(imagename,mslist,outsols,clusterfile):
    with open(mslist,'r') as f:
        mslistname = [l for l in (line.strip() for line in f) if l][-1]  # last non-empty line
    checkname=mslistname+'/killMS.'+outsols+'.sols.npz'
    if o['restart'] and os.path.isfile(checkname):
        warn('Solutions file '+checkname+' already exists, not running killMS step')
    else:
        if imagename != '':
            runcommand = "killMS.py --MSName %s --SolverType KAFCA --PolMode Scalar --BaseImageName %s --dt %i --Weighting Natural --BeamMode LOFAR --LOFARBeamMode=A --NIterKF 6 --CovQ 0.1 --LambdaKF=%f --NCPU %i --OutSolsName %s --NChanSols %i --InCol CORRECTED_DATA"%(mslist,imagename,o['dt'], o['LambdaKF'], o['NCPU_killms'], outsols, o['NChanSols'])
            if clusterfile != '':
                runcommand+=' --NodesFile '+clusterfile
        else:
            # in current code, not used
            runcommand = "killMS.py --MSName %s --SolverType KAFCA --PolMode Scalar --SkyModel %s --dt %i --Weighting Natural --BeamMode LOFAR --LOFARBeamMode=A --NIterKF 6 --CovQ 0.1 --NCPU %i --OutSolsName %s --NChanSols %i --InCol CORRECTED_DATA"%(mslist,skymodel, o['dt'] ,o['NCPU_killms'],outsols, o['NChanSols'])
        run(runcommand,dryrun=o['dryrun'],log=logfilename('KillMS-'+outsols+'.log'),quiet=o['quiet'])

def make_model(maskname,imagename):
    fname=imagename+'.npy'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeModel step')
    else:
        runcommand = "MakeModel.py --MaskName=%s --BaseImageName=%s --NCluster=%i --DoPlot=0"%(maskname,imagename,o['ndir'])
        run(runcommand,dryrun=o['dryrun'],log=logfilename('MakeModel-'+maskname+'.log'),quiet=o['quiet'])


if __name__=='__main__':
    # Main loop

    colname='CORRECTED_DATA'

    o=options(sys.argv[1])
    if o['mslist'] is None:
        die('MS list must be specified')

    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])

    # Clear the shared memory
    run('CleanSHM.py',dryrun=o['dryrun'])

    # Image full bandwidth to create a model
    ddf_image('image_dirin_MSMF',o['mslist'],'','MSMF','','',50E-3,3,'',o['robust'])
    make_mask('image_dirin_MSMF.app.restored.fits',o['ga'])
    #imagenoise = find_imagenoise('image_dirin_MSMF.restored.fits',1E-3)
    ddf_image('image_dirin_GAm',o['mslist'],'image_dirin_MSMF.app.restored.fits.mask.fits','GA','','','',4,'',o['robust'])
    make_mask('image_dirin_GAm.app.restored.fits',o['ga'])

    # Calibrate off the model
    make_model('image_dirin_GAm.restored.fits.mask.fits','image_dirin_GAm')
    killms_data('image_dirin_GAm',o['mslist'],'killms_p1','image_dirin_GAm.npy.ClusterCat.npy')

    # now if bootstrapping has been done then change the column name

    if o['bootstrap']:
        run('bootstrap.py '+sys.argv[1],dryrun=o['dryrun'],log=None)
        colname='SCALED_DATA'

    # Apply phase solutions and image again
    ddf_image('image_phase1',o['mslist'],'image_dirin_GAm.app.restored.fits.mask.fits','GA','killms_p1','P','',3,'',o['robust'])
    make_mask('image_phase1.app.restored.fits',o['phase'])
    ddf_image('image_phase1m',o['mslist'],'image_phase1.restored.fits.mask.fits','GA','killms_p1','P','',2,'image_phase1.DicoModel',o['robust'])
    make_mask('image_phase1m.app.restored.fits',o['phase'])

    # Calibrate off the model
    killms_data('image_phase1m',o['mslist'],'killms_ap1','')

    # Apply phase and amplitude solutions and image again
    ddf_image('image_ampphase1',o['mslist'],'image_phase1m.app.restored.fits.mask.fits','GA','killms_ap1','AP','',3,'',o['robust'])
    make_mask('image_ampphase1.app.restored.fits',o['ampphase'])
    ddf_image('image_ampphase1m',o['mslist'],'image_ampphase1.app.restored.fits.mask.fits','GA','killms_ap1','AP','',2,'image_ampphase1.DicoModel',o['robust'])
    make_mask('image_ampphase1m.app.restored.fits',o['ampphase'])

    # Now move to the full dataset, if it exists

    if o['full_mslist'] is None:
        warn('No full MS list supplied, stopping here')
    else:
        # single AP cal of full dataset and final image. Is this enough?
        killms_data('image_ampphase1m',o['full_mslist'],'killms_f_ap1','')
        ddf_image('image_full_ampphase1',o['full_mslist'],'image_ampphase1m.app.restored.fits.mask.fits','GA','killms_f_ap1','AP','',3,'',o['final_robust'])
        make_mask('image_full_ampphase1.app.restored.fits',o['full'])
        ddf_image('image_full_ampphase1m',o['full_mslist'],'image_full_ampphase1.app.restored.fits.mask.fits','GA','killms_f_ap1','AP','',3,'image_full_ampphase1.DicoModel',o['final_robust'])

        
