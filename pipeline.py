#!/usr/bin/python

# Routine is to use killms/ddf to self-calibrate the data
import os,sys
import os.path
from auxcodes import run,find_imagenoise,warn,die
from options import options
from shutil import copyfile,rmtree

def logfilename(s):
    if o['logging'] is not None:
        return o['logging']+'/'+s 
    else:
        return None

"""
DDF2 changes to ddf_image:
1) new parameters for DDF2
2) changes to filenames of output files
3) re-use cache second time round
4) Force resolution to something sensible.
"""

def ddf_image(imagename,mslist,cleanmask=None,cleanmode='MSMF',ddsols=None,applysols=None,threshold=None,majorcycles=3,previous_image=None,use_dicomodel=False,robust=0,beamsize=None,reuse_psf=False,reuse_dirty=False,verbose=False,saveimages=None,imsize=None,cellsize=None,uvrange=None,colname='CORRECTED_DATA'):
    # saveimages lists _additional_ images to save
    if saveimages is None:
        saveimages=''
    saveimages+='onNeds'

    if beamsize is None:
        beamsize=o['psf_arcsec']
    if imsize is None:
        imsize=o['imsize']
    if cellsize is None:
        cellsize=o['cellsize']

    fname=imagename+'.app.restored.fits'

    runcommand = "DDF.py --ImageName=%s --MSName=%s --NFreqBands=2 --ColName %s --NCPU=%i --Mode=Clean --CycleFactor=0 --MaxMinorIter=1000000 --MaxMajorIter=%s --MinorCycleMode %s --BeamMode=LOFAR --LOFARBeamMode=A --SaveIms [Residual_i] --Robust %f --Npix=%i --wmax 50000 --Nw 100 --SaveImages %s --Cell %f --NFacets=11 --NEnlargeData 0 --NChanDegridPerMS 1 --RestoringBeam %f"%(imagename,mslist,colname,o['NCPU_DDF'],majorcycles,cleanmode,robust,imsize,saveimages,cellsize,beamsize)
    if cleanmode == 'GA':
        runcommand += ' --GASolvePars [S,Alpha] --BICFactor 0'
    if cleanmask is not None:
        runcommand += ' --CleanMaskImage=%s'%cleanmask
    if applysols is not None:
        runcommand += ' --GlobalNorm=MeanAbs --DDModeGrid=%s --DDModeDeGrid=%s --DDSols=%s'%(applysols,applysols,ddsols)
    if previous_image is not None:
        runcommand += ' --InitDicoModel=%s.DicoModel' % previous_image
    if threshold is not None:
        runcommand += ' --FluxThreshold=%f'%threshold
    if uvrange is not None:
        runcommand += ' --UVRangeKM=[%f,%f]' % (uvrange[0],uvrange[1])
    if reuse_dirty:
        # possible that crashes could destroy the cache, so need to check
        if os.path.exists(mslist+'.ddfcache/Dirty'):
            runcommand += ' --ResetDirty=-1'
    if reuse_psf:
        if os.path.exists(mslist+'.ddfcache/PSF'):
            runcommand += ' --ResetPSF=-1'

    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF step')
        if verbose:
            print 'would have run',runcommand
    else:
         run(runcommand,dryrun=o['dryrun'],log=logfilename('DDF-'+imagename+'.log'),quiet=o['quiet'])

    if previous_image is not None and (reuse_dirty or reuse_psf):
        # this may mean that the Norm and dirty images won't be
        # written. They're the same as for the previous image, so just
        # make links:
        try:
            os.symlink(previous_image+'.dirty.fits',imagename+'.dirty.fits')
            os.symlink(previous_image+'.Norm.fits',imagename+'.Norm.fits')
        except OSError:
            # if the symlink failed the files were there already, so pass
            pass

def make_mask(imagename,thresh,verbose=False):
    fname=imagename+'.mask.fits'
    runcommand = "MakeMask.py --RestoredIm=%s --Th=%s --Box=50,2"%(imagename,thresh)
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeMask step')
        if verbose:
            print 'Would have run',runcommand
    else:
        run(runcommand,dryrun=o['dryrun'],log=logfilename('MM-'+imagename+'.log'),quiet=o['quiet'])

def killms_data(imagename,mslist,outsols,clusterfile,colname='CORRECTED_DATA'):
    # run killms individually on each MS -- allows restart if it failed in the middle
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    for f in filenames:
        checkname=f+'/killMS.'+outsols+'.sols.npz'
        if o['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running killMS step')
        else:
            runcommand = "killMS.py --MSName %s --SolverType KAFCA --PolMode Scalar --BaseImageName %s --dt %i --Weighting Natural --BeamMode LOFAR --LOFARBeamMode=A --NIterKF 6 --CovQ 0.1 --LambdaKF=%f --NCPU %i --OutSolsName %s --NChanSols %i --InCol %s"%(f,imagename,o['dt'], o['LambdaKF'], o['NCPU_killms'], outsols, o['NChanSols'],colname)
            if clusterfile != '':
                runcommand+=' --NodesFile '+clusterfile
            run(runcommand,dryrun=o['dryrun'],log=logfilename('KillMS-'+f+'_'+outsols+'.log'),quiet=o['quiet'])


def make_model(maskname,imagename):
    # returns True if the step was run, False if skipped
    fname=imagename+'.npy'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeModel step')
        return False
    else:
        runcommand = "MakeModel.py --MaskName=%s --BaseImageName=%s --NCluster=%i --DoPlot=0"%(maskname,imagename,o['ndir'])
        run(runcommand,dryrun=o['dryrun'],log=logfilename('MakeModel-'+maskname+'.log'),quiet=o['quiet'])
        return True

def clearcache(mslist):
    try:
        rmtree(mslist+'.ddfcache')
    except OSError:
        pass

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
    if o['clearcache']:
        # Clear the cache, we don't know where it's been
        clearcache(o['mslist'])
        if o['full_mslist'] is not None:
            clearcache(o['full_mslist'])

    # Image full bandwidth to create a model
    ddf_image('image_dirin_MSMF',o['mslist'],cleanmode='MSMF',threshold=50e-3,majorcycles=3,robust=o['robust'])
    make_mask('image_dirin_MSMF.app.restored.fits',o['ga'])
    #imagenoise = find_imagenoise('image_dirin_MSMF.restored.fits',1E-3)
    ddf_image('image_dirin_GAm',o['mslist'],cleanmask='image_dirin_MSMF.app.restored.fits.mask.fits',cleanmode='GA',majorcycles=4,robust=o['robust'],previous_image='image_dirin_MSMF',reuse_psf=True,reuse_dirty=True)
    make_mask('image_dirin_GAm.app.restored.fits',o['ga'])

    # Calibrate off the model
    if make_model('image_dirin_GAm.app.restored.fits.mask.fits','image_dirin_GAm'):
        # if this step runs, clear the cache to remove facet info
        clearcache(o['mslist'])
    
    killms_data('image_dirin_GAm',o['mslist'],'killms_p1','image_dirin_GAm.npy.ClusterCat.npy')

    # now if bootstrapping has been done then change the column name

    if o['bootstrap']:
        run('bootstrap.py '+sys.argv[1],dryrun=o['dryrun'],log=None)
        colname='SCALED_DATA'

    # Apply phase solutions and image again
    ddf_image('image_phase1',o['mslist'],cleanmask='image_dirin_GAm.app.restored.fits.mask.fits',cleanmode='GA',ddsols='killms_p1',applysols='P',majorcycles=3,robust=o['robust'],colname=colname)
    make_mask('image_phase1.app.restored.fits',o['phase'])
    ddf_image('image_phase1m',o['mslist'],cleanmask='image_phase1.app.restored.fits.mask.fits',cleanmode='GA',ddsols='killms_p1',applysols='P',majorcycles=3,previous_image='image_phase1',robust=o['robust'],reuse_psf=True,use_dicomodel=True,colname=colname)
    make_mask('image_phase1m.app.restored.fits',o['phase'])

    # Calibrate off the model
    killms_data('image_phase1m',o['mslist'],'killms_ap1','',colname=colname)

    # Apply phase and amplitude solutions and image again
    ddf_image('image_ampphase1',o['mslist'],cleanmask='image_phase1m.app.restored.fits.mask.fits',cleanmode='GA',ddsols='killms_ap1',applysols='AP',majorcycles=3,robust=o['robust'],colname=colname)
    make_mask('image_ampphase1.app.restored.fits',o['ampphase'])
    ddf_image('image_ampphase1m',o['mslist'],cleanmask='image_ampphase1.app.restored.fits.mask.fits',cleanmode='GA',ddsols='killms_ap1',applysols='AP',majorcycles=2,previous_image='image_ampphase1',use_dicomodel=True,robust=o['robust'],reuse_psf=True,colname=colname)
    make_mask('image_ampphase1m.app.restored.fits',o['ampphase'])

    # Now move to the full dataset, if it exists

    if o['full_mslist'] is None:
        warn('No full MS list supplied, stopping here')
    else:
        # single AP cal of full dataset and final image. Is this enough?
        killms_data('image_ampphase1m',o['full_mslist'],'killms_f_ap1','',colname=colname)
        ddf_image('image_full_ampphase1',o['full_mslist'],cleanmask='image_ampphase1m.app.restored.fits.mask.fits',cleanmode='GA',ddsols='killms_f_ap1',applysols='AP',majorcycles=3,robust=o['final_robust'],colname=colname)
        make_mask('image_full_ampphase1.app.restored.fits',o['full'])
        ddf_image('image_full_ampphase1m',o['full_mslist'],cleanmask='image_full_ampphase1.app.restored.fits.mask.fits',cleanmode='GA',ddsols='killms_f_ap1',applysols='AP',majorcycles=3,previous_image='image_full_ampphase1',use_dicomodel=True,robust=o['final_robust'],beamsize=o['final_psf_arcsec'],reuse_psf=True,saveimages='H',colname=colname)
