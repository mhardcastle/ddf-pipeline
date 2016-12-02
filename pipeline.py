#!/usr/bin/python

# Routine is to use killms/ddf to self-calibrate the data
import os,sys
import os.path
from auxcodes import report,run,find_imagenoise,warn,die
from options import options,print_options
from shutil import copyfile,rmtree
import pyrap.tables as pt
from modify_mask import modify_mask

def logfilename(s,options=None):
    if options is None:
        options=o
    if options['logging'] is not None:
        return options['logging']+'/'+s 
    else:
        return None

def check_imaging_weight(mslist_name):

    report('Checking for IMAGING_WEIGHT in input MSS')
    mslist=[s.strip() for s in open(mslist_name).readlines()]
    for ms in mslist:
        t = pt.table(ms)
        try:
            dummy=t.getcoldesc('IMAGING_WEIGHT')
        except RuntimeError:
            dummy=None
        t.close()
        if dummy is not None:
            warn('Table '+ms+' already has imaging weights')
        else:
            pt.addImagingColumns(ms)

def ddf_image(imagename,mslist,cleanmask=None,cleanmode='MSMF',ddsols=None,applysols=None,threshold=None,majorcycles=3,previous_image=None,use_dicomodel=False,robust=0,beamsize=None,reuse_psf=False,reuse_dirty=False,verbose=False,saveimages=None,imsize=None,cellsize=None,uvrange=None,colname='CORRECTED_DATA',peakfactor=0.1,dicomodel_base=None,options=None,singlefreq=False,do_decorr=None,donorm=True,dirty_from_resid=False,clusterfile=None):
    # saveimages lists _additional_ images to save
    if saveimages is None:
        saveimages=''
    saveimages+='onNeds'
    if options is None:
        options=o # attempt to get global if it exists

    if do_decorr is None:
        do_decorr=options['do_decorr']
    if beamsize is None:
        beamsize=options['psf_arcsec']
    if imsize is None:
        imsize=options['imsize']
    if cellsize is None:
        cellsize=options['cellsize']

    fname=imagename+'.app.restored.fits'

    runcommand = "DDF.py --ImageName=%s --MSName=%s --PeakFactor %f --NFreqBands=%i --ColName %s --NCPU=%i --Mode=Clean --CycleFactor=0 --MaxMinorIter=1000000 --MaxMajorIter=%s --MinorCycleMode %s --BeamMode=LOFAR --LOFARBeamMode=A --SaveIms [Residual_i] --Robust %f --Npix=%i --wmax 50000 --Nw 100 --SaveImages %s --Cell %f --NFacets=11 --NEnlargeData 0 --NChanDegridPerMS 1 --RestoringBeam %f"%(imagename,mslist,peakfactor,1 if singlefreq else 2,colname,options['NCPU_DDF'],majorcycles,cleanmode,robust,imsize,saveimages,cellsize,beamsize)
    if do_decorr:
        runcommand += ' --DecorrMode=FT'
    if cleanmode == 'SSD':
        if singlefreq:
            runcommand += ' --SSDSolvePars [S] --BICFactor 0'
        else:
            runcommand += ' --SSDSolvePars [S,Alpha] --BICFactor 0'
    if clusterfile is not None:
        runcommand += ' --CatNodes=%s' % clusterfile
    if cleanmask is not None:
        runcommand += ' --CleanMaskImage=%s'%cleanmask
    if applysols is not None:
        if donorm:
            runcommand += ' --GlobalNorm=MeanAbs'
        runcommand += ' --DDModeGrid=%s --DDModeDeGrid=%s --DDSols=%s'%(applysols,applysols,ddsols)
    if use_dicomodel:
        if dicomodel_base is not None:
            runcommand += ' --InitDicoModel=%s.DicoModel' % dicomodel_base
        elif previous_image is not None:
            runcommand += ' --InitDicoModel=%s.DicoModel' % previous_image
        else:
            raise RuntimeError('use_dicomodel is set but no dicomodel supplied')
    if threshold is not None:
        runcommand += ' --FluxThreshold=%f'%threshold
    if uvrange is not None:
        runcommand += ' --UVRangeKm=[%f,%f]' % (uvrange[0],uvrange[1])
    if dirty_from_resid:
        # possible that crashes could destroy the cache, so need to check
        if os.path.exists(mslist+'.ddfcache/LastResidual'):
            runcommand += ' --DirtyFromLastResid=1 --ResetDirty=-1'
    if reuse_dirty:
        if os.path.exists(mslist+'.ddfcache/Dirty'):
            runcommand += ' --ResetDirty=-1'
    if reuse_psf:
        if os.path.exists(mslist+'.ddfcache/PSF'):
            runcommand += ' --ResetPSF=-1'

    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF step')
        if verbose:
            print 'would have run',runcommand
    else:
         run(runcommand,dryrun=options['dryrun'],log=logfilename('DDF-'+imagename+'.log',options=options),quiet=options['quiet'])

def make_mask(imagename,thresh,verbose=False,use_tgss=False,options=None):
    if options is None:
        options=o # attempt to get global
    fname=imagename+'.mask.fits'
    runcommand = "MakeMask.py --RestoredIm=%s --Th=%s --Box=50,2"%(imagename,thresh)
    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeMask step')
        if verbose:
            print 'Would have run',runcommand
    else:
        run(runcommand,dryrun=options['dryrun'],log=logfilename('MM-'+imagename+'.log',options=options),quiet=options['quiet'])
        if use_tgss and options['tgss'] is not None:
            report('Merging the mask with TGSS catalogue')
            # TGSS path is provided, this means we want to add the positions of bright TGSS sources to the mask
            modify_mask(fname,fname,options['tgss'],options['tgss_radius'],options['tgss_flux'])

def killms_data(imagename,mslist,outsols,clusterfile=None,colname='CORRECTED_DATA',stagedir=None,dicomodel=None):
    # run killms individually on each MS -- allows restart if it failed in the middle
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    for f in filenames:
        checkname=f+'/killMS.'+outsols+'.sols.npz'
        if o['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running killMS step')
        else:
            dostage=False
            if stagedir is not None and not(os.path.exists(f+'.ddfcache')):
                # here we assume that if ddfcache exists it's because killms has already been run -- so we don't stage if so
                dostage=True
                print 'Staging to',stagedir
                os.system('rsync -a --progress %s %s' % (f,stagedir))
                    
            runcommand = "killMS.py --MSName %s%s --SolverType KAFCA --PolMode Scalar --BaseImageName %s --dt %i --Weighting Natural --BeamMode LOFAR --LOFARBeamMode=A --NIterKF 6 --CovQ 0.1 --LambdaKF=%f --NCPU %i --OutSolsName %s --NChanSols %i --InCol %s"%(stagedir+'/' if dostage else '',f,imagename,o['dt'], o['LambdaKF'], o['NCPU_killms'], outsols, o['NChanSols'],colname)
            if clusterfile is not None:
                runcommand+=' --NodesFile '+clusterfile
            if dicomodel is not None:
                runcommand+=' --DicoModel '+dicomodel
            run(runcommand,dryrun=o['dryrun'],log=logfilename('KillMS-'+f+'_'+outsols+'.log'),quiet=o['quiet'])
            if dostage:
                print 'Staging back'
                os.system('rsync -a --progress %s/%s .' % (stagedir,f))
                os.system('rsync -a --progress %s/%s.ddfcache .' % (stagedir,f))
                os.rmtree(stagedir+'/'+f)

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

def mask_dicomodel(indico,maskname,outdico):
    if o['restart'] and os.path.isfile(outdico):
        warn('File '+outdico+' already exists, skipping MaskDicoModel step')
        return False
    else:
        runcommand = "MaskDicoModel.py --MaskName=%s --InDicoModel=%s --OutDicoModel=%s"%(maskname,indico,outdico) 
        run(runcommand,dryrun=o['dryrun'],log=logfilename('MaskDicoModel-'+maskname+'.log'),quiet=o['quiet'])
        return True

def clearcache(mslist):
    report('Clearing cache for '+mslist)
    try:
        rmtree(mslist+'.ddfcache')
    except OSError:
        pass
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    for f in filenames:
        try:
            rmtree(f+'.ddfcache')
        except OSError:
            pass


if __name__=='__main__':
    # Main loop
    if len(sys.argv)<2:
        warn('pipeline.py must be called with a parameter file.\nSee below for a complete list of possible options with their default values.')
        print_options()
        sys.exit(1)

    o=options(sys.argv[1])
    if o['mslist'] is None:
        die('MS list must be specified')

    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])

    # Set column name for first steps
    colname=o['colname']

    # Clear the shared memory
    run('CleanSHM.py',dryrun=o['dryrun'])    
    if o['clearcache']:
        # Clear the cache, we don't know where it's been
        clearcache(o['mslist'])
        if o['full_mslist'] is not None:
            clearcache(o['full_mslist'])

    # Check imaging weights -- needed before DDF
    check_imaging_weight(o['mslist'])

    # Image full bandwidth to create a model
    ddf_image('image_dirin_MSMF',o['mslist'],cleanmode='MSMF',threshold=50e-3,majorcycles=3,robust=o['robust'],colname=colname)
    make_mask('image_dirin_MSMF.app.restored.fits',o['ga'],use_tgss=True)

    # cluster to get facets
    if make_model('image_dirin_MSMF.app.restored.fits.mask.fits','image_dirin_MSMF'):
        # if this step runs, clear the cache to remove facet info
        clearcache(o['mslist'])

    # Now SSD clean with the new facets
    ddf_image('image_dirin_SSDm',o['mslist'],cleanmask='image_dirin_MSMF.app.restored.fits.mask.fits',cleanmode='SSD',majorcycles=4,robust=o['robust'],previous_image='image_dirin_MSMF',reuse_psf=True,reuse_dirty=True,peakfactor=0.05,colname=colname,clusterfile='image_dirin_MSMF.npy.ClusterCat.npy')
    make_mask('image_dirin_SSDm.app.restored.fits',o['ga'],use_tgss=True)

    # now remove old, bad components from the DicoModel -- these are not in the new mask
    mask_dicomodel('image_dirin_SSDm.DicoModel','image_dirin_SSDm.app.restored.fits.mask.fits','image_dirin_SSDm_masked.DicoModel')

    killms_data('image_dirin_SSDm',o['mslist'],'killms_p1',colname=colname,dicomodel='image_dirin_SSDm_masked.DicoModel')

    # now if bootstrapping has been done then change the column name
    if o['bootstrap']:
#        from bootstrap import run_bootstrap
        report('Running bootstrap')
        run('bootstrap.py '+sys.argv[1],log=None)
#        run_bootstrap(o)
        colname='SCALED_DATA'

    # Apply phase solutions and image again
    ddf_image('image_phase1',o['mslist'],cleanmask='image_dirin_SSDm.app.restored.fits.mask.fits',cleanmode='SSD',ddsols='killms_p1',applysols='P',majorcycles=2,robust=o['robust'],colname=colname,use_dicomodel=True,dicomodel_base='image_dirin_SSDm_masked',peakfactor=0.01)
    make_mask('image_phase1.app.restored.fits',o['phase'],use_tgss=True)
    ddf_image('image_phase1m',o['mslist'],cleanmask='image_phase1.app.restored.fits.mask.fits',cleanmode='SSD',ddsols='killms_p1',applysols='P',majorcycles=3,previous_image='image_phase1',robust=o['robust'],reuse_psf=True,dirty_from_resid=True,use_dicomodel=True,colname=colname,peakfactor=0.01)
    make_mask('image_phase1m.app.restored.fits',o['phase'],use_tgss=True)
    mask_dicomodel('image_phase1m.DicoModel','image_phase1m.app.restored.fits.mask.fits','image_phase1m_masked.DicoModel')
    # Calibrate off the model
    killms_data('image_phase1m',o['mslist'],'killms_ap1',colname=colname,dicomodel='image_phase1m_masked.DicoModel')

    # Apply phase and amplitude solutions and image again
    ddf_image('image_ampphase1',o['mslist'],cleanmask='image_phase1m.app.restored.fits.mask.fits',cleanmode='SSD',ddsols='killms_ap1',applysols='AP',majorcycles=2,robust=o['robust'],colname=colname,use_dicomodel=True,dicomodel_base='image_phase1m_masked',peakfactor=0.01)
    make_mask('image_ampphase1.app.restored.fits',o['ampphase'],use_tgss=True)
    ddf_image('image_ampphase1m',o['mslist'],cleanmask='image_ampphase1.app.restored.fits.mask.fits',cleanmode='SSD',ddsols='killms_ap1',applysols='AP',majorcycles=2,previous_image='image_ampphase1',use_dicomodel=True,robust=o['robust'],reuse_psf=True,dirty_from_resid=True,colname=colname,peakfactor=0.01)
    make_mask('image_ampphase1m.app.restored.fits',o['ampphase'],use_tgss=True)
    mask_dicomodel('image_ampphase1m.DicoModel','image_ampphase1m.app.restored.fits.mask.fits','image_ampphase1m_masked.DicoModel')
    # Now move to the full dataset, if it exists

    if o['full_mslist'] is None:
        warn('No full MS list supplied, stopping here')
    else:
        # Check imaging weights -- needed before DDF
        check_imaging_weight(o['full_mslist'])
        # single AP cal of full dataset and final image. Is this enough?
        killms_data('image_ampphase1m',o['full_mslist'],'killms_f_ap1',colname=colname,clusterfile='image_dirin_SSDm.NodesCat.npy',stagedir=o['stagedir'],dicomodel='image_ampphase1m_masked.DicoModel')
        ddf_image('image_full_ampphase1',o['full_mslist'],cleanmask='image_ampphase1m.app.restored.fits.mask.fits',cleanmode='SSD',ddsols='killms_f_ap1',applysols='AP',majorcycles=2,beamsize=o['final_psf_arcsec'],robust=o['final_robust'],colname=colname,use_dicomodel=True,dicomodel_base='image_ampphase1m_masked')
        make_mask('image_full_ampphase1.app.restored.fits',o['full'],use_tgss=True)
        ddf_image('image_full_ampphase1m',o['full_mslist'],cleanmask='image_full_ampphase1.app.restored.fits.mask.fits',cleanmode='SSD',ddsols='killms_f_ap1',applysols='AP',majorcycles=3,previous_image='image_full_ampphase1',use_dicomodel=True,robust=o['final_robust'],beamsize=o['final_psf_arcsec'],reuse_psf=True,dirty_from_resid=True,saveimages='H',colname=colname,peakfactor=0.001)

        if o['low_psf_arcsec'] is not None:
            # low-res reimage requested
            uvrange=[0.1,2.5*206.0/o['low_psf_arcsec']]
            if o['low_imsize'] is not None:
                low_imsize=o['low_imsize'] # allow over-ride
            else:
                low_imsize=o['imsize']*o['cellsize']/o['low_cell']
            # make an MSMF from one dataset as an initial mask. Use
            # the same name as bootstrap does, so if that's run, we
            # have the mask already (but need to make sure these match!)
            mslist=[s.strip() for s in open(o['mslist']).readlines()]
            ddf_image('image_low_initial_MSMF',mslist[0],cleanmode='MSMF',ddsols='killms_f_ap1',applysols='AP',majorcycles=3,threshold=5e-2,robust=o['low_robust'],uvrange=uvrange,beamsize=o['low_psf_arcsec'],imsize=low_imsize,cellsize=o['low_cell'],singlefreq=True)
            make_mask('image_low_initial_MSMF.app.restored.fits',20)
            make_mask('image_full_ampphase1m.app.restored.fits',o['full'])
            ddf_image('image_full_low',o['full_mslist'],cleanmask='image_low_initial_MSMF.app.restored.fits.mask.fits',cleanmode='SSD',ddsols='killms_f_ap1',applysols='AP',majorcycles=2,robust=o['low_robust'],uvrange=uvrange,beamsize=o['low_psf_arcsec'],imsize=low_imsize,cellsize=o['low_cell'],peakfactor=0.05)
            make_mask('image_full_low.app.restored.fits',o['full'])
            ddf_image('image_full_low_m',o['full_mslist'],cleanmask='image_full_low.app.restored.fits.mask.fits',cleanmode='SSD',ddsols='killms_f_ap1',applysols='AP',majorcycles=3,robust=o['low_robust'],uvrange=uvrange,beamsize=o['low_psf_arcsec'],imsize=low_imsize,cellsize=o['low_cell'],peakfactor=0.001,previous_image='image_full_low',use_dicomodel=True,reuse_psf=True,saveimages='H',donorm=False)
