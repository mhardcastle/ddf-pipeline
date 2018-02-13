#!/usr/bin/env python
"""
ddf-pipeline, a pipeline for LOFAR data reduction
Copyright (C) 2017 Martin Hardcastle (mjh@extragalactic.info) and others

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

import os,sys
import os.path
from auxcodes import report,run,find_imagenoise,warn,die,Catcher,dotdict,separator
from parset import option_list
from options import options,print_options
from shutil import copyfile,rmtree,move
import glob
import pyrap.tables as pt
from modify_mask import modify_mask
from make_extended_mask import make_extended_mask,merge_mask,add_manual_mask
from histmsamp import find_uvmin,sumdico
import numpy as np
from astropy.io import fits
from pipeline_version import version
__version__=version()
import datetime
import threading
from archive_old_solutions import do_archive
from remove_bootstrap import remove_columns
from killMS.Other import MyPickle

def summary(o):
    with open('summary.txt','w') as f:
        ts='{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
        f.write('ddf-pipeline completed at '+ts+'\n')
        f.write('ddf-pipeline version was '+__version__+'\n')
        from DDFacet.DDF import report_version as ddf_version
        f.write('DDF version was '+ddf_version()+'\n')
        from killMS.Other.logo import report_version as killms_version
        f.write('killMS version was '+killms_version()+'\n\n')
        f.write('Options dictionary was as follows:\n')
        for k in o:
            f.write("%-20s : %s\n" % (k,str(o[k])))

def logfilename(s,options=None):
    if options is None:
        options=o
    if options['logging'] is not None:
        return options['logging']+'/'+s 
    else:
        return None

def find_cache_dir(options):
    cache_dir=options['cache_dir']

    # allow cache_dir that only exists on some machines to be specified,
    # fall back to working directory otherwise
    if cache_dir is None:
        cache_dir='.'
    elif not os.path.isdir(cache_dir):
        cache_dir='.'
    return cache_dir

def check_imaging_weight(mslist_name):

    # returns a boolean that says whether it did something
    result=False
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
            result=True
    return result

def ddf_shift(imagename,shiftfile,catcher=None,options=None,verbose=False):
    if catcher: catcher.check()
    if options is None:
        options=o # attempt to get global if it exists

    cache_dir=find_cache_dir(options)

    runcommand='DDF.py '+imagename+'.parset --Output-Name='+imagename+'_shift --Output-Mode=RestoreAndShift --Output-ShiftFacetsFile='+shiftfile+' --Predict-InitDicoModel '+imagename+'.DicoModel --Cache-SmoothBeam=force --Cache-Dir='+cache_dir
    
    fname=imagename+'_shift.app.facetRestored.fits'
    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF-shift step')
        if verbose:
            print 'would have run',runcommand
    else:
         run(runcommand,dryrun=options['dryrun'],log=logfilename('DDF-'+imagename+'_shift.log',options=options),quiet=options['quiet'])


def ddf_image(imagename,mslist,cleanmask=None,cleanmode='HMP',ddsols=None,applysols=None,threshold=None,majorcycles=3,use_dicomodel=False,robust=0,beamsize=None,beamsize_minor=None,beamsize_pa=None,reuse_psf=False,reuse_dirty=False,verbose=False,saveimages=None,imsize=None,cellsize=None,uvrange=None,colname='CORRECTED_DATA',peakfactor=0.1,dicomodel_base=None,options=None,do_decorr=None,normalization=None,dirty_from_resid=False,clusterfile=None,HMPsize=None,automask=True,automask_threshold=10.0,smooth=False,noweights=False,cubemode=False,apply_weights=True,catcher=None,rms_factor=3.0,predict_column=None,conditional_clearcache=False,PredictSettings=None,RMSFactorInitHMP=3.):

    if catcher: catcher.check()

    # saveimages lists _additional_ images to save
    if saveimages is None:
        saveimages=''
    saveimages+='onNeds'
    if options is None:
        options=o # attempt to get global if it exists

    if HMPsize is None:
        HMPsize=options['HMPsize']
    if do_decorr is None:
        do_decorr=options['do_decorr']
    if beamsize is None:
        beamsize=options['psf_arcsec']
    if imsize is None:
        imsize=options['imsize']
    if cellsize is None:
        cellsize=options['cellsize']
        
    cache_dir=find_cache_dir(options)

    if majorcycles>0:
        fname=imagename+'.app.restored.fits'
    else:
        fname=imagename+'.dirty.fits'

    if PredictSettings is not None:
        fname="_has_predicted_OK.%s.info"%imagename

    runcommand = "DDF.py --Output-Name=%s --Data-MS=%s --Deconv-PeakFactor %f --Data-ColName %s --Parallel-NCPU=%i --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=%s --Deconv-Mode %s --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust %f --Image-NPix=%i --CF-wmax 50000 --CF-Nw 100 --Output-Also %s --Image-Cell %f --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=%f --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=%s"%(imagename,mslist,peakfactor,colname,options['NCPU_DDF'],majorcycles,cleanmode,robust,imsize,saveimages,float(cellsize),rms_factor,cache_dir)

    runcommand += " --GAClean-RMSFactorInitHMP %f"%RMSFactorInitHMP
    
    if PredictSettings is None:
        runcommand += " --Output-Mode=Clean"
    else:
        runcommand += " --Output-Mode=%s --Predict-ColName %s"%PredictSettings
        
        
    if beamsize_minor is not None:
        runcommand += ' --Output-RestoringBeam %f,%f,%f'%(beamsize,beamsize_minor,beamsize_pa)
    elif beamsize is not None:
        runcommand += ' --Output-RestoringBeam %f'%(beamsize)
    
    if apply_weights:
        runcommand+=' --Weight-ColName="IMAGING_WEIGHT"'
    else:
        runcommand+=' --Weight-ColName="None"'

    if cubemode:
        # number of channels equals number of distinct freqs in data
        freqs=[]
        mss=[l.rstrip() for l in open(mslist).readlines()]
        for ms in mss:
            t = pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False)
            freq=t[0]['REF_FREQUENCY']
            if freq not in freqs:
                freqs.append(freq)
        channels=len(freqs)
        runcommand+=' --Output-Cubes I --Freq-NBand=%i' % channels
    else:
        runcommand+=' --Freq-NBand=2'
    
    if do_decorr:
        runcommand += ' --RIME-DecorrMode=FT'

    if cleanmode == 'SSD':
        runcommand += ' --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0'
    if clusterfile is not None:
        runcommand += ' --Facets-CatNodes=%s' % clusterfile
    if automask:
        runcommand += ' --Mask-Auto=1 --Mask-SigTh=%.2f' % automask_threshold
    if cleanmask is not None:
        runcommand += ' --Mask-External=%s'%cleanmask
    if applysols is not None:
        if normalization is not None:
            if normalization[:3]=='Abs':
                normalization='Mean'+normalization # backward compat. hack
            runcommand += ' --DDESolutions-GlobalNorm='+normalization
        runcommand += ' --DDESolutions-DDModeGrid=%s --DDESolutions-DDModeDeGrid=%s --DDESolutions-DDSols=%s'%(applysols,applysols,ddsols)
    if use_dicomodel:
        if dicomodel_base is not None:
            runcommand += ' --Predict-InitDicoModel=%s.DicoModel' % dicomodel_base
        else:
            raise RuntimeError('use_dicomodel is set but no dicomodel supplied')
    if threshold is not None:
        runcommand += ' --Deconv-FluxThreshold=%f'%threshold
    if uvrange is not None:
        runcommand += ' --Selection-UVRangeKm=[%f,%f]' % (uvrange[0],uvrange[1])
    if dirty_from_resid and reuse_dirty:
        raise RuntimeError('Cannot combine reuse_dirty and dirty_from_resid')
    if dirty_from_resid:
        # possible that crashes could destroy the cache, so need to check
        if os.path.exists(cache_dir+'/'+mslist+'.ddfcache/LastResidual'):
            runcommand += ' --Cache-Dirty forceresidual'
    if reuse_dirty:
        if os.path.exists(cache_dir+'/'+mslist+'.ddfcache/Dirty'):
            runcommand += ' --Cache-Dirty forcedirty'
    if reuse_psf:
        if os.path.exists(cache_dir+'/'+mslist+'.ddfcache/PSF'):
            runcommand += ' --Cache-PSF force'

    if HMPsize is not None:
        runcommand += ' --GAClean-MinSizeInit=%i' % HMPsize

    if options['nobar']:
        runcommand += ' --Log-Boring=1'

    if smooth:
        runcommand += ' --Beam-Smooth=1'

    if predict_column is not None:
        runcommand += ' --Predict-ColName=%s' % predict_column
        
    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF step')
        if verbose:
            print 'would have run',runcommand
    else:
        if conditional_clearcache:
            clearcache(mslist,options)
        run(runcommand,dryrun=options['dryrun'],log=logfilename('DDF-'+imagename+'.log',options=options),quiet=options['quiet'])

        # Ugly way to see if predict has been already done
        if PredictSettings is not None:
            fname=os.system("touch %s"%fname)
        
        
def make_external_mask(fname,templatename,use_tgss=True,options=None,extended_use=None,clobber=False,cellsize='cellsize'):
    # cellsize specifies which option value to get this from
    if options is None:
        options=o # attempt to get global

    if options['dryrun']: return

    if options['restart'] and os.path.isfile(fname) and not clobber:
        warn('External mask already exists, not creating it')
    else:
        report('Make blank external mask')
        hdus=fits.open(templatename)
        hdus[0].data=np.zeros_like(hdus[0].data,dtype=np.int32)
        hdus.writeto(fname,clobber=True)
        hdus.close()
        if use_tgss and options['tgss'] is not None:
            report('Merging the mask with TGSS catalogue')
            # TGSS path is provided, this means we want to add the positions of bright TGSS sources to the mask
            modify_mask(fname,fname,options['tgss'],options['tgss_radius'],options['tgss_flux'],do_extended=options['tgss_extended'],cellsize=options[cellsize],pointsize=options['tgss_pointlike'])

        if options['region'] is not None:
            report('Merging with mask with user-specified region')
            add_manual_mask(fname,options['region'],fname)

        if options['extended_size'] is not None and extended_use is not None:
            report('Merging with automatic extended mask')
            merge_mask(fname,extended_use,fname)

def make_mask(imagename,thresh,verbose=False,options=None,external_mask=None,catcher=None):
    if catcher: catcher.check()

    # mask_use specifies a mask file to use
    if options is None:
        options=o # attempt to get global

    if options['dryrun']: return
    fname=imagename+'.mask.fits'
    runcommand = "MakeMask.py --RestoredIm=%s --Th=%s --Box=50,2"%(imagename,thresh)
    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeMask step')
        if verbose:
            print 'Would have run',runcommand
    else:
        run(runcommand,dryrun=options['dryrun'],log=logfilename('MM-'+imagename+'.log',options=options),quiet=options['quiet'])
        if external_mask is not None:
            merge_mask(fname,external_mask,fname)
    return fname
            
def killms_data(imagename,mslist,outsols,clusterfile=None,colname='CORRECTED_DATA',niterkf=6,dicomodel=None,uvrange=None,wtuv=None,robust=None,catcher=None,dt=None,options=None,SolverType="KAFCA",PolMode="Scalar",MergeSmooth=False,NChanSols=1,DISettings=None,EvolutionSolFile=None,CovQ=0.1):

    if options is None:
        options=o # attempt to get global if it exists

    cache_dir=find_cache_dir(options)

    # run killms individually on each MS -- allows restart if it failed in the middle
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    for f in filenames:
        if catcher: catcher.check()
        checkname=f+'/killMS.'+outsols+'.sols.npz'
        if o['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running killMS step')
        else:
            
            runcommand = "kMS.py --MSName %s --SolverType %s --PolMode %s --BaseImageName %s --dt %f --NIterKF %i --CovQ %f --LambdaKF=%f --NCPU %i --OutSolsName %s --NChanSols %i --PowerSmooth=%f --InCol %s"%(f,SolverType,PolMode,imagename,dt,niterkf, CovQ, o['LambdaKF'], o['NCPU_killms'], outsols, NChanSols,o['PowerSmooth'],colname)
            if robust is None:
                runcommand+=' --Weighting Natural'
            else:
                runcommand+=' --Weighting Briggs --Robust=%f' % robust
            if uvrange is not None:
                if wtuv is not None:
                    runcommand+=' --WTUV=%f --WeightUVMinMax=%f,%f' % (wtuv, uvrange[0], uvrange[1])
                else:
                    runcommand+=' --UVMinMax=%f,%f' % (uvrange[0], uvrange[1])
            if o['nobar']:
                runcommand+=' --DoBar=0'

            if DISettings is None:
                runcommand+='  --BeamMode LOFAR --LOFARBeamMode=A --DDFCacheDir=%s'%cache_dir
                if clusterfile is not None:
                    runcommand+=' --NodesFile '+clusterfile
                if dicomodel is not None:
                    runcommand+=' --DicoModel '+dicomodel
                if EvolutionSolFile is not None:
                    runcommand+=' --EvolutionSolFile '+EvolutionSolFile
                    
            else:
                runcommand+=" --SolverType %s --PolMode %s --SkyModelCol %s --OutCol %s --ApplyCal 0"%DISettings
                
                
            rootfilename=outsols.split('/')[-1]
            f=f.replace("/","_")
            run(runcommand,dryrun=o['dryrun'],log=logfilename('KillMS-'+f+'_'+rootfilename+'.log'),quiet=o['quiet'])
    if MergeSmooth:
        return smooth_solutions(mslist,outsols,catcher=None,dryrun=o['dryrun'])
    else:
        return outsols

    
def make_model(maskname,imagename,catcher=None):
    # returns True if the step was run, False if skipped
    if catcher: catcher.check()

    fname=imagename+'.npy'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeModel step')
        return False
    else:
        runcommand = "MakeModel.py --MaskName=%s --BaseImageName=%s --NCluster=%i --DoPlot=0"%(maskname,imagename,o['ndir'])
        run(runcommand,dryrun=o['dryrun'],log=logfilename('MakeModel-'+maskname+'.log'),quiet=o['quiet'])
        return True

def mask_dicomodel(indico,maskname,outdico,catcher=None):
    if catcher: catcher.check()

    if o['restart'] and os.path.isfile(outdico):
        warn('File '+outdico+' already exists, skipping MaskDicoModel step')
    else:
        runcommand = "MaskDicoModel.py --MaskName=%s --InDicoModel=%s --OutDicoModel=%s"%(maskname,indico,outdico) 
        run(runcommand,dryrun=o['dryrun'],log=logfilename('MaskDicoModel-'+maskname+'.log'),quiet=o['quiet'])
    return outdico.split(".")[0]

def rmtglob(path):
    g=glob.glob(path)
    for f in g:
        print 'Removing',f
        rmtree(f)

def _basename(path):
    return os.path.basename(path.rstrip(os.path.sep))

def mvglob(path,dest):
    g=glob.glob(path)
    for f in g:
        print 'Moving',f,'to',dest
        # work round shutil non-overwriting behaviour
        real_dst = os.path.join(dest, _basename(f))
        print 'Target is',real_dst
        if os.path.exists(real_dst):
            os.remove(real_dst)
        move(f,dest)

def clearcache(mslist,options):
    cachedir=find_cache_dir(options)

    report('Clearing cache for '+mslist)
    if os.path.isfile(mslist):
        filenames=[l.strip() for l in open(mslist,'r').readlines()]
    else:
        filenames=[]

    try:
        rmtglob(cachedir+'/'+mslist+'*.ddfcache')
        rmtglob(mslist+'*.ddfcache')
    except OSError:
        pass
    for f in filenames:
        try:
            rmtglob(cachedir+'/'+f+'*.ddfcache')
        except OSError:
            pass

def optimize_uvmin(rootname,mslist,colname,uvmin_limit=None):
    uvminfile=rootname+'_uvmin.txt'
    report('Optimizing uvmin for self-cal')
    if os.path.isfile(uvminfile):
        result=float(open(uvminfile).readlines()[0].rstrip())
    else:
        level=sumdico(rootname)
        result=find_uvmin(mslist,level,colname=colname)*1.1
        print 'Will use shortest baseline of',result,'km'
        with open(uvminfile,'w') as f:
            f.write('%f\n' % result)
    if uvmin_limit is not None and result<uvmin_limit:
        result=uvmin_limit
    return result

def smooth_solutions(mslist,ddsols,catcher=None,dryrun=False):
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    with open('solslist.txt','w') as f:
        for fname in filenames:
            solname =fname+'/killMS.'+ddsols+'.sols.npz'
            f.write('%s\n'%(solname))

    checkname='%s_merged.npz'%ddsols
    if o['restart'] and os.path.isfile(checkname):
        warn('Solutions file '+checkname+' already exists, not running MergeSols step')
    else:
        run('MergeSols.py --SolsFilesIn=solslist.txt --SolFileOut=%s_merged.npz'%ddsols,dryrun=dryrun)
    checkname='%s_smoothed.npz'%ddsols
    if o['restart'] and os.path.isfile(checkname):
        warn('Solutions file '+checkname+' already exists, not running SmoothSols step')
    else:
        run('SmoothSols.py --SolsFileIn=%s_merged.npz --SolsFileOut=%s_smoothed.npz --InterpMode=TEC,PolyAmp'%(ddsols,ddsols),dryrun=dryrun)

    return '%s_smoothed.npz'%ddsols

def full_clearcache(o):
    clearcache(o['mslist'],o)
    clearcache('temp_mslist.txt',o)
    if o['full_mslist'] is not None:
        clearcache(o['full_mslist'],o)


def subtract_data(mslist,col1,col2):
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    for f in filenames:
        print 'Subtracting',f
        t = pt.table(f,readonly=False)
        desc=t.getcoldesc(col1)
        desc['name']='SUBTRACTED_DATA'
        t.addcols(desc)
        d1=t.getcol(col1)
        d2=t.getcol(col2)
        t.putcol('SUBTRACTED_DATA',d1-d2)
        t.close()

def main(o=None):
    if o is None:
        o=MyPickle.Load("ddf-pipeline.last")

    if "DDF_PIPELINE_CATALOGS" not in os.environ.keys():
        warn("You need to define the environment variable DDF_PIPELINE_CATALOGS where your catalogs are located")
        sys.exit(2)

    o["tgss"]=o["tgss"].replace("$$",os.environ["DDF_PIPELINE_CATALOGS"])
    o["catalogues"]=[l.replace("$$",os.environ["DDF_PIPELINE_CATALOGS"]) for l in o["catalogues"]]
    lCat=o["catalogues"]+[o["tgss"]]
    for fCat in lCat:
        if not os.path.isfile(fCat):
            warn("Catalog %s does not exist"%fCat)
            sys.exit(2)
          
    
    
    if o['catch_signal']:
        catcher=Catcher()
    else:
        catcher=None

    uvrange=[o['image_uvmin'],o['uvmax']]
    killms_uvrange=[0,1000]
    if o['solutions_uvmin'] is not None:
        killms_uvrange[0]=o['solutions_uvmin']
    if o['mslist'] is None:
        die('MS list must be specified')

    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])

    # Set column name for first steps
    colname=o['colname']

    # Clear the shared memory
    run('CleanSHM.py',dryrun=o['dryrun'])    

    # Check imaging weights -- needed before DDF
    new=check_imaging_weight(o['mslist'])

    if o['clearcache'] or new or o['redofrom']:
        # Clear the cache, we don't know where it's been. If this is a
        # completely new dataset it is always safe (and required) to
        # clear the cache -- solves problems where the cache is not
        # stored per dataset. If we are redoing, cache needs to be removed
        full_clearcache(o)

    if o['redofrom']:

        if not os.path.isdir(o['archive_dir']):
            os.mkdir(o['archive_dir'])

        # If redoing from a specified point, move all relevant files somewhere
        # else, and archive relevant killms solutions

        # we list the stages and the templates of the files that have to
        # be removed to restart from after each one
        stages = [('start', ('image_dirin*','external_mask.fits'), None),
                  ('dirin', ('*bootstrap*', 'image_phase1*', '*crossmatch*', 'external_mask_ext.fits'), 'p1'),
                  ('phase', 'image_ampphase1*', 'ap1'),
                  ('ampphase', ('image_full_low*','full-mask*.fits','external_mask_ext-deep.fits'), 'f_ap1'),
                  ('fulllow', 'image_full_ampphase1*', None),
                  ('full', 'image_full_ampphase2*', 'f_ap2'),
                  ('full2', ('panstarrs-*', 'astromap.fits', 'facet-offset.txt', 'summary.txt'), None)]

        after=False
        bootstrap_removed=False
        alist=[]
        for i,stage in enumerate(stages):
            sname,files,sols=stage
            if sname==o['redofrom']:
                after=True
            if after:
                if not isinstance(files,(list,tuple)):
                    files=(files,)
                for f in files:
                    mvglob(f,o['archive_dir'])
                if sols:
                    alist.append(sols)
                if i<2 and not(bootstrap_removed):
                    warn('Removing bootstrap')
                    if o['full_mslist'] is not None:
                        remove_columns(o['full_mslist'])
                    else:
                        remove_columns(o['mslist'])
                    bootstrap_removed=True
    
        if not after:
            die('Redofrom option not supported')
        else:
            do_archive(o,alist)

    # ##########################################################
    # Initial dirty image to allow an external (TGSS) mask to be made
    separator("Initial dirty")
    ddf_image('image_dirin_SSD_init',o['mslist'],cleanmask=None,cleanmode='SSD',majorcycles=0,robust=o['image_robust'],
              reuse_psf=False,reuse_dirty=False,peakfactor=0.05,colname=colname,clusterfile=None,
              apply_weights=o['apply_weights'][0],uvrange=uvrange,catcher=catcher)

    separator("External mask")
    external_mask='external_mask.fits'
    make_external_mask(external_mask,'image_dirin_SSD_init.dirty.fits',use_tgss=True,clobber=False)

    # Deep SSD clean with this external mask and automasking
    separator("DI Deconv")
    ddf_image('image_dirin_SSD',o['mslist'],cleanmask=external_mask,cleanmode='SSD',
              majorcycles=3,robust=o['image_robust'],reuse_psf=True,reuse_dirty=True,
              peakfactor=0.001,rms_factor=0,
              colname=colname,clusterfile=None,automask=True,
              automask_threshold=o['thresholds'][0],apply_weights=o['apply_weights'][0],
              uvrange=uvrange,catcher=catcher)

    # make a mask from the final image
    separator("MakeMask")
    CurrentMaskName=make_mask('image_dirin_SSD.app.restored.fits',o['thresholds'][0],external_mask=external_mask,catcher=catcher)
    CurrentBaseDicoModelName=mask_dicomodel('image_dirin_SSD.DicoModel','image_dirin_SSD.app.restored.fits.mask.fits','image_dirin_SSD_masked.DicoModel',catcher=catcher)

    # ##########################################################
    # cluster to get facets
    separator("Cluster")
    if not os.path.exists('image_dirin_SSD.Norm.fits') and not os.path.islink("image_dirin_SSD.Norm.fits"):
        os.symlink('image_dirin_SSD_init.Norm.fits','image_dirin_SSD.Norm.fits')
    if not os.path.exists('image_dirin_SSD.dirty.fits') and not os.path.islink("image_dirin_SSD.dirty.fits"):
        os.symlink('image_dirin_SSD_init.dirty.fits','image_dirin_SSD.dirty.fits')
    if make_model('image_dirin_SSD.app.restored.fits.mask.fits','image_dirin_SSD',catcher=catcher):
        # if this step runs, clear the cache to remove facet info
        clearcache(o['mslist'],o)

    # ##########################################################
    if o['auto_uvmin']:
        killms_uvrange[0]=optimize_uvmin('image_dirin_SSD',o['mslist'],colname,o['solutions_uvmin'])

    if o['exitafter'] == 'dirin':
        warn('User specified exit after image_dirin.')
        sys.exit(2)
        
    separator("DD calibration")
    CurrentDDkMSSolName=killms_data('image_dirin_SSD',o['mslist'],'DDS0',colname=colname,
                                    dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                    clusterfile='image_dirin_SSD.npy.ClusterCat.npy',
                                    CovQ=0.02,
                                    niterkf=o['NIterKF'][0],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],dt=o['dt'],
                                    catcher=catcher,NChanSols=o['NChanSols'],
                                    MergeSmooth=True)


    separator("PhaseOnly deconv")
    ddf_image('image_phase1',o['mslist'],
              cleanmask=CurrentMaskName,
              cleanmode='SSD',
              ddsols=CurrentDDkMSSolName,applysols='P',majorcycles=2,robust=o['image_robust'],
              colname=colname,peakfactor=0.001,automask=True,
              automask_threshold=o['thresholds'][1],
              normalization=o['normalize'][0],apply_weights=o['apply_weights'][1],uvrange=uvrange,
              use_dicomodel=True,dicomodel_base=CurrentBaseDicoModelName,catcher=catcher)

    CurrentBaseDicoModelName='image_phase1'
    
    # separator("Mask for deeper deconv")
    # CurrentMaskName=make_mask('image_phase1.app.restored.fits',o['thresholds'][1],external_mask=external_mask,catcher=catcher)
    # CurrentBaseDicoModelName=mask_dicomodel('image_phase1.DicoModel',CurrentMaskName,'image_phase1_masked.DicoModel',catcher=catcher)

    # separator("Phase deeper deconv")
    # ddf_image('image_phase1m',o['mslist'],
    #           cleanmask=CurrentMaskName,
    #           reuse_psf=True,#reuse_dirty=True,dirty_from_resid=True,
    #           cleanmode='SSD',
    #           ddsols=CurrentDDkMSSolName,applysols='P',majorcycles=2,robust=o['image_robust'],
    #           colname=colname,peakfactor=0.001,automask=True,
    #           automask_threshold=o['thresholds'][1],
    #           normalization=o['normalize'][0],apply_weights=o['apply_weights'][1],uvrange=uvrange,
    #           use_dicomodel=True,dicomodel_base=CurrentBaseDicoModelName,catcher=catcher)
    
    separator("DD calibration")
    CurrentDDkMSSolName=killms_data('image_phase1',o['mslist'],'DDS1',colname=colname,
                                    dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                    CovQ=0.02,
                                    clusterfile='image_dirin_SSD.npy.ClusterCat.npy',
                                    niterkf=o['NIterKF'][0],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],dt=o['dt'],
                                    catcher=catcher,NChanSols=o['NChanSols'],
                                    EvolutionSolFile=CurrentDDkMSSolName,
                                    MergeSmooth=True)


    separator("AmpPhase deconv")
    ddf_image('image_ampphase1',o['mslist'],
              cleanmask=CurrentMaskName,cleanmode='SSD',
              ddsols=CurrentDDkMSSolName,applysols='AP',majorcycles=1,robust=o['image_robust'],
              colname=colname,peakfactor=0.001,automask=True,
              automask_threshold=o['thresholds'][1],
              normalization=o['normalize'][0],apply_weights=o['apply_weights'][1],uvrange=uvrange,
              use_dicomodel=True,dicomodel_base='image_phase1',catcher=catcher)


    separator("Update Mask")
    CurrentMaskName=make_mask('image_ampphase1.app.restored.fits',7,external_mask=external_mask,catcher=catcher)
    CurrentBaseDicoModelName=mask_dicomodel('image_ampphase1.DicoModel',CurrentMaskName,'image_ampphase1m_masked.DicoModel',catcher=catcher)
    # CurrentMaskName,CurrentBaseDicoModelName='image_ampphase1m.app.restored.fits.mask.fits','image_ampphase1m_masked.DicoModel'

    # # Compute the DD predict
    # separator("Compute DD predict")
    # ddf_image('DDS0_Predict',o['mslist'],
    #           #cleanmask=external_mask,
    #           cleanmode='SSD',
    #           #majorcycles=1,robust=o['image_robust'],peakfactor=0.01,
    #           colname=colname,uvrange=uvrange,
    #           #automask=True,automask_threshold=o['thresholds'][1],normalization=o['normalize'][0],
    #           #apply_weights=o['apply_weights'][1],
    #           use_dicomodel=True,applysols='AP',
    #           dicomodel_base='image_dirin_SSD_masked',catcher=catcher,
    #           ddsols=CurrentDDkMSSolName, PredictSettings=("Predict","DD_PREDICT"))
    # separator("Calibrate DI")
    # killms_data('image_dirin_SSD',o['mslist'],'DIS0',colname=colname,
    #             dicomodel='image_dirin_SSD_masked.DicoModel',clusterfile='image_dirin_SSD.npy.ClusterCat.npy',
    #             niterkf=o['NIterKF'][0],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
    #             catcher=catcher,
    #             dt=o['dt_di'],
    #             NChanSols=o['NChanSols_di'],
    #             DISettings=("CohJones","IFull","DD_PREDICT","DATA_DI_CORRECTED"))
    # colname_dical="DATA_DI_CORRECTED"

    # ##########################################################
    # run bootstrap, and change the column name if it runs
    if o['bootstrap']:
        separator("Bootstrap")
        report('Running bootstrap')
        run('bootstrap.py '+' '.join(sys.argv[1:]),log=None,dryrun=o["dryrun"])
        colname='SCALED_DATA'

    # ##########################################################
    # make the extended mask if required and possible
    if os.path.isfile('image_bootstrap.app.mean.fits') and o['extended_size'] is not None:
        separator("MakeMask")
        if o['restart'] and os.path.isfile('bootstrap-mask-high.fits'):
            warn('Extended source mask already exists, using existing version')
        else:
            report('Making the extended source mask')
            mask_base_image='image_bootstrap.app.mean.fits'
            make_extended_mask(mask_base_image,'image_dirin_SSD.app.restored.fits',rmsthresh=o['extended_rms'],sizethresh=o['extended_size'],rootname='bootstrap',rmsfacet=o['rmsfacet'])
        external_mask='external_mask_ext.fits'
        make_external_mask(external_mask,'image_dirin_SSD_init.dirty.fits',use_tgss=True,clobber=False,extended_use='bootstrap-mask-high.fits')


    if o['exitafter'] == 'phase':
        warn('User specified exit after image_phase.')
        sys.exit(2)

    # #########################################################################
    # ###############                  BIG MSLIST               ###############
    # #########################################################################

    # ##########################################################
    # Calibrate off the model
    if o['auto_uvmin']:
        killms_uvrange[0]=optimize_uvmin('image_phase1',o['mslist'],colname,o['solutions_uvmin'])

    separator("DD Calibration (full mslist)")
    CurrentDDkMSSolName=killms_data('image_ampphase1',o['full_mslist'],'DDS2_full',
                                    colname=colname,
                                    CovQ=0.02,
                                    dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                    niterkf=o['NIterKF'][1],uvrange=killms_uvrange,wtuv=o['wtuv'],
                                    robust=o['solutions_robust'],dt=o['dt'],
                                    catcher=catcher,
                                    MergeSmooth=True,NChanSols=o['NChanSols'])
    # Compute the DD predict
    separator("Compute DD Predict (full mslist)")
    ddf_image('PredictDSS2',o['full_mslist'],cleanmode='SSD',
              applysols='AP',majorcycles=1,robust=o['image_robust'],colname=colname,peakfactor=0.01,
              automask=True,automask_threshold=o['thresholds'][1],normalization=o['normalize'][0],
              apply_weights=o['apply_weights'][1],uvrange=uvrange,use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              catcher=catcher,
              ddsols=CurrentDDkMSSolName, PredictSettings=("Predict","DD_PREDICT"))

    separator("Compute DI calibration (full mslist)")
    killms_data('PredictDSS2',o['full_mslist'],'DIS0',colname=colname,
                dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                clusterfile='image_dirin_SSD.npy.ClusterCat.npy',
                niterkf=o['NIterKF'][0],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                catcher=catcher,
                dt=o['dt_di'],
                NChanSols=o['NChanSols_di'],
                DISettings=("CohJones","IFull","DD_PREDICT","DATA_DI_CORRECTED"))
    colname="DATA_DI_CORRECTED"

    # ###############################################
    # Apply phase and amplitude solutions and image again
    separator("Deconvolution AP (full mslist)")
    ddf_kw={}
    if o['msss_mode']:
        ddf_kw['cubemode']=True
        ddf_kw['smooth']=True

    ddf_image('image_full_ampphase_di',o['full_mslist'],
              cleanmask=CurrentMaskName,
              cleanmode='SSD',ddsols=CurrentDDkMSSolName,
              applysols='AP',majorcycles=0,robust=o['image_robust'],
              colname=colname,use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][2],
              normalization=o['normalize'][1],uvrange=uvrange,
              apply_weights=o['apply_weights'][2],catcher=catcher,**ddf_kw)

    separator("MakeMask")
    CurrentMaskName=make_mask('image_full_ampphase_di.app.restored.fits',10,external_mask=external_mask,catcher=catcher)

    separator("Finish Deconvolution AP (full mslist)")
    ddf_image('image_full_ampphase_di_m',o['full_mslist'],
              cleanmask=CurrentMaskName,
              reuse_psf=True,dirty_from_resid=True,
              cleanmode='SSD',ddsols=CurrentDDkMSSolName,
              applysols='AP',majorcycles=1,robust=o['image_robust'],
              colname=colname,use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][2],
              normalization=o['normalize'][1],uvrange=uvrange,
              apply_weights=o['apply_weights'][2],catcher=catcher,
              RMSFactorInitHMP=1.,**ddf_kw)

    if o['exitafter'] == 'ampphase':
        warn('User specified exit after image_ampphase.')
        sys.exit(2)

    # Now move to the full dataset, if it exists
    # Check imaging weights -- needed before DDF
    check_imaging_weight(o['full_mslist'])

    if o['auto_uvmin']:
        killms_uvrange[0]=optimize_uvmin('image_full_ampphase1',o['mslist'],colname,o['solutions_uvmin'])

    separator("MakeMask")
    CurrentMaskName=make_mask('image_full_ampphase_di_m.app.restored.fits',7,external_mask=external_mask,catcher=catcher)
    CurrentBaseDicoModelName=mask_dicomodel('image_full_ampphase_di_m.DicoModel',CurrentMaskName,'image_full_ampphase_di_m_masked.DicoModel',catcher=catcher)

    separator("DD Calibration (full mslist)")
    CurrentDDkMSSolName=killms_data('image_full_ampphase_di_m',
                                    o['full_mslist'],'DDS3_full',
                                    colname=colname,
                                    clusterfile='image_dirin_SSD.npy.ClusterCat.npy',
                                    dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                    niterkf=6,
                                    CovQ=0.1,
                                    uvrange=killms_uvrange,
                                    wtuv=o['wtuv'],
                                    robust=o['solutions_robust'],
                                    dt=o['dt'],catcher=catcher)#,EvolutionSolFile=CurrentDDkMSSolName)

    # here we do only image the residuals, and restore so use majorcycles=0
    # (psf is not used so we set reuse_psf=True, so that DDFacet does not recompute it)
    ddf_image('image_full_ampphase_di_m.NS',o['full_mslist'],
              cleanmask=CurrentMaskName,
              reuse_psf=True,
              cleanmode='SSD',ddsols=CurrentDDkMSSolName,
              applysols='AP',majorcycles=0,robust=o['final_robust'],
              colname=colname,use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              beamsize=o['final_psf_arcsec'],
              peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][2],
              normalization=o['normalize'][1],uvrange=uvrange,
              apply_weights=o['apply_weights'][2],catcher=catcher,RMSFactorInitHMP=1.,**ddf_kw)

    return
    us_ddsols=ddsols=CurrentDDkMSSolName
    
    # Do the low-res image first so we can use a mask from it on
    # the high-res image
    
    if o['do_lowres_imaging']:
        # low-res reimage requested
        low_uvrange=[o['image_uvmin'],2.5*206.0/o['low_psf_arcsec']]
        if o['low_imsize'] is not None:
            low_imsize=o['low_imsize'] # allow over-ride
        else:
            low_imsize=o['imsize']*o['cellsize']/o['low_cell']
        # if mask-low exists then use it
        if os.path.isfile('bootstrap-mask-low.fits') and low_imsize==o['bsimsize']:
            extmask='bootstrap-mask-low.fits'
        else:
            extmask=None
        ddf_image('image_full_low',o['full_mslist'],cleanmask=extmask,cleanmode='SSD',ddsols=ddsols,applysols='AP',majorcycles=2,robust=o['low_robust'],uvrange=low_uvrange,beamsize=o['low_psf_arcsec'],imsize=low_imsize,cellsize=o['low_cell'],peakfactor=0.001,smooth=True,automask=True,automask_threshold=5,normalization=o['normalize'][2],colname=colname,catcher=catcher)
        make_mask('image_full_low.app.restored.fits',3.0,external_mask=extmask,catcher=catcher)
        ddf_image('image_full_low_im',o['full_mslist'],cleanmask='image_full_low.app.restored.fits.mask.fits',cleanmode='SSD',ddsols=ddsols,applysols='AP',majorcycles=1,robust=o['low_robust'],uvrange=low_uvrange,beamsize=o['low_psf_arcsec'],imsize=low_imsize,cellsize=o['low_cell'],peakfactor=0.001,smooth=True,automask=True,automask_threshold=5,normalization=o['normalize'][2],colname=colname,reuse_psf=True,dirty_from_resid=True,use_dicomodel=True,dicomodel_base='image_full_low',catcher=catcher)
        if o['restart'] and os.path.isfile('full-mask-low.fits'):
            warn('Full-bw mask exists, not making it')
        else:
            report('Making the full-bw extended source mask')
            make_extended_mask('image_full_low_im.app.restored.fits','image_dirin_SSD.app.restored.fits',rmsthresh=o['extended_rms'],sizethresh=1500,rootname='full',rmsfacet=o['rmsfacet'])
        extmask='full-mask-low.fits'
        make_mask('image_full_low_im.app.restored.fits',3.0,external_mask=extmask,catcher=catcher)
        ddf_image('image_full_low_m',o['full_mslist'],cleanmask='image_full_low_im.app.restored.fits.mask.fits',cleanmode='SSD',ddsols=ddsols,applysols='AP',majorcycles=1,robust=o['low_robust'],uvrange=low_uvrange,beamsize=o['low_psf_arcsec'],imsize=low_imsize,cellsize=o['low_cell'],peakfactor=0.001,smooth=True,automask=True,automask_threshold=4,normalization=o['normalize'][2],colname=colname,reuse_psf=True,dirty_from_resid=True,use_dicomodel=True,dicomodel_base='image_full_low_im',catcher=catcher,rms_factor=o['final_rmsfactor'])
        external_mask='external_mask_ext-deep.fits'
        if os.path.isfile(external_mask):
            warn('Deep external mask already exists, skipping creation')
        else:
            report('Make deep external mask')
            make_external_mask(external_mask,'image_dirin_SSD_init.dirty.fits',use_tgss=True,clobber=False,extended_use='full-mask-high.fits')

    # make mask from the previous run, will use new external mask if it exists
    make_mask('image_full_ampphase1.app.restored.fits',o['thresholds'][2],external_mask=external_mask,catcher=catcher)
    mask_dicomodel('image_full_ampphase1.DicoModel','image_full_ampphase1.app.restored.fits.mask.fits','image_ampphase1_full_masked.DicoModel',catcher=catcher)

    if o['exitafter'] == 'fulllow':
        warn('User specified exit after image_full_low.')
        sys.exit(2)

        # # before starting the final image, run the download thread if needed
        # if o['method'] is not None:
        #     report('Checking if optical catalogue download is required')
        #     from get_cat import get_cat, download_required
        #     if download_required(o['method']):
        #         download_thread = threading.Thread(target=get_cat, args=('panstarrs',))
        #         download_thread.start()
        #     else:
        #         warn('All data present, skipping download')
        #         download_thread = None

    # final image
    separator("DD Imaging (full mslist)")
    ddf_kw={}
    if o['final_psf_arcsec'] is not None:
        ddf_kw['beamsize']=o['final_psf_arcsec']
        if o['final_psf_minor_arcsec'] is not None:
            if o['final_psf_pa_deg'] is None:
                die('If final minor axis is supplied, position angle must be supplied too')
            ddf_kw['beamsize_minor']=o['final_psf_minor_arcsec']
            ddf_kw['beamsize_pa']=o['final_psf_pa_deg']

    ddf_image('image_full_ampphase2',o['full_mslist'],cleanmask='image_full_ampphase1.app.restored.fits.mask.fits',
              cleanmode='SSD',
              ddsols=ddsols,applysols='AP',
              majorcycles=1,robust=o['final_robust'],colname=colname,use_dicomodel=True,dicomodel_base='image_full_ampphase1_masked',
              peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][3],smooth=True,
              normalization=o['normalize'][2],uvrange=uvrange,
              apply_weights=o['apply_weights'][3],catcher=catcher,**ddf_kw)
    make_mask('image_full_ampphase2.app.restored.fits',o['thresholds'][3],external_mask=external_mask,catcher=catcher)
    mask_dicomodel('image_full_ampphase2.DicoModel','image_full_ampphase2.app.restored.fits.mask.fits','image_full_ampphase2_masked.DicoModel',catcher=catcher)

    if o['do_dynspec']:
        ddf_kw['predict_column']='PREDICT_DATA'

    separator("Deeper DD deconv (full mslist)")
    ddf_image('image_full_ampphase2m',o['full_mslist'],cleanmask='image_full_ampphase2.app.restored.fits.mask.fits',cleanmode='SSD',ddsols=ddsols,applysols='AP',majorcycles=1,robust=o['final_robust'],colname=colname,use_dicomodel=True,dicomodel_base='image_full_ampphase2_masked',peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][3],smooth=True,normalization=o['normalize'][2],reuse_psf=True,dirty_from_resid=True,uvrange=uvrange,apply_weights=o['apply_weights'][3],catcher=catcher,rms_factor=o['final_rmsfactor'],**ddf_kw)

    LastImage="image_full_ampphase2m.app.restored.fits"

        
    #     if o['second_selfcal']:
    #         ddsols="killms_f_ap2"
    #         us_ddsols=ddsols
    #         LastImage="image_full_ampphase2.app.restored.fits"
    #         if not os.path.exists('image_full_ampphase1m.Norm.fits'):
    #             os.symlink('image_full_ampphase1.Norm.fits','image_full_ampphase1m.Norm.fits')
    #         if o['auto_uvmin']:
    #             killms_uvrange[0]=optimize_uvmin('image_full_ampphase1m',o['mslist'],colname,o['solutions_uvmin'])
    #         make_mask('image_full_ampphase1m.app.restored.fits',o['thresholds'][3],external_mask=external_mask,catcher=catcher)
    #         mask_dicomodel('image_full_ampphase1m.DicoModel','image_full_ampphase1m.app.restored.fits.mask.fits','image_full_ampphase1m_masked.DicoModel',catcher=catcher)
    #         killms_data('image_full_ampphase1m',o['full_mslist'],'killms_f_ap2',colname=colname,clusterfile='image_dirin_SSD.npy.ClusterCat.npy',dicomodel='image_full_ampphase1m_masked.DicoModel',niterkf=o['NIterKF'][2],dt=o['final_dt'],catcher=catcher)
    #         ddf_image('image_full_ampphase2',o['full_mslist'],cleanmask='image_full_ampphase1m.app.restored.fits.mask.fits',cleanmode='SSD',ddsols='killms_f_ap2',applysols='AP',majorcycles=1,robust=o['final_robust'],colname=colname,use_dicomodel=True,dicomodel_base='image_full_ampphase1m_masked',peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][3],smooth=True,uvrange=uvrange,apply_weights=o['apply_weights'][3],catcher=catcher,rms_factor=o['final_rmsfactor'],**ddf_kw)

    #     if o['do_dynspec']:
    #         g=glob.glob('DynSpecs_*')
    #         if len(g)>0:
    #             warn('DynSpecs results directory %s already exists, skipping DynSpecs' % g[0])
    #         else:
    #             runcommand="ms2dynspec.py --ms big-mslist.txt --data SCALED_DATA --model PREDICT_DATA --sols killMS.%s.sols.npz --rad 2. --image %s --LogBoring %i"%(us_ddsols,LastImage,o['nobar'])
    #             run(runcommand,dryrun=o['dryrun'],log=logfilename('ms2dynspec.log'),quiet=o['quiet'])
            
    #     if o['method'] is not None:
    #         # have we got the catalogue?
    #         if download_thread is not None and download_thread.isAlive():
    #             warn('Waiting for background download thread to finish...')
    #             download_thread.join()
    #         # maybe the thread died, check the files are there
    #         if download_required(o['method']):
    #             warn('Retrying download for some or all of the catalogue')
    #             get_cat(o['method'])

    #         facet_offset_file='facet-offset.txt'
    #         if o['restart'] and os.path.isfile(facet_offset_file):
    #             warn('Offset file already exists, not running offsets.py')
    #         else:
    #             run('offsets.py '+' '.join(sys.argv[1:]),log=None)

    #         last_image_root='image_full_ampphase1m'
    #         if o['second_selfcal']:
    #             last_image_root='image_full_ampphase2'

    #         # check for LastResidual in cache. In case of a restart,
    #         # this may not be present, in which case we have to
    #         # remake.
    #         cachedir=find_cache_dir(o)
    #         full_mslist_file = os.path.basename(o['full_mslist'])
    #         if not(os.path.isfile(cachedir+'/'+full_mslist_file+'.ddfcache/LastResidual')) or not(os.path.isfile(cachedir+'/'+full_mslist_file+'.ddfcache/PSF')):
    #             ddf_image('image_full_ampphase1m_reimage',full_mslist_file,cleanmask='image_full_ampphase1.app.restored.fits.mask.fits',cleanmode='SSD',ddsols=ddsols,applysols='AP',majorcycles=0,robust=o['final_robust'],colname=colname,use_dicomodel=True,dicomodel_base='image_full_ampphase1m',peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][3],smooth=True,normalization=o['normalize'][2],reuse_psf=False,dirty_from_resid=False,uvrange=uvrange,apply_weights=o['apply_weights'][3],catcher=catcher,**ddf_kw)
    #             os.symlink('Dirty',cachedir+'/'+full_mslist_file+'.ddfcache/LastResidual')
    #             os.symlink('Dirty.hash',cachedir+'/'+full_mslist_file+'LastResidual.hash')

    #         ddf_shift(last_image_root,facet_offset_file,options=o,catcher=catcher)

            
    # # we got to the end, write a summary file
    
    # summary(o)
    # if o['clearcache_end']:
    #     full_clearcache(o)


if __name__=='__main__':
    # Main loop
    report('Welcome to ddf-pipeline, version '+__version__)
    if len(sys.argv)<2:
        warn('pipeline.py must be called with at least one parameter file or a command-line\noption list.\nE.g "pipeline.py example.cfg second_example.cfg --solutions-robust=0.1"\nSee below for a complete list of possible options with their default values.')
        print_options(option_list)
        sys.exit(1)

    o=options(sys.argv[1:],option_list)
    MyPickle.Save(o, "ddf-pipeline.last")

    main(o)
    

    
