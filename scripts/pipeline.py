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

import sys,os
if "PYTHONPATH_FIRST" in os.environ.keys() and int(os.environ["PYTHONPATH_FIRST"]):
    sys.path = os.environ["PYTHONPATH"].split(":") + sys.path
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

def get_solutions_timerange(sols):
    t = np.load(sols)['BeamTimes']
    return np.min(t),np.max(t)

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

    runcommand='DDF.py '+imagename+'.parset --Output-Name='+imagename+'_shift --Output-Mode=RestoreAndShift --Output-ShiftFacetsFile='+shiftfile+' --Predict-InitDicoModel '+imagename+'.DicoModel --Cache-SmoothBeam=force --Log-Memory 1 --Cache-Dir='+cache_dir
    
    fname=imagename+'_shift.app.facetRestored.fits'
    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF-shift step')
        if verbose:
            print 'would have run',runcommand
    else:
         run(runcommand,dryrun=options['dryrun'],log=logfilename('DDF-'+imagename+'_shift.log',options=options),quiet=options['quiet'])


def ddf_image(imagename,mslist,cleanmask=None,cleanmode='HMP',ddsols=None,applysols=None,threshold=None,majorcycles=3,use_dicomodel=False,robust=0,beamsize=None,beamsize_minor=None,beamsize_pa=None,reuse_psf=False,reuse_dirty=False,verbose=False,saveimages=None,imsize=None,cellsize=None,uvrange=None,colname='CORRECTED_DATA',peakfactor=0.1,dicomodel_base=None,options=None,do_decorr=None,normalization=None,dirty_from_resid=False,clusterfile=None,HMPsize=None,automask=True,automask_threshold=10.0,smooth=False,noweights=False,cubemode=False,apply_weights=True,catcher=None,rms_factor=3.0,predict_column=None,conditional_clearcache=False,PredictSettings=None,RMSFactorInitHMP=1.,MaxMinorIterInitHMP=10000,OuterSpaceTh=None,AllowNegativeInitHMP=False):

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

    if PredictSettings is not None and PredictSettings[0]=="Predict":
        fname="_has_predicted_OK.%s.info"%imagename

    runcommand = "DDF.py --Output-Name=%s --Data-MS=%s --Deconv-PeakFactor %f --Data-ColName %s --Parallel-NCPU=%i --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=%s --Deconv-Mode %s --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust %f --Image-NPix=%i --CF-wmax 50000 --CF-Nw 100 --Output-Also %s --Image-Cell %f --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=%f --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=%s --Log-Memory 1"%(imagename,mslist,peakfactor,colname,options['NCPU_DDF'],majorcycles,cleanmode,robust,imsize,saveimages,float(cellsize),rms_factor,cache_dir)

    runcommand += " --GAClean-RMSFactorInitHMP %f"%RMSFactorInitHMP
    runcommand += " --GAClean-MaxMinorIterInitHMP %f"%MaxMinorIterInitHMP
    if AllowNegativeInitHMP:
        runcommand += " --GAClean-AllowNegativeInitHMP True"
    if OuterSpaceTh is not None:
        runcommand += " --HMP-OuterSpaceTh %f"%OuterSpaceTh
        
    runcommand+=' --DDESolutions-SolsDir=%s'%options["SolsDir"]
    runcommand+=' --Cache-Weight=reset'

    
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
    if dicomodel_base is None and use_dicomodel:
        raise RuntimeError('that s wrong')
        
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
    return imagename
        
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


def clusterGA(imagename="image_dirin_SSD_m.app.restored.fits",OutClusterCat=None,options=None,use_makemask_products=False):

    if os.path.isfile(OutClusterCat):
        warn('File %s already exists, skipping clustering step'%OutClusterCat)
        return

    if not ".app.restored.fits" in imagename:
        raise RuntimeError('Input image should be an apparant restored image')

    if options is None:
        options=o # attempt to get global if it exists

    if use_makemask_products:
        runcommand="MakeCatalog.py --RestoredIm %s --rmsmean_map [Noise.mean.fits,Noise.fits]"%imagename
    else:
        runcommand="MakeCatalog.py --RestoredIm %s"%imagename 
    run(runcommand,dryrun=options['dryrun'],log=logfilename('MakeCatalog-'+imagename+'.log',options=options),quiet=options['quiet'])

    Name=imagename.split(".app.restored.fits")[0]

    #runcommand="ClusterCat.py --SourceCat %s.app.restored.pybdsm.srl.fits --AvoidPolygons MaskDiffuse.pickle --NGen 100 --FluxMin 0.1"%Name
    if use_makemask_products:
        runcommand="ClusterCat.py --SourceCat %s.app.restored.pybdsm.srl.fits --AvoidPolygons MaskDiffuse.pickle --DoPlot=0 --NGen 100 --NCPU %i"%(Name,options['NCPU_DDF'])
    else:
        runcommand="ClusterCat.py --SourceCat %s.app.restored.pybdsm.srl.fits --DoPlot=0 --NGen 100 --NCPU %i"%(Name,options['NCPU_DDF'])
    if OutClusterCat is not None:
        runcommand+=" --OutClusterCat %s"%OutClusterCat
    runcommand+=" --NCluster %i"%o['ndir']
    run(runcommand,dryrun=options['dryrun'],log=logfilename('MakeCluster-'+imagename+'.log',options=options),quiet=options['quiet'])


def make_mask(imagename,thresh,verbose=False,options=None,external_mask=None,catcher=None,OutMaskExtended=None):
    if catcher: catcher.check()

    # mask_use specifies a mask file to use
    if options is None:
        options=o # attempt to get global

    if options['dryrun']: return
    fname=imagename+'.mask.fits'
    runcommand = "MakeMask.py --RestoredIm=%s --Th=%s --Box=50,2"%(imagename,thresh)
    if OutMaskExtended is not None:
        runcommand += " --OutMaskExtended %s --OutNameNoiseMap Noise"%(OutMaskExtended)


        
    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeMask step')
        if verbose:
            print 'Would have run',runcommand
    else:
        run(runcommand,dryrun=options['dryrun'],log=logfilename('MM-'+imagename+'.log',options=options),quiet=options['quiet'])
        if external_mask is not None:
            if type(external_mask) is str:
                merge_mask(fname,external_mask,fname)
            else:
                for mask in external_mask:
                    merge_mask(fname,mask,fname)
    return fname
            


def killms_data(imagename,mslist,outsols,clusterfile=None,colname='CORRECTED_DATA',niterkf=6,dicomodel=None,
                uvrange=None,wtuv=None,robust=None,catcher=None,dt=None,options=None,
                SolverType="KAFCA",PolMode="Scalar",MergeSmooth=False,NChanSols=1,
                DISettings=None,EvolutionSolFile=None,CovQ=0.1,InterpToMSListFreqs=None):

    if options is None:
        options=o # attempt to get global if it exists

    cache_dir=find_cache_dir(options)

    # run killms individually on each MS -- allows restart if it failed in the middle
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    for f in filenames:
        if catcher: catcher.check()

        SolsDir=options["SolsDir"]
        fname=f
        if SolsDir is None or SolsDir=="":
            solname =fname+'/killMS.'+outsols+'.sols.npz'
        else:
            MSName=os.path.abspath(f).split("/")[-1]
            solname =os.path.abspath(SolsDir)+"/"+MSName+'/killMS.'+outsols+'.sols.npz'
        checkname=solname



        #checkname=f+'/killMS.'+outsols+'.sols.npz'
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

            runcommand+=' --SolsDir=%s'%options["SolsDir"]
            
                
            if DISettings is None:
                runcommand+='  --BeamMode LOFAR --LOFARBeamMode=A --DDFCacheDir=%s'%cache_dir
                if clusterfile is not None:
                    runcommand+=' --NodesFile '+clusterfile
                if dicomodel is not None:
                    runcommand+=' --DicoModel '+dicomodel
                if EvolutionSolFile is not None:
                    runcommand+=' --EvolutionSolFile '+EvolutionSolFile
                    
            else:
                runcommand+=" --SolverType %s --PolMode %s --SkyModelCol %s --OutCol %s --ApplyToDir 0"%DISettings
                _,_,ModelColName,_=DISettings
                _,dt,_,n_df=give_dt_dnu(f,
                                        DataCol=colname,
                                        ModelCol=ModelColName,
                                        T=10.)
                runcommand+=" --dt %f --NChanSols %i"%(dt+1e-4,n_df)
                
                
                
            rootfilename=outsols.split('/')[-1]
            f_=f.replace("/","_")
            run(runcommand,dryrun=o['dryrun'],log=logfilename('KillMS-'+f_+'_'+rootfilename+'.log'),quiet=o['quiet'])

            # Clip anyway - on IMAGING_WEIGHT by default
            if DISettings is not None:
                ClipCol=DISettings[-1]
            else:
                ClipCol=colname
            runcommand="ClipCal.py --MSName %s --ColName %s"%(f,ClipCol)
            run(runcommand,dryrun=o['dryrun'],log=logfilename('ClipCal-'+f_+'_'+rootfilename+'.log'),quiet=o['quiet'])

    if MergeSmooth:
        outsols=smooth_solutions(mslist,outsols,catcher=None,dryrun=o['dryrun'],InterpToMSListFreqs=InterpToMSListFreqs)



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

def smooth_solutions(mslist,ddsols,catcher=None,dryrun=False,InterpToMSListFreqs=None):
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    full_sollist = []
    start_times = []
    SolsDir=o["SolsDir"]
    if SolsDir is None or SolsDir=="":
        for fname in filenames:
            solname =fname+'/killMS.'+ddsols+'.sols.npz'
            t0,t1 = get_solutions_timerange(solname)
            start_times.append(t0)
            full_sollist.append(solname)
            f.write('%s\n'%(solname))
    else:
        for fname in filenames:
            MSName=os.path.abspath(fname).split("/")[-1]
            solname =os.path.abspath(SolsDir)+"/"+MSName+'/killMS.'+ddsols+'.sols.npz'
            t0,t1 = get_solutions_timerange(solname)
            start_times.append(t0)
            full_sollist.append(solname)

    Ustart_times = np.unique(start_times)

    for start_time in Ustart_times:
        with open('solslist_%s.txt'%start_time,'w') as f:
            for i in range(0,len(full_sollist)):
                if start_times[i] == start_time:
                    solname = full_sollist[i]
                    f.write('%s\n'%(solname))
        
        checkname='%s_%s_merged.npz'%(ddsols,start_time)
        if o['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running MergeSols step')
        else:
            run('MergeSols.py --SolsFilesIn=solslist_%s.txt --SolFileOut=%s_%s_merged.npz'%(start_time,ddsols,start_time),dryrun=dryrun)
        checkname='%s_%s_smoothed.npz'%(ddsols,start_time)
        if o['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running SmoothSols step')
        else:
            run('SmoothSols.py --SolsFileIn=%s_%s_merged.npz --SolsFileOut=%s_%s_smoothed.npz --InterpMode=TEC,PolyAmp'%(ddsols,start_time,ddsols,start_time),dryrun=dryrun)

        smoothoutname='%s_%s_smoothed.npz'%(ddsols,start_time)
        if InterpToMSListFreqs:
            interp_outname="%s_%s_interp.npz"%(smoothoutname,start_time)
            checkname=interp_outname
            if o['restart'] and os.path.isfile(checkname):
                warn('Solutions file '+checkname+' already exists, not running MergeSols step')
            else:
                command="InterpSols.py --SolsFileIn %s --SolsFileOut %s --MSOutFreq %s"%(smoothoutname,interp_outname,InterpToMSListFreqs)
                run(command,dryrun=dryrun)
        
        for i in range(0,len(full_sollist)):
            if start_times[i] == start_time:
		symsolname = full_sollist[i].replace(ddsols,ddsols+'_smoothed')
                if o['restart'] and os.path.isfile(symsolname):
	            warn('Symlink ' + symsolname + ' already exists')
		else:
	            warn('Symlink ' + symsolname + ' does not exist -- creating')
                    os.symlink(os.path.abspath('%s_%s_smoothed.npz'%(ddsols,start_time)),symsolname)
        outname = ddsols + '_smoothed'

    return outname

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

def give_dt_dnu(msname,DataCol="DATA",ModelCol="DI_PREDICT",T=10.):
    t=pt.table(msname,ack=False)
    d=t.getcol(DataCol)
    dt_bin_sec=t.getcol("INTERVAL",0,1,1)[0]
    _,nch,_=d.shape
    f=t.getcol("FLAG")
    p=t.getcol(ModelCol)
    t.close()
    fp=f[:,:,np.array([1,2])]
    dp=d[:,:,np.array([1,2])]
    dps=dp[fp==0]
    da=np.abs(d[:,:,0][f[:,:,0]==0])
    S=np.std(dps)
    M=np.mean(da)
    nb=T**2/(M/S)**2

    # find the size of the channel step  
    nch_step=int(round(np.sqrt(nb)))
    nch_step=np.max([1,nch_step])
    nch_step=np.min([nch,nch_step])
    warn('nch_step=%i'%(nch_step))

    # find the step to have equal interval size
    #nch_bin=int(nch/nch_step)+1
    #nch_step=int(nch/float(nch_bin))
    lDiv=np.array([i for i in range(1,nch+1) if nch%i==0])
    inch=np.argmin(np.abs(lDiv-nch_step))
    nch_step=lDiv[inch]
    nch_step=np.max([1,nch_step])
    nch_step=np.min([nch,nch_step])
    
    nt_step=int(round(nb/float(nch_step)))
    nt_step=np.max([1,nt_step])

    SNR=np.sqrt(nt_step*nch_step)*M/S
    warn('Using (dt,df)=(%i,%i) for CubiCal run of %s with (<|model|>,std)=(%.2f,%.2f) giving SNR=%.2f'%(nt_step,nch_step,msname,M,S,SNR))
    
    return nt_step, nt_step*dt_bin_sec/60.0, nch_step, nch/nch_step
    
def cubical_data(mslist,
                 NameSol="DI0",
                 n_dt=1,
                 n_df=2,
                 n_DT=None,
                 DataColName="DATA",
                 ModelColName="DD_PREDICT",
                 OutColName="DATA_DI_CORRECTED",
                 options=None,
                 ReinitWeights=False):
    if n_DT is None:
        n_DT=10*n_dt
        
    if options is None:
        options=o # attempt to get global if it exists

    
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    for f in filenames:
        ThisMSName=os.path.abspath(f)
        SolsDir=options["SolsDir"]
        
        MSName=ThisMSName.split("/")[-1]
        if SolsDir is None or SolsDir=="":
            solname ="%s/CubiCal_%s"%(MSName,NameSol)
        else:
            DirName=os.path.abspath(SolsDir)+"/"+MSName
            solname =os.path.abspath(SolsDir)+"/"+MSName+'/CubiCal_%s'%NameSol
            if not os.path.isdir(DirName):
                os.makedirs(DirName)
        checkname="%s.noise.antchan.png"%solname

        if o['restart'] and os.path.isfile(checkname):
            warn('File '+checkname+' already exists, not running CubiCal step')
            continue

        n_dt,_,n_df,_=give_dt_dnu(ThisMSName,
                                DataCol=DataColName,
                                ModelCol=ModelColName,
                                T=10.)

        n_DT=10*n_dt

        command="gocubical --data-ms %s --out-mode sc --g-time-int %i --g-freq-int %i --data-time-chunk %i --data-freq-chunk 0 --data-column %s --model-list %s --out-column %s --dist-ncpu %i --weight-column None --out-casa-gaintables 0 --flags-reinit-bitflag 1 --flags-save None --out-name %s --g-max-prior-error 1e10 --g-max-post-error 1e10"%(ThisMSName,n_dt,n_df,n_DT,DataColName,ModelColName,OutColName,o['NCPU_DDF'],solname)

        run(command,dryrun=o['dryrun'])#,log=logfilename('CubiCal-'+f_+'_'+rootfilename+'.log'),quiet=o['quiet'])

        runcommand="ClipCal.py --MSName %s --ColName %s"%(ThisMSName,OutColName)
        if ReinitWeights:
            runcommand+=" --ReinitWeights 1"
            
        run(runcommand,dryrun=o['dryrun'])#,log=logfilename('ClipCal-'+f_+'_'+rootfilename+'.log'),quiet=o['quiet'])
        

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
    separator("DI Deconv (externally defined sources)")
    CurrentBaseDicoModelName=ddf_image('image_dirin_SSD',o['mslist'],cleanmask=external_mask,cleanmode='SSD',
                                       majorcycles=1,robust=o['image_robust'],reuse_psf=True,reuse_dirty=True,
                                       peakfactor=0.01,rms_factor=3,
                                       colname=colname,clusterfile=None,automask=True,
                                       automask_threshold=o['thresholds'][0],apply_weights=o['apply_weights'][0],
                                       uvrange=uvrange,catcher=catcher)

    
    separator("Make the diffuse emission mask")
    # Make the diffuse emission mask
    _=make_mask('image_dirin_SSD.residual01.fits',
                o['thresholds'][0],
                external_mask=external_mask,
                catcher=catcher,
                OutMaskExtended="MaskDiffuse")
    separator("Merge diffuse emission mask into external mask")
    merge_mask(external_mask,"MaskDiffuse.fits",external_mask)

    # make a mask from the final image
    separator("Make mask for next iteration")
    CurrentMaskName=make_mask('image_dirin_SSD.app.restored.fits',
                              o['thresholds'][0],
                              external_mask=external_mask,
                              catcher=catcher)
    
    
    separator("Continue deconvolution")
    CurrentBaseDicoModelName=ddf_image('image_dirin_SSD_m',o['mslist'],
                                       cleanmask=CurrentMaskName,cleanmode='SSD',
                                       majorcycles=2,robust=o['image_robust'],
                                       reuse_psf=True,
                                       dicomodel_base=CurrentBaseDicoModelName,
                                       use_dicomodel=True,
                                       dirty_from_resid=True,
                                       peakfactor=0.001,rms_factor=0,
                                       colname=colname,clusterfile=None,
                                       automask=True,
                                       automask_threshold=o['thresholds'][0],apply_weights=o['apply_weights'][0],
                                       uvrange=uvrange,catcher=catcher,
                                       RMSFactorInitHMP=1.,
                                       MaxMinorIterInitHMP=10000,
                                       PredictSettings=("Clean","DD_PREDICT"))



    #########################
    if o['clusterfile'] is None:
        separator("Cluster the sky model")
        ClusterFile='image_dirin_SSD_m.npy.ClusterCat.npy'
        clusterGA(imagename="image_dirin_SSD_m.app.restored.fits",
                  OutClusterCat=ClusterFile,
                  use_makemask_products=True)
    else:
        ClusterFile=o['clusterfile']
        warn('Using user-specifed cluster file '+ClusterFile)

    #########################
    clearcache(o['mslist'],o)
    separator("Deconv clustered DI image")
    CurrentBaseDicoModelName=ddf_image('image_dirin_SSD_m_c',o['mslist'],
                                       cleanmask=CurrentMaskName,
                                       cleanmode='SSD',
                                       majorcycles=1,robust=o['image_robust'],
                                       #reuse_psf=True,
                                       #reuse_dirty=True,
                                       dicomodel_base=CurrentBaseDicoModelName,
                                       use_dicomodel=True,
                                       #dirty_from_resid=True,
                                       peakfactor=0.001,rms_factor=0,
                                       colname=colname,
                                       clusterfile=ClusterFile,
                                       automask=True,
                                       automask_threshold=o['thresholds'][0],
                                       apply_weights=o['apply_weights'][0],
                                       uvrange=uvrange,catcher=catcher,
                                       RMSFactorInitHMP=1.,
                                       MaxMinorIterInitHMP=10000,
                                       PredictSettings=("Clean","DD_PREDICT"))

    if o['exitafter'] == 'dirin':
        warn('User specified exit after image_dirin.')
        sys.exit(2)

        separator("DI CAL")
    ########################
    killms_data('PredictDI_0',o['mslist'],'DIS0',colname=colname,
                dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                niterkf=o['NIterKF'][0],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                catcher=catcher,
                dt=o['dt_di'],
                NChanSols=o['NChanSols_di'],
                DISettings=("CohJones","IFull","DD_PREDICT","DATA_DI_CORRECTED"))
    # cubical_data(o['mslist'],
    #              NameSol="DIS0",
    #              n_dt=1,
    #              n_df=2,
    #              n_DT=None,
    #              DataColName=colname,
    #              ModelColName="DD_PREDICT",
    #              OutColName="DATA_DI_CORRECTED",
    #              ReinitWeights=True)
    
    colname="DATA_DI_CORRECTED"


    _=ddf_image('image_dirin_SSD_m_c_di',o['mslist'],
                cleanmask=CurrentMaskName,cleanmode='SSD',
                majorcycles=0,robust=o['image_robust'],
                #reuse_psf=True,
                dicomodel_base=CurrentBaseDicoModelName,
                use_dicomodel=True,
                #dirty_from_resid=True,
                peakfactor=0.001,rms_factor=0,
                colname=colname,clusterfile=None,
                automask=True,
                automask_threshold=o['thresholds'][0],
                apply_weights=True,#o['apply_weights'][0],
                uvrange=uvrange,catcher=catcher,
                RMSFactorInitHMP=1.,
                MaxMinorIterInitHMP=10000,
                PredictSettings=("Clean","DD_PREDICT"))


    CurrentBaseDicoModelName=ddf_image('image_dirin_SSD_m_c_di_m',o['mslist'],
                                       cleanmask=CurrentMaskName,cleanmode='SSD',
                                       majorcycles=1,robust=o['image_robust'],
                                       reuse_psf=True,
                                       reuse_dirty=True,
                                       dicomodel_base=CurrentBaseDicoModelName,
                                       use_dicomodel=True,
                                       #dirty_from_resid=True,
                                       peakfactor=0.001,rms_factor=0,
                                       colname=colname,clusterfile=None,
                                       automask=True,
                                       automask_threshold=o['thresholds'][0],
                                       apply_weights=True,#o['apply_weights'][0],
                                       uvrange=uvrange,catcher=catcher,
                                       RMSFactorInitHMP=1.,
                                       MaxMinorIterInitHMP=10000,
                                       PredictSettings=("Clean","DD_PREDICT"))


    if o['exitafter'] == 'dirin_di':
        warn('User specified exit after image_dirin with DI calibration.')
        sys.exit(2)


    separator("DD calibration")
    CurrentDDkMSSolName=killms_data('image_dirin_SSD_m_c_di_m',o['mslist'],'DDS0',colname=colname,
                                    dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                    clusterfile=ClusterFile,
                                    CovQ=0.02,
                                    niterkf=o['NIterKF'][0],
                                    #CovQ=0.1,
                                    #niterkf=6,
                                    uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],dt=o['dt_slow'],
                                    catcher=catcher,NChanSols=o['NChanSols'],
                                    MergeSmooth=True)

    # ##########################################################
    # run bootstrap, and change the column name if it runs
    if o['bootstrap']:
        separator("Bootstrap")
        report('Running bootstrap')
        run('bootstrap.py '+' '.join(sys.argv[1:]),log=None,dryrun=o["dryrun"])
        colname='SCALED_DATA'
        if o['exitafter'] == 'bootstrap':
            warn('User specified exit after phase-only deconvolution.')
            sys.exit(2)

    # make a mask from the full-res image
    separator("Make mask for next iteration")
    CurrentMaskName=make_mask('image_dirin_SSD_m_c_di_m.app.restored.fits',
                              o['thresholds'][1],
                              external_mask=external_mask,
                              catcher=catcher)

    separator("PhaseOnly deconv")
    CurrentBaseDicoModelName=ddf_image('image_phase1',o['mslist'],
                                       cleanmask=CurrentMaskName,
                                       cleanmode='SSD',
                                       ddsols=CurrentDDkMSSolName,applysols='P',majorcycles=2,robust=o['image_robust'],
                                       colname=colname,peakfactor=0.001,automask=True,
                                       automask_threshold=o['thresholds'][1],
                                       normalization=o['normalize'][0],apply_weights=o['apply_weights'][1],uvrange=uvrange,
                                       use_dicomodel=True,
                                       dicomodel_base=CurrentBaseDicoModelName,
                                       catcher=catcher,
                                       RMSFactorInitHMP=1.,
                                       MaxMinorIterInitHMP=10000,
                                       PredictSettings=("Clean","DD_PREDICT"))

    if o['exitafter'] == 'phase':
        warn('User specified exit after phase-only deconvolution.')
        sys.exit(2)

    separator("Mask for deeper deconv")
    CurrentMaskName=make_mask('image_phase1.app.restored.fits',o['thresholds'][1],external_mask=external_mask,catcher=catcher)
    CurrentBaseDicoModelName=mask_dicomodel('image_phase1.DicoModel',CurrentMaskName,'image_phase1_masked.DicoModel',catcher=catcher)

    separator("DD calibration")
    CurrentDDkMSSolName=killms_data('image_phase1',o['mslist'],'DDS1',colname=colname,
                                    dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                    CovQ=0.02,
                                    clusterfile=ClusterFile,
                                    niterkf=o['NIterKF'][0],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                                    dt=o['dt_slow'],
                                    catcher=catcher,NChanSols=o['NChanSols'],
                                    EvolutionSolFile=CurrentDDkMSSolName,
                                    MergeSmooth=True)
    ##############################################

    separator("AmpPhase deconv")
    CurrentBaseDicoModelName=ddf_image('image_ampphase1',o['mslist'],
                                       cleanmask=CurrentMaskName,cleanmode='SSD',
                                       ddsols=CurrentDDkMSSolName,applysols='AP',majorcycles=1,robust=o['image_robust'],
                                       colname=colname,peakfactor=0.001,automask=True,
                                       automask_threshold=o['thresholds'][1],
                                       normalization=o['normalize'][0],apply_weights=o['apply_weights'][1],uvrange=uvrange,
                                       use_dicomodel=True,
                                       dicomodel_base=CurrentBaseDicoModelName,
                                       catcher=catcher,
                                       RMSFactorInitHMP=1.,
                                       #AllowNegativeInitHMP=True,
                                       MaxMinorIterInitHMP=10000,
                                       PredictSettings=("Clean","DD_PREDICT"))

    if o['exitafter'] == 'ampphase':
        warn('User specified exit after amp-phase deconvolution.')
        sys.exit(2)

    separator("Make Mask")
    CurrentMaskName=make_mask('image_ampphase1.app.restored.fits',o['thresholds'][1],external_mask=external_mask,catcher=catcher)
    CurrentBaseDicoModelName=mask_dicomodel('image_ampphase1.DicoModel',CurrentMaskName,'image_ampphase1m_masked.DicoModel',catcher=catcher)

    separator("Second DI calibration")
    ddf_image('Predict_DI1',o['mslist'],
              cleanmask=CurrentMaskName,cleanmode='SSD',
              ddsols=CurrentDDkMSSolName,applysols='AP',majorcycles=1,robust=o['image_robust'],
              colname=colname,peakfactor=0.001,automask=True,
              automask_threshold=o['thresholds'][1],
              normalization=o['normalize'][0],apply_weights=o['apply_weights'][1],uvrange=uvrange,
              use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              catcher=catcher,
              RMSFactorInitHMP=1.,
              #AllowNegativeInitHMP=True,
              MaxMinorIterInitHMP=10000,
              PredictSettings=("Predict","DD_PREDICT"))

    separator("Another DI step")
    killms_data('PredictDI_1',o['mslist'],'DIS1',colname=o['colname'],
                dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                #clusterfile=ClusterFile,
                niterkf=o['NIterKF'][0],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                catcher=catcher,
                dt=o['dt_di'],
                NChanSols=o['NChanSols_di'],
                DISettings=("CohJones","IFull","DD_PREDICT","DATA_DI_CORRECTED"))
    # cubical_data(o['mslist'],
    #              NameSol="DIS1",
    #              n_dt=1,
    #              n_df=2,
    #              n_DT=None,
    #              DataColName=o['colname'],
    #              ModelColName="DD_PREDICT",
    #              OutColName="DATA_DI_CORRECTED")

    CurrentBaseDicoModelName=ddf_image('image_ampphase1_di',o['mslist'],
                                       cleanmask=CurrentMaskName,cleanmode='SSD',
                                       ddsols=CurrentDDkMSSolName,applysols='AP',majorcycles=1,robust=o['image_robust'],
                                       colname=colname,peakfactor=0.001,automask=True,
                                       automask_threshold=o['thresholds'][1],
                                       normalization=o['normalize'][0],apply_weights=o['apply_weights'][1],uvrange=uvrange,
                                       use_dicomodel=True,
                                       dicomodel_base=CurrentBaseDicoModelName,
                                       catcher=catcher,
                                       RMSFactorInitHMP=1.,
                                       #AllowNegativeInitHMP=True,
                                       MaxMinorIterInitHMP=10000,
                                       PredictSettings=("Clean","DD_PREDICT"))

    if o['exitafter'] == 'ampphase_di':
        warn('User specified exit after amp-phase plus DI deconvolution.')
        sys.exit(2)

    if o['full_mslist'] is None:
        warn('No full mslist provided, stopping here')
        summary(o)
        sys.exit(3)
        
    separator("DD calibration of full mslist")
    
    CurrentDDkMSSolName=killms_data('image_ampphase1_di',o['full_mslist'],'DDS2_full',
                                    colname=o['colname'],
                                    dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                    CovQ=0.1,
                                    clusterfile=ClusterFile,
                                    niterkf=6,
                                    # o['NIterKF'][0],
                                    uvrange=killms_uvrange,
                                    wtuv=o['wtuv'],
                                    robust=o['solutions_robust'],
                                    dt=o['dt_slow'],
                                    catcher=catcher,
                                    NChanSols=o['NChanSols'],
                                    # EvolutionSolFile=CurrentDDkMSSolName,
                                    MergeSmooth=True)
    
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


    # #########################################################################
    # ###############                  BIG MSLIST               ###############
    # #########################################################################


    # check full mslist imaging weights
    check_imaging_weight(o['full_mslist'])
        
    # Compute the DD predict
    colname=o['colname']
    separator("Compute DD Predict (full mslist)")
    ddf_image('Predict_DDS2',o['full_mslist'],cleanmode='SSD',
              applysols='AP',majorcycles=1,robust=o['image_robust'],colname=colname,peakfactor=0.01,
              automask=True,automask_threshold=o['thresholds'][1],normalization=o['normalize'][0],
              apply_weights=o['apply_weights'][0],uvrange=uvrange,use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              catcher=catcher,
              ddsols=CurrentDDkMSSolName, PredictSettings=("Predict","DD_PREDICT"))

    separator("Compute DI calibration (full mslist)")
    # cubical_data(o['full_mslist'],
    #              NameSol="DIS2_full",
    #              n_dt=1,
    #              n_df=2,
    #              n_DT=None,
    #              DataColName=o['colname'],
    #              ModelColName="DD_PREDICT",
    #              OutColName="DATA_DI_CORRECTED")
    killms_data('Predict_DDS2',o['full_mslist'],'DIS2_full',colname=colname,
                dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                clusterfile=ClusterFile,
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
              applysols='AP',
              majorcycles=0,
              #robust=o['image_robust'],
              robust=o['final_robust'],
              colname=colname,use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              AllowNegativeInitHMP=True,
              peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][2],
              normalization=o['normalize'][1],uvrange=uvrange,smooth=True,
              apply_weights=o['apply_weights'][2],catcher=catcher)

    separator("MakeMask")
    CurrentMaskName=make_mask('image_full_ampphase_di.app.restored.fits',10,external_mask=external_mask,catcher=catcher)

    separator("Finish Deconvolution AP (full mslist)")
    CurrentBaseDicoModelName=ddf_image('image_full_ampphase_di_m',o['full_mslist'],
                                       cleanmask=CurrentMaskName,
                                       reuse_psf=True,
                                       reuse_dirty=True,
                                       cleanmode='SSD',ddsols=CurrentDDkMSSolName,
                                       applysols='AP',majorcycles=1,robust=o['final_robust'],
                                       colname=colname,use_dicomodel=True,
                                       dicomodel_base=CurrentBaseDicoModelName,
                                       peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][2],
                                       normalization=o['normalize'][1],uvrange=uvrange,
                                       apply_weights=o['apply_weights'][2],catcher=catcher,
                                       AllowNegativeInitHMP=True,
                                       RMSFactorInitHMP=.5,
                                       MaxMinorIterInitHMP=10000,smooth=True)

    if o['exitafter'] == 'ampphase':
        warn('User specified exit after image_ampphase.')
        sys.exit(2)

    separator("DD Calibration (full mslist)")
    CurrentDDkMSSolName=killms_data('image_full_ampphase_di_m',
                                    o['full_mslist'],'DDS3_full',
                                    colname=colname,
                                    clusterfile=ClusterFile,
                                    dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                    niterkf=6,
                                    CovQ=0.1,
                                    uvrange=killms_uvrange,
                                    wtuv=o['wtuv'],
                                    robust=o['solutions_robust'],
                                    MergeSmooth=True,
                                    dt=o['dt_fast'],catcher=catcher)#,EvolutionSolFile=CurrentDDkMSSolName)


    if o['low_psf_arcsec'] is not None:
        # low-res image requested
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

        ddf_image('image_full_low',o['full_mslist'],
                  cleanmask=extmask,
                  cleanmode='SSD',ddsols=CurrentDDkMSSolName,
                  applysols='AP',
		  AllowNegativeInitHMP=True,
                  majorcycles=2,robust=o['low_robust'],
                  colname=colname,use_dicomodel=False,
                  uvrange=low_uvrange,beamsize=o['low_psf_arcsec'],
                  imsize=low_imsize,cellsize=o['low_cell'],peakfactor=0.001,
                  smooth=True,automask=True,automask_threshold=5,normalization=o['normalize'][2],
                  catcher=catcher)

        make_mask('image_full_low.app.restored.fits',3.0,external_mask=extmask,catcher=catcher)

        ddf_image('image_full_low_im',o['full_mslist'],
              cleanmask='image_full_low.app.restored.fits.mask.fits',
              cleanmode='SSD',ddsols=CurrentDDkMSSolName,
              applysols='AP',
	      AllowNegativeInitHMP=True,
              majorcycles=1,robust=o['low_robust'],
              uvrange=low_uvrange,beamsize=o['low_psf_arcsec'],
              imsize=low_imsize,cellsize=o['low_cell'],peakfactor=0.001,
              smooth=True,automask=True,automask_threshold=5,normalization=o['normalize'][2],colname=colname,
              reuse_psf=True,dirty_from_resid=True,use_dicomodel=True,dicomodel_base='image_full_low',
              catcher=catcher)


        if o['restart'] and os.path.isfile('full-mask-low.fits'):
            warn('Full-bw mask exists, not making it')
        else:
            report('Making the full-bw extended source mask')
            make_extended_mask('image_full_low_im.app.restored.fits','image_dirin_SSD.app.restored.fits',rmsthresh=o['extended_rms'],sizethresh=1500,rootname='full',rmsfacet=o['rmsfacet'])
            report('Make_extended_mask returns')
            extmask='full-mask-low.fits'
            make_mask('image_full_low_im.app.restored.fits',3.0,external_mask=extmask,catcher=catcher)

        ddf_image('image_full_low_m',o['full_mslist'],
              cleanmask='image_full_low_im.app.restored.fits.mask.fits',
              cleanmode='SSD',ddsols=CurrentDDkMSSolName,
              applysols='AP',
	      AllowNegativeInitHMP=True,
              majorcycles=1,robust=o['low_robust'],
              uvrange=low_uvrange,beamsize=o['low_psf_arcsec'],
              imsize=low_imsize,cellsize=o['low_cell'],peakfactor=0.001,
              smooth=True,automask=True,automask_threshold=4,normalization=o['normalize'][2],colname=colname,
              reuse_psf=True,dirty_from_resid=True,use_dicomodel=True,dicomodel_base='image_full_low_im',
              catcher=catcher,rms_factor=o['final_rmsfactor'])
        external_mask='external_mask_ext-deep.fits'
        if os.path.isfile(external_mask):
            warn('Deep external mask already exists, skipping creation')
        else:
            report('Make deep external mask')
            make_external_mask(external_mask,'image_dirin_SSD_init.dirty.fits',use_tgss=True,clobber=False,extended_use='full-mask-high.fits')

    # ##########################################################
    if o['exitafter'] == 'fulllow':
        warn('User specified exit after full low.')
        sys.exit(2)



    separator("MakeMask")
    CurrentMaskName=make_mask('image_full_ampphase_di_m.app.restored.fits',o['thresholds'][2],external_mask=external_mask,catcher=catcher)
    CurrentBaseDicoModelName=mask_dicomodel('image_full_ampphase_di_m.DicoModel',CurrentMaskName,'image_full_ampphase_di_m_masked.DicoModel',catcher=catcher)
            
    # full resolution, one iter of deconvolution
    separator("DD imaging (full resolution)")
    ddf_kw={}
    if o['final_psf_arcsec'] is not None:
        ddf_kw['beamsize']=o['final_psf_arcsec']
        if o['final_psf_minor_arcsec'] is not None:
            if o['final_psf_pa_deg'] is None:
                die('If final minor axis is supplied, position angle must be supplied too')
            ddf_kw['beamsize_minor']=o['final_psf_minor_arcsec']
            ddf_kw['beamsize_pa']=o['final_psf_pa_deg']
    ddf_image('image_full_ampphase_di_m.NS',o['full_mslist'],
              cleanmask=CurrentMaskName,
              reuse_psf=False,
              cleanmode='SSD',ddsols=CurrentDDkMSSolName,
              applysols='AP',majorcycles=1,robust=o['final_robust'],
              colname=colname,use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              AllowNegativeInitHMP=True,
              peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][2],
              normalization=o['normalize'][1],uvrange=uvrange,smooth=True,
              apply_weights=o['apply_weights'][2],catcher=catcher,RMSFactorInitHMP=1.,**ddf_kw)

    separator('Write summary and tidy up')
    summary(o)
    if o['clearcache_end']:
         full_clearcache(o)

    return

    # ##########################################################
    # EVERYTHING BELOW HERE IS OLD CODE                        #
    # CONTAINS DYNSPEC AND ASTROMETRY CODE
    # DO NOT DELETE UNTIL MERGED BACK IN                       #
    # ##########################################################


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
    

    
