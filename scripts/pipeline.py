#!/usr/bin/env python
"""
ddf-pipeline, a pipeline for LOFAR data reduction
Copyright (C) 2017-2022 Martin Hardcastle (mjh@extragalactic.info) and others

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
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from future import standard_library
standard_library.install_aliases()
from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div
import sys,os
if "PYTHONPATH_FIRST" in list(os.environ.keys()) and int(os.environ["PYTHONPATH_FIRST"]):
    sys.path = os.environ["PYTHONPATH"].split(":") + sys.path
import os.path
from auxcodes import report,run,warn,die,Catcher,dotdict,separator,MSList
from parset import option_list
from options import options,print_options
from shutil import copyfile,rmtree,move
import glob
import pyrap.tables as pt
from redo_dppp_di import redo_dppp_di
from modify_mask import modify_mask
from make_extended_mask import make_extended_mask,merge_mask,add_manual_mask
from histmsamp import find_uvmin,sumdico
import numpy as np
from astropy.io import fits
from pipeline_version import version
__version__=version()
import datetime
import threading

try:
    from killMS.Other import MyPickle
except ImportError:
    MyPickle=None
from surveys_db import use_database,update_status,SurveysDB

def summary(o):
    with open('summary.txt','w') as f:
        ts='{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
        f.write('ddf-pipeline completed at '+ts+'\n')
        f.write('ddf-pipeline version was '+__version__+'\n')
        from DDFacet.DDF import report_version as ddf_version
        f.write('DDF version was '+ddf_version()+'\n')
        from killMS.Other.logo import report_version as killms_version
        f.write('killMS version was '+killms_version()+'\n')
        if o['do_dynspec']:
            from DynSpecMS import dynspecms_version
            f.write('DynSpecMS version was '+dynspecms_version.version()+'\n\n')
        f.write('Options dictionary was as follows:\n')
        for k in o:
            f.write("%-20s : %s\n" % (k,str(o[k])))

def stop(v=2):
    if use_database():
        update_status(None,'Stopped')
    sys.exit(v)
            
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
    error=False
    report('Checking for IMAGING_WEIGHT in input MSS')
    mslist=[s.strip() for s in open(mslist_name).readlines()]
    for ms in mslist:
        try:
            t = pt.table(ms)
        except RuntimeError:
            print('Failed to open table',ms,'-- table may be missing or corrupt')
            error=True
        else:
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
    if error:
        raise RuntimeError('One or more tables failed to open')
    return result

def parse_parset(parsets,use_headings=False):
    keywords={}
    for parset in parsets:
        if os.path.isfile(parset):
            break
    else:
        parset=None
    if parset is not None:
        with open(parset) as infile:
            lines=infile.readlines()
        prefix=''
        for l in lines:
            bits=l.split()
            if use_headings and l[0]=='[':
                prefix=bits[0][1:-1]+'-'
            if len(bits)>0 and l[0]!='#' and l[0]!='_' and not(l[0].isspace()) and l[0]!='[':
                if len(bits)>2:
                    content=bits[2]
                    if content[0]=='#':
                        content=''
                else:
                    content=None
                keywords[prefix+bits[0]]=content

        return keywords
    
    else:
        warn('Cannot find parset, some features may not work')
        return {}

def ddf_shift(imagename,shiftfile,catcher=None,options=None,dicomodel=None,verbose=False):
    if catcher: catcher.check()
    if options is None:
        options=o # attempt to get global if it exists

    keywords=parse_parset([os.environ['DDF_DIR']+'/DDFacet/DDFacet/Parset/DefaultParset.cfg'],use_headings=True)
        
    cache_dir=find_cache_dir(options)
    if dicomodel is None:
        dicomodel=imagename+'.DicoModel'
    runcommand='DDF.py '+imagename+'.parset --Misc-ConserveMemory=1  --Output-Name='+imagename+'_shift --Output-Mode=RestoreAndShift --Output-ShiftFacetsFile='+shiftfile+' --Predict-InitDicoModel '+dicomodel+' --Cache-SmoothBeam=force --Log-Memory 1 --Cache-Dir='+cache_dir
    if 'Misc-IgnoreDeprecationMarking' in keywords:
        runcommand+=' --Misc-IgnoreDeprecationMarking=1'
    
    fname=imagename+'_shift.app.facetRestored.fits'
    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF-shift step')
        if verbose:
            print('would have run',runcommand)
    else:
         run(runcommand,dryrun=options['dryrun'],log=logfilename('DDF-'+imagename+'_shift.log',options=options),quiet=options['quiet'])


def ddf_image(imagename,mslist,cleanmask=None,cleanmode='HMP',ddsols=None,applysols=None,threshold=None,majorcycles=3,use_dicomodel=False,robust=0,beamsize=None,beamsize_minor=None,beamsize_pa=None,reuse_psf=False,reuse_dirty=False,verbose=False,saveimages=None,imsize=None,cellsize=None,uvrange=None,colname='CORRECTED_DATA',peakfactor=0.1,dicomodel_base=None,options=None,do_decorr=None,normalization=None,dirty_from_resid=False,clusterfile=None,HMPsize=None,automask=True,automask_threshold=10.0,smooth=False,noweights=False,cubemode=False,apply_weights=True,use_weightspectrum=False,catcher=None,rms_factor=3.0,predict_column=None,conditional_clearcache=False,PredictSettings=None,RMSFactorInitHMP=1.,MaxMinorIterInitHMP=10000,OuterSpaceTh=None,AllowNegativeInitHMP=False,phasecenter=None,polcubemode=False,channels=None,startchan=None,endchan=None,stokes=None):

    if catcher: catcher.check()

    # saveimages lists _additional_ images to save
    if saveimages is None:
        saveimages=''
    saveimages+='onNeds'
    if options is None:
        options=o # attempt to get global if it exists

    keywords=parse_parset([os.environ['DDF_DIR']+'/DDFacet/DDFacet/Parset/DefaultParset.cfg'],use_headings=True)
        
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

    runcommand = "DDF.py --Misc-ConserveMemory=1 --Output-Name=%s --Data-MS=%s --Deconv-PeakFactor %f --Data-ColName %s --Parallel-NCPU=%i --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=%s --Deconv-Mode %s --Beam-Model=LOFAR --Beam-PhasedArrayMode=A --Weight-Robust %f --Image-NPix=%i --CF-wmax 50000 --CF-Nw 100 --Output-Also %s --Image-Cell %f --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=%f --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=%s --Cache-DirWisdomFFTW=%s --Debug-Pdb=never --Log-Memory 1"%(imagename,mslist,peakfactor,colname,options['NCPU_DDF'],majorcycles,cleanmode,robust,imsize,saveimages,float(cellsize),rms_factor,cache_dir,cache_dir)

    runcommand += " --GAClean-RMSFactorInitHMP %f"%RMSFactorInitHMP
    runcommand += " --GAClean-MaxMinorIterInitHMP %f"%MaxMinorIterInitHMP
    if AllowNegativeInitHMP:
        runcommand += " --GAClean-AllowNegativeInitHMP True"
    if OuterSpaceTh is not None:
        runcommand += " --HMP-OuterSpaceTh %f"%OuterSpaceTh
    if options['use_splitisland']:
       runcommand += " --SSDClean-MaxIslandSize 100"

    runcommand+=' --DDESolutions-SolsDir=%s'%options["SolsDir"]
    runcommand+=' --Cache-Weight=reset'

    if 'Misc-IgnoreDeprecationMarking' in keywords:
        runcommand+=' --Misc-IgnoreDeprecationMarking=1'

    if 'Beam-At' in keywords:
        runcommand+=' --Beam-At=%s'%options['beam_at']
        
    if PredictSettings is None:
        runcommand += " --Output-Mode=Clean"
    else:
        if len(PredictSettings) == 2:
            runcommand += " --Output-Mode=%s --Predict-ColName %s"%PredictSettings
        elif len(PredictSettings) == 3:
            runcommand += " --Output-Mode=%s --Predict-ColName %s --Predict-MaskSquare [0,%i]"%PredictSettings
        else:
            raise RuntimeError('PredictSettings has the wrong dimensions %s '%PredictSettings)

    if beamsize_minor is not None:
        runcommand += ' --Output-RestoringBeam %f,%f,%f'%(beamsize,beamsize_minor,beamsize_pa)
    elif beamsize is not None:
        runcommand += ' --Output-RestoringBeam %f'%(beamsize)
    
    if apply_weights:
        runcommand+=' --Weight-ColName="IMAGING_WEIGHT"'
    else:
        if not use_weightspectrum:
            runcommand+=' --Weight-ColName="None"'
        else:
            runcommand+=' --Weight-ColName="WEIGHT_SPECTRUM"'

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

    if polcubemode:
        runcommand+=' --Output-Cubes=dD --RIME-PolMode=QU --Output-Mode=Dirty  --Freq-NBand=%i --Selection-ChanStart=%s --Selection-ChanEnd=%s' % (channels,startchan,endchan)

    if not cubemode and not polcubemode:
        runcommand+=' --Freq-NBand=2'
    if stokes:
        runcommand +=' --RIME-PolMode=%s --Output-Mode=Dirty'%stokes


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
        
    if phasecenter is not None:
        runcommand += " --Image-PhaseCenterRADEC=[%s,%s]"%(phasecenter[0],phasecenter[1])
    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF step')
        if verbose:
            print('would have run',runcommand)
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
        hdus.writeto(fname,overwrite=True)
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
    filenames=[Name+'.app.restored.pybdsm.srl.fits',Name+'.app.restored.pybdsf.srl.fits']
    for filename in filenames:
        if os.path.isfile(filename):
            break
    else:
        die('Catalogue file does not exist!')
    if use_makemask_products:
        runcommand="ClusterCat.py --SourceCat %s --AvoidPolygons MaskDiffuse.pickle --DoPlot=0 --NGen 100 --NCPU %i"%(filename,options['NCPU_DDF'])
    else:
        runcommand="ClusterCat.py --SourceCat %s --DoPlot=0 --NGen 100 --NCPU %i"%(filename,options['NCPU_DDF'])
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
            print('Would have run',runcommand)
    else:
        run(runcommand,dryrun=options['dryrun'],log=logfilename('MM-'+imagename+'.log',options=options),quiet=options['quiet'])
        if external_mask is not None:
            if isinstance(external_mask,list) or isinstance(external_mask,tuple):
                for mask in external_mask:
                    merge_mask(fname,mask,fname)
            else:
                merge_mask(fname,external_mask,fname)
    return fname
            


def killms_data(imagename,mslist,outsols,clusterfile=None,colname='CORRECTED_DATA',niterkf=6,dicomodel=None,
                uvrange=None,wtuv=None,robust=None,catcher=None,dt=None,options=None,
                SolverType="KAFCA",PolMode="Scalar",MergeSmooth=False,NChanSols=1,
                DISettings=None,EvolutionSolFile=None,CovQ=0.1,InterpToMSListFreqs=None,
                SkipSmooth=False,PreApplySols=None,SigmaFilterOutliers=None):

    if options is None:
        options=o # attempt to get global if it exists

    cache_dir=find_cache_dir(options)

    if 'KILLMS_DIR' in os.environ:
        # different versions have different locations for the parset, so check them all
        keywords=parse_parset([os.environ['KILLMS_DIR']+'/killMS/killMS/Parset/DefaultParset.cfg',os.environ['KILLMS_DIR']+'/killMS/Parset/DefaultParset.cfg'])
    else:
        keywords={}
    
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
            runcommand = "kMS.py --MSName %s --SolverType %s --PolMode %s --BaseImageName %s --dt %f --NIterKF %i --CovQ %f --LambdaKF=%f --NCPU %i --OutSolsName %s --InCol %s"%(f,SolverType,PolMode,imagename,dt,niterkf, CovQ, o['LambdaKF'], o['NCPU_killms'], outsols,colname)

            # check for option to stop pdb call and use it if present
            
            if 'DebugPdb' in keywords:
                runcommand+=' --DebugPdb=0'
                
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
            
            if PreApplySols:
                runcommand+=' --PreApplySols=[%s]'%PreApplySols

                
            if DISettings is None:
                runcommand+=' --NChanSols %i' % NChanSols
                runcommand+=' --BeamMode LOFAR --PhasedArrayMode=A --DDFCacheDir=%s'%cache_dir
                if 'BeamAt' in keywords:
                    runcommand+=' --BeamAt=%s'%options['beam_at']

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
        outsols=smooth_solutions(mslist,outsols,catcher=None,dryrun=o['dryrun'],InterpToMSListFreqs=InterpToMSListFreqs,
                                 SkipSmooth=SkipSmooth,SigmaFilterOutliers=SigmaFilterOutliers)
        



    return outsols

def compress_fits(filename,q):
    command='fpack -q %i %s' % (q,filename)
    run(command,dryrun=o['dryrun'])
    
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
        print('Removing',f)
        rmtree(f)

def _basename(path):
    return os.path.basename(path.rstrip(os.path.sep))

def mvglob(path,dest):
    g=glob.glob(path)
    for f in g:
        print('Moving',f,'to',dest)
        # work round shutil non-overwriting behaviour
        real_dst = os.path.join(dest, _basename(f))
        print('Target is',real_dst)
        if os.path.exists(real_dst):
            if os.path.isdir(real_dst):
                rmtree(real_dst)
            else:
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

def smooth_solutions(mslist,ddsols,catcher=None,dryrun=False,InterpToMSListFreqs=None,SkipSmooth=False,SigmaFilterOutliers=None):
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
        with open('solslist_%.2f.txt'%start_time,'w') as f:
            for i in range(0,len(full_sollist)):
                if start_times[i] == start_time:
                    solname = full_sollist[i]
                    f.write('%s\n'%(solname))
        
        checkname='%s_%.2f_merged.npz'%(ddsols,start_time)
        if o['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running MergeSols step')
        else:
            ss='MergeSols.py --SolsFilesIn=solslist_%.2f.txt --SolFileOut=%s '%(start_time,checkname)
            if SigmaFilterOutliers:
                ss+=" --SigmaFilterOutliers %f"%SigmaFilterOutliers
            run(ss,dryrun=dryrun)
            
        checkname='%s_%.2f_smoothed.npz'%(ddsols,start_time)
        if o['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running SmoothSols step')
        elif SkipSmooth:
            warn('Skipping smoothing Solutions file')
        else:
            run('SmoothSols.py --SolsFileIn=%s_%.2f_merged.npz --SolsFileOut=%s --InterpMode=%s --NCPU=%s'%(ddsols,start_time,checkname,o['smoothingtype'],o['NCPU_killms']),dryrun=dryrun)

        smoothoutname='%s_%.2f_smoothed.npz'%(ddsols,start_time)

        if InterpToMSListFreqs:
            interp_outname="%s_%.2f_interp.npz"%(smoothoutname,start_time)
            checkname=interp_outname
            if o['restart'] and os.path.isfile(checkname):
                warn('Solutions file '+checkname+' already exists, not running InterpSols step')
            else:
                command="InterpSols.py --SolsFileIn %s --SolsFileOut %s --MSOutFreq %s --NCPU=%s"%(smoothoutname,interp_outname,InterpToMSListFreqs,o['NCPU_killms'])
                run(command,dryrun=dryrun)
        
        for i in range(0,len(full_sollist)):
            if start_times[i] == start_time:
                if not SkipSmooth:
                    symsolname = full_sollist[i].replace(ddsols,ddsols+'_smoothed')
                else:
                    symsolname = full_sollist[i].replace(ddsols,ddsols+'_merged')                 
                # always overwrite the symlink to allow the dataset to move -- costs nothing
                if os.path.islink(symsolname):
                    warn('Symlink ' + symsolname + ' already exists, recreating')
                    os.unlink(symsolname)

                if not SkipSmooth:
                    os.symlink(os.path.abspath('%s_%.2f_smoothed.npz'%(ddsols,start_time)),symsolname)
                else:
                    os.symlink(os.path.abspath('%s_%.2f_merged.npz'%(ddsols,start_time)),symsolname)
                    
                    
        if SkipSmooth:
            outname = ddsols + '_merged'
        else:
            outname = ddsols + '_smoothed'

    return outname

def full_clearcache(o,extras=None):
    clearcache(o['mslist'],o)
    clearcache('temp_mslist.txt',o)
    if o['full_mslist'] is not None:
        clearcache(o['full_mslist'],o)
    if extras is not None:
        for mslist in extras:
            clearcache(mslist,o)

def subtract_data(mslist,col1,col2):
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    for f in filenames:
        print('Subtracting',f)
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
    warn('Using (dt,df)=(%i,%i) for self-cal run of %s with (<|model|>,std)=(%.2f,%.2f) giving SNR=%.2f'%(nt_step,nch_step,msname,M,S,SNR))
    
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

def ingest_dynspec(obsid='*'):
    report('Ingesting dynamic spectra (%s) into the database' % obsid)
    with SurveysDB() as sdb:
        sdb.cur.execute('lock table spectra write')
        field=os.path.basename(os.getcwd())
        g=glob.glob('DynSpecs_'+obsid)
        for f in g:
            if '.tgz' in f:
                continue
            bits=f.split('_')
            obsid=bits[1]
            CatName=f+'/Catalog.npy'
            print("Loading %s"%CatName)
            try:
                catalogue=np.load(CatName)
            except:
                print("   %s does not exist"%CatName)
                continue
            # match filenames to names
            fd={}
            for r in catalogue:
                name=r['Name']
                if isinstance(name,np.bytes_):
                    name=name.decode('utf-8')
                fd[name]=''
            gf=glob.glob(f+'/TARGET/*.fits')+glob.glob(f+'/OFF/*.fits')
            for ff in gf:
                hdu=fits.open(ff)
                name=hdu[0].header['NAME']
                assert(name in fd)
                fd[name]=ff
                hdu.close()
            sdb.cur.execute('delete from spectra where obsid="%s"' % obsid)
            for i,r in enumerate(catalogue):
                name=r['Name']
                if isinstance(name,np.bytes_):
                    name=name.decode('utf-8')
                sExec='insert into spectra values ( "%s", "%s", "%s", "%s", "%s", "%s", %.7f, %.7f, %g, %g, %g, %g )' % (field+'_'+obsid+'_'+str(i), name, r['Type'], field, obsid, fd[name], r['ra']*180.0/np.pi, r['dec']*180.0/np.pi, r['FluxI'], r['FluxV'], r['sigFluxI'], r['sigFluxV'])
                print(sExec)
                sdb.cur.execute(sExec)
        
    

def subtract_vis(mslist=None,colname_a="CORRECTED_DATA",colname_b="DATA_SUB",out_colname="DATA_SUB"):
    from pyrap.tables import table
    f=open(mslist)
    mslist=f.readlines()
    mslist=[msname.replace("\n","") for msname in mslist]
    for msname in mslist:
        report('Subtracting: %s = %s - %s'%(out_colname,colname_a,colname_b))
        t=table(msname,readonly=False)
        d=t.getcol(colname_a)
        p=t.getcol(colname_b)
        d-=p
        if out_colname not in t.colnames():
            report('Adding column %s in %s'%(out_colname,msname))
            desc=t.getcoldesc(colname_a)
            desc["name"]=out_colname
            desc['comment']=desc['comment'].replace(" ","_")
            t.addcols(desc)
        t.putcol(out_colname,d)
        t.close()
    

def subtractOuterSquare(o):
    
    wide_imsize=o['wide_imsize']
    NPixSmall=o['imsize'] #int(NPixLarge/float(o['fact_reduce_field']))
    colname=o['colname']


    wide_uvrange=[o['image_uvmin'],2.5*206.0/o['wide_psf_arcsec']]
    
    killms_uvrange=[0,1000]
    if o['solutions_uvmin'] is not None:
        killms_uvrange[0]=o['solutions_uvmin']


    if o['catch_signal']:
        catcher=Catcher()
    else:
        catcher=None
    
    
    #if o['wide_psf_arcsec'] is not None:
    # wide-res image requested
    #if o['wide_imsize'] is not None:
        #wide_imsize=o['wide_imsize'] # allow over-ride
    #else:
        #wide_imsize=o['imsize']*o['cellsize']/o['wide_cell']
    extmask=None

    ddf_image('image_full_wide',o['mslist'],
                cleanmask=extmask,
                cleanmode='SSD',
                AllowNegativeInitHMP=True,
                majorcycles=2,robust=o['wide_robust'],
                colname=colname,use_dicomodel=False,
                uvrange=wide_uvrange,beamsize=o['wide_psf_arcsec'],
                imsize=o['wide_imsize'],cellsize=o['wide_cell'],peakfactor=0.001,
                apply_weights=False,use_weightspectrum=o['use_weightspectrum'],
                smooth=True,automask=True,automask_threshold=o['thresholds'][0],normalization=o['normalize'][2],
                catcher=catcher)


    external_mask='wide_external_mask.fits'
    make_external_mask(external_mask,'image_full_wide.dirty.fits',use_tgss=True,clobber=False)
    
    make_mask('image_full_wide.app.restored.fits',3.0,external_mask=external_mask,catcher=catcher)
    
    ddf_image('image_full_wide_im',o['mslist'],
            cleanmask='image_full_wide.app.restored.fits.mask.fits',
            cleanmode='SSD',
            AllowNegativeInitHMP=True,
            majorcycles=1,robust=o['wide_robust'],
            uvrange=wide_uvrange,beamsize=o['wide_psf_arcsec'],
            imsize=o['wide_imsize'],cellsize=o['wide_cell'],peakfactor=0.001,
            apply_weights=False,use_weightspectrum=o['use_weightspectrum'],
            smooth=True,automask=True,automask_threshold=o['thresholds'][0],normalization=o['normalize'][2],colname=colname,
            reuse_psf=True,dirty_from_resid=True,use_dicomodel=True,dicomodel_base='image_full_wide',
            catcher=catcher)


    # predict outside the central rectangle
    
    NpixMaskSquare = np.floor(0.95*o['imsize']*o['cellsize']/o['wide_cell'])
    
    FileHasPredicted='image_full_wide_predict.HasPredicted'
    if o['restart'] and os.path.isfile(FileHasPredicted):
        warn('File %s already exists, skipping Predict step'%FileHasPredicted)
    else:
        ddf_image('image_full_wide_predict',o['full_mslist'],colname=colname,robust=o['wide_robust'],
            cleanmask='image_full_wide.app.restored.fits.mask.fits',
                  cleanmode='SSD',
                  #majorcycles=1,automask=True,automask_threshold=o['thresholds'][1],
                  #ddsols='wide_killms_p1',
                  #applysols='AP',#normalization=o['normalize'][0],
                  peakfactor=0.001,
                  apply_weights=False,use_weightspectrum=o['use_weightspectrum'],
                  #uvrange=wide_uvrange,beamsize=o['wide_psf_arcsec'],
                  beamsize=o['wide_psf_arcsec'],
                  imsize=o['wide_imsize'],cellsize=o['wide_cell'],
                  use_dicomodel=True,catcher=catcher,
                  PredictSettings=("Predict","DATA_SUB",NpixMaskSquare),
                  dicomodel_base='image_full_wide_im')
        os.system("touch %s"%FileHasPredicted)



    # subtract predicted visibilities
    FileHasSubtracted='image_full_wide_predict.HasSubtracted'
    if o['restart'] and os.path.isfile(FileHasSubtracted):
        warn('File %s already exists, skipping subtract vis step'%FileHasSubtracted)
    else:
        subtract_vis(mslist=o['full_mslist'],colname_a=colname,colname_b="DATA_SUB",out_colname="DATA_SUB")
        os.system("touch %s"%FileHasSubtracted)


    ## test subtracted...
    ## sanity check
    ddf_image('image_full_wide_im_sub',o['mslist'],
            cleanmask='image_full_wide.app.restored.fits.mask.fits',
            cleanmode='SSD',
            AllowNegativeInitHMP=True,
            majorcycles=1,robust=o['wide_robust'],
            uvrange=wide_uvrange,beamsize=o['wide_psf_arcsec'],
            imsize=o['wide_imsize'],cellsize=o['wide_cell'],peakfactor=0.001,
            apply_weights=False,use_weightspectrum=o['use_weightspectrum'],
            smooth=True,automask=True,automask_threshold=o['thresholds'][0],normalization=o['normalize'][2],colname='DATA_SUB',
            reuse_psf=True,dirty_from_resid=False,use_dicomodel=False,
            catcher=catcher)


def main(o=None):
    if o is None and MyPickle is not None:
        o=MyPickle.Load("ddf-pipeline.last")

    if '$$' in o['tgss'] or np.any(['$$' in l for l in o['catalogues']]):
        if "DDF_PIPELINE_CATALOGS" not in list(os.environ.keys()):
            die("You need to define the environment variable DDF_PIPELINE_CATALOGS where your catalogs are located")
        o["tgss"]=o["tgss"].replace("$$",os.environ["DDF_PIPELINE_CATALOGS"])
        o["catalogues"]=[l.replace("$$",os.environ["DDF_PIPELINE_CATALOGS"]) for l in o["catalogues"]]
        
    lCat=o["catalogues"]+[o["tgss"]]
    for fCat in lCat:
        if not os.path.isfile(fCat):
            die("Catalog %s does not exist"%fCat)

    if o['catch_signal']:
        catcher=Catcher()
    else:
        catcher=None

    if o['remove_columns']:
        warn('Removing all pipeline-created columns')
        run('remove_columns.py '+o['full_mslist'],log=None,dryrun=o['dryrun'])
        
    uvrange=[o['image_uvmin'],o['uvmax']]
    killms_uvrange=[0,1000]
    if o['solutions_uvmin'] is not None:
        killms_uvrange[0]=o['solutions_uvmin']
    if o['mslist'] is None:
        die('MS list must be specified')

    # Set column name for first steps
    colname=o['colname']

    # Check if the column exists in one MS. Important to do this
    # before we check imaging weights, because that will create empty
    # versions of e.g. CORRECTED_DATA
    mslist=[s.strip() for s in open(o['mslist']).readlines()]
    t = pt.table(mslist[0])
    try:
        dummy=t.getcoldesc(colname)
    except RuntimeError:
        dummy=None
    t.close()
    if dummy is None:
        die('Dataset does not contain the column "%s"' % colname)
    
    # Clear the shared memory
    run('CleanSHM.py',dryrun=o['dryrun'])    

    # Pipeline started!
    if use_database():
        update_status(None,'Running',time='start_date')
    
    if o['redofrom']:

        if not os.path.isdir(o['archive_dir']):
            os.mkdir(o['archive_dir'])

        # Redofrom as a concept no longer really works because of the
        # re-use of columns in the DI steps.  Hence the only options
        # here are 'start', which removes all but the MS and mslist
        # files, and 'dirin' which retains the first
        # direction-independent images. Both of these will strip out
        # all of the extra columns in the MS and (by removing SOLSDIR)
        # remove all of the solutions.

        report('Removing old files for a redo from '+o['redofrom'])
        files=glob.glob('*')
        keep=glob.glob('*.ms')+[o['mslist'],o['full_mslist'],o['archive_dir']]+glob.glob('*.cfg')
        if o['clusterfile'] is not None:
            keep.append(o['clusterfile'])
        if o['redofrom']=='start':
            pass
        elif o['redofrom']=='dirin':
            keep+=glob.glob('image_dirin_SSD_init.*') + glob.glob('image_dirin_SSD.*') + glob.glob('image_dirin_SSD_m.*') + glob.glob('MaskDiffuse*') + glob.glob('Noise*.fits')
        else:
            die('Redofrom option not implemented')
            
        if o['full_mslist'] is not None:
            run('remove_columns.py '+o['full_mslist'],log=None,dryrun=o['dryrun'])
        else:
            run('remove_columns.py '+o['mslist'],log=None,dryrun=o['dryrun'])
        for f in files:
            if f not in keep:
                mvglob(f,o['archive_dir'])

        if o['exitafter'] == 'cleanup':
            warn('User specified exit after cleanup')
            stop(2)

                
    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])
       
    # Check imaging weights -- needed before DDF
    new=check_imaging_weight(o['mslist'])
    if o['clearcache'] or new or o['redofrom']:
        # Clear the cache, we don't know where it's been. If this is a
        # completely new dataset it is always safe (and required) to
        # clear the cache -- solves problems where the cache is not
        # stored per dataset. If we are redoing, cache needs to be removed
        full_clearcache(o)

    # ##########################################################
    if o['redo_DI']:
        separator('Redo DI correction')
        redo_dppp_di(o)

    # ##########################################################
    # subtract outer square
    if o['do_wide']:
        subtractOuterSquare(o)
        colname="DATA_SUB"
        #ReduceFactor=o['fact_reduce_field']
        #NPixSmall=int(o['imsize']/float(ReduceFactor))
        #o['imsize']=NPixSmall
        #o['ndir']=int(o['ndir']/float(ReduceFactor))

    # start of 'Big If' for reducing multiple datasets with a pre-made sky model
    if o['basedicomodel'] is None:
        # ##########################################################
        # Initial dirty image to allow an external (TGSS) mask to be made
        separator("Initial dirty")
        ddf_image('image_dirin_SSD_init',o['mslist'],cleanmask=None,cleanmode='SSD',majorcycles=0,robust=o['image_robust'],
                  reuse_psf=False,reuse_dirty=False,peakfactor=0.05,colname=colname,clusterfile=None,
                  apply_weights=o['apply_weights'][0], use_weightspectrum=o['use_weightspectrum'], uvrange=uvrange,catcher=catcher)

        separator("External mask")
        external_mask='external_mask.fits'
        make_external_mask(external_mask,'image_dirin_SSD_init.dirty.fits',use_tgss=True,clobber=False)
        
        if o['external_fits_mask'] is not None:
            merge_mask(external_mask,o['external_fits_mask'],external_mask)

        # Deep SSD clean with this external mask and automasking
        separator("DI Deconv (externally defined sources)")
        CurrentBaseDicoModelName=ddf_image('image_dirin_SSD',o['mslist'],cleanmask=external_mask,cleanmode='SSD',
                                           majorcycles=1,robust=o['image_robust'],reuse_psf=True,reuse_dirty=True,
                                           peakfactor=0.01,rms_factor=3,
                                           colname=colname,clusterfile=None,automask=True,
                                           automask_threshold=o['thresholds'][0],apply_weights=o['apply_weights'][0], use_weightspectrum=o['use_weightspectrum'],
                                           uvrange=uvrange,catcher=catcher)
    
        separator("Make the diffuse emission mask")
        # Make the diffuse emission mask
        _=make_mask('image_dirin_SSD.residual01.fits',
                    o['thres_outmaskextended'],
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
                                           automask_threshold=o['thresholds'][0],apply_weights=o['apply_weights'][0],use_weightspectrum=o['use_weightspectrum'],
                                           uvrange=uvrange,catcher=catcher,
                                           RMSFactorInitHMP=1.,
                                           MaxMinorIterInitHMP=10000,
                                           PredictSettings=("Clean","DD_PREDICT"))

        if o['exitafter'] == 'initial':
            warn('User specified exit after initial image')
            stop(2)


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
        if o['clearcache'] or new or o['redofrom']:
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
                                           apply_weights=o['apply_weights'][0],use_weightspectrum=o['use_weightspectrum'],
                                           uvrange=uvrange,catcher=catcher,
                                           RMSFactorInitHMP=1.,
                                           MaxMinorIterInitHMP=10000,
                                           PredictSettings=("Clean","DD_PREDICT"))

        if o['exitafter'] == 'dirin':
            warn('User specified exit after image_dirin.')
            stop(2)

        if not o['skip_di']:
            separator("DI CAL")
            ########################
            killms_data('PredictDI_0',o['mslist'],'DIS0',colname=colname,
                        dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                        niterkf=o['NIterKF'][0],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                        catcher=catcher,
                        dt=o['dt_di'],
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

            # make a mask from the full-res image
            separator("Make mask for next iteration")
            CurrentMaskName=make_mask('image_dirin_SSD_m_c_di_m.app.restored.fits',
                                    o['thresholds'][1],
                                    external_mask=external_mask,
                                    catcher=catcher)

            if o['exitafter'] == 'dirin_di':
                warn('User specified exit after image_dirin with DI calibration.')
                stop(2)


        separator("DD calibration")
        CurrentDDkMSSolName=killms_data(CurrentBaseDicoModelName,o['mslist'],'DDS0',colname=colname,
                                        dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                        clusterfile=ClusterFile,
                                        CovQ=0.02,
                                        niterkf=o['NIterKF'][1],
                                        #CovQ=0.1,
                                        #niterkf=6,
                                        uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],dt=o['dt_slow'],
                                        catcher=catcher,NChanSols=o['NChanSols'],
                                        MergeSmooth=o['smoothing'])

        # ##########################################################
        # run bootstrap, and change the column name if it runs
        if o['bootstrap']:
            separator("Bootstrap")
            report('Running bootstrap')
            run('bootstrap.py '+' '.join(sys.argv[1:]),log=None,dryrun=o["dryrun"])
            colname=colname+'_SCALED' # DI corrected, scaled
            if o['exitafter'] == 'bootstrap':
                warn('User specified exit after bootstrap.')
                stop(2)


        separator("PhaseOnly deconv")
        print('Smoothing is',o['smoothing'],'Current DDkMS name is',CurrentDDkMSSolName)
        CurrentBaseDicoModelName=ddf_image('image_phase1',o['mslist'],
                                           cleanmask=CurrentMaskName,
                                           cleanmode='SSD',
                                           ddsols=CurrentDDkMSSolName,applysols=o['apply_sols'][0],majorcycles=2,robust=o['image_robust'],
                                           colname=colname,peakfactor=0.001,automask=True,
                                           automask_threshold=o['thresholds'][1],
                                           normalization=o['normalize'][0],apply_weights=o['apply_weights'][1],use_weightspectrum=o['use_weightspectrum'],uvrange=uvrange,
                                           use_dicomodel=True,
                                           dicomodel_base=CurrentBaseDicoModelName,
                                           catcher=catcher,
                                           RMSFactorInitHMP=1.,
                                           MaxMinorIterInitHMP=10000,
                                           PredictSettings=("Clean","DD_PREDICT"))

        if o['exitafter'] == 'phase':
            warn('User specified exit after phase-only deconvolution.')
            stop(2)

        separator("Mask for deeper deconv")
        CurrentMaskName=make_mask('image_phase1.app.restored.fits',o['thresholds'][1],external_mask=external_mask,catcher=catcher)
        CurrentBaseDicoModelName=mask_dicomodel('image_phase1.DicoModel',CurrentMaskName,'image_phase1_masked.DicoModel',catcher=catcher)

        separator("DD calibration")
        CurrentDDkMSSolName=killms_data('image_phase1',o['mslist'],'DDS1',colname=colname,
                                        dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                        CovQ=0.02,
                                        clusterfile=ClusterFile,
                                        niterkf=o['NIterKF'][2],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                                        dt=o['dt_slow'],
                                        catcher=catcher,NChanSols=o['NChanSols'],
                                        EvolutionSolFile=CurrentDDkMSSolName,
                                        MergeSmooth=o['smoothing'])
        ##############################################

        separator("AmpPhase deconv")
        CurrentBaseDicoModelName=ddf_image('image_ampphase1',o['mslist'],
                                       cleanmask=CurrentMaskName,cleanmode='SSD',
                                       ddsols=CurrentDDkMSSolName,applysols=o['apply_sols'][1],majorcycles=1,robust=o['image_robust'],
                                       colname=colname,peakfactor=0.001,automask=True,
                                       automask_threshold=o['thresholds'][1],
                                       normalization=o['normalize'][0],apply_weights=o['apply_weights'][1],use_weightspectrum=o['use_weightspectrum'],uvrange=uvrange,
                                       use_dicomodel=True,
                                       dicomodel_base=CurrentBaseDicoModelName,
                                       catcher=catcher,
                                       RMSFactorInitHMP=1.,
                                       #AllowNegativeInitHMP=True,
                                       MaxMinorIterInitHMP=10000,
                                       PredictSettings=("Clean","DD_PREDICT"))

        if o['exitafter'] == 'ampphase':
            warn('User specified exit after amp-phase deconvolution.')
            stop(2)

        separator("Make Mask")
        CurrentMaskName=make_mask('image_ampphase1.app.restored.fits',o['thresholds'][1],external_mask=external_mask,catcher=catcher)
        CurrentBaseDicoModelName=mask_dicomodel('image_ampphase1.DicoModel',CurrentMaskName,'image_ampphase1m_masked.DicoModel',catcher=catcher)

        if not o['skip_di']:
            separator("Second DI calibration")
            ddf_image('Predict_DI1',o['mslist'],
                    cleanmask=CurrentMaskName,cleanmode='SSD',
                    ddsols=CurrentDDkMSSolName,applysols=o['apply_sols'][2],majorcycles=1,robust=o['image_robust'],
                    colname=colname,peakfactor=0.001,automask=True,
                    automask_threshold=o['thresholds'][1],
                    normalization=o['normalize'][0],apply_weights=o['apply_weights'][1],use_weightspectrum=o['use_weightspectrum'],uvrange=uvrange,
                    use_dicomodel=True,
                    dicomodel_base=CurrentBaseDicoModelName,
                    catcher=catcher,
                    RMSFactorInitHMP=1.,
                    #AllowNegativeInitHMP=True,
                    MaxMinorIterInitHMP=10000,
                    PredictSettings=("Predict","DD_PREDICT"))

            
            separator("Another DI step")
            if o['bootstrap']:
                colname='SCALED_DATA'
            else:
                if o['do_wide']:
                    colname='DATA_SUB'
                else:
                    colname=o['colname']
            killms_data('PredictDI_1',o['mslist'],'DIS1',colname=colname,
                        dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                        #clusterfile=ClusterFile,
                        niterkf=o['NIterKF'][3],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                        catcher=catcher,
                        dt=o['dt_di'],
                        DISettings=("CohJones","IFull","DD_PREDICT","DATA_DI_CORRECTED"))
            # cubical_data(o['mslist'],
            #              NameSol="DIS1",
            #              n_dt=1,
            #              n_df=2,
            #              n_DT=None,
            #              DataColName=o['colname'],
            #              ModelColName="DD_PREDICT",
            #              OutColName="DATA_DI_CORRECTED")

            colname='DATA_DI_CORRECTED' # again
            CurrentBaseDicoModelName=ddf_image('image_ampphase1_di',o['mslist'],
                                        cleanmask=CurrentMaskName,cleanmode='SSD',
                                        ddsols=CurrentDDkMSSolName,applysols=o['apply_sols'][3],
                                        majorcycles=1,robust=o['image_robust'],
                                        colname=colname,peakfactor=0.001,automask=True,
                                        automask_threshold=o['thresholds'][1],
                                        normalization=o['normalize'][0],
                                        apply_weights=o['apply_weights'][1],use_weightspectrum=o['use_weightspectrum'],uvrange=uvrange,
                                        use_dicomodel=True,
                                        dicomodel_base=CurrentBaseDicoModelName,
                                        catcher=catcher,
                                        RMSFactorInitHMP=1.,
                                        #AllowNegativeInitHMP=True,
                                        MaxMinorIterInitHMP=10000,
                                        PredictSettings=("Clean","DD_PREDICT"))

            if o['exitafter'] == 'ampphase_di':
                warn('User specified exit after amp-phase plus DI deconvolution.')
                stop(2)


        if o['full_mslist'] is None:
            warn('No full mslist provided, stopping here')
            summary(o)
            stop(3)
        
        # #########################################################################
        # ###############                  BIG MSLIST               ###############
        # #########################################################################

        # check full mslist imaging weights
        check_imaging_weight(o['full_mslist'])

        if o['bootstrap']:
            colname='SCALED_DATA'
        else:
            if o['do_wide']:
                colname='DATA_SUB'
            else:
                colname=o['colname']

        if not o['skip_di']:
            separator("Make Mask")
            CurrentMaskName=make_mask('image_ampphase1_di.app.restored.fits',o['thresholds'][1],external_mask=external_mask,catcher=catcher)
            CurrentBaseDicoModelName=mask_dicomodel('image_ampphase1_di.DicoModel',CurrentMaskName,'image_ampphase1_di_masked.DicoModel',catcher=catcher)
            CurrentImageName= 'image_ampphase1_di'
        else:
            CurrentImageName = 'image_ampphase1'

    else:
        # alternative branch of massive if!
        if o['clusterfile'] is None:
            warn('No clusterfile provided, stopping here')
            summary(o)
            stop(4)
        if o['baseimagename'] is None:
            warn('No baseimage provided, stopping here')
            summary(o)
            stop(4)
        if o['basemaskname'] is None:
            warn('No mask file provided, stopping here')
            summary(o)
            stop(4)
        ClusterFile=o['clusterfile']
        CurrentMaskName = o['basemaskname']
        CurrentBaseDicoModelName = o['basedicomodel']
        CurrentImageName = o['baseimagename']
        external_mask = CurrentMaskName

    separator("DD calibration of full mslist")
    CurrentDDkMSSolName=killms_data(CurrentImageName,o['full_mslist'],'DDS2_full',
                                    colname=colname,
                                    dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                    CovQ=0.1,
                                    clusterfile=ClusterFile,
                                    niterkf=o['NIterKF'][4],
                                    uvrange=killms_uvrange,
                                    wtuv=o['wtuv'],
                                    robust=o['solutions_robust'],
                                    dt=o['dt_slow'],
                                    catcher=catcher,
                                    NChanSols=o['NChanSols'],
                                    # EvolutionSolFile=CurrentDDkMSSolName,
                                    MergeSmooth=o['smoothing'])
    
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
        
    if not o['skip_di']:
        # Compute the DD predict
        colname=o['colname']
        if o['do_wide']:
            colname ='DATA_SUB'
        separator("Compute DD Predict (full mslist)")
        ddf_image('Predict_DDS2',o['full_mslist'],cleanmode='SSD',
                applysols=o['apply_sols'][4],majorcycles=1,robust=o['image_robust'],colname=colname,peakfactor=0.01,
                automask=True,automask_threshold=o['thresholds'][1],normalization=o['normalize'][0],
                apply_weights=o['apply_weights'][0],use_weightspectrum=o['use_weightspectrum'],uvrange=uvrange,use_dicomodel=True,
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
                    niterkf=o['NIterKF'][5],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                    catcher=catcher,
                    dt=o['dt_di'],
                    DISettings=("CohJones","IFull","DD_PREDICT","DATA_DI_CORRECTED"))
        colname="DATA_DI_CORRECTED"

    # ###############################################
    # Apply phase and amplitude solutions and image again
    separator("Deconvolution AP (full mslist)")
    ddf_kw={}
    if o['msss_mode']:
        ddf_kw['cubemode']=True
        ddf_kw['smooth']=True

    if o['final_psf_arcsec'] is not None:
        ddf_kw['beamsize']=o['final_psf_arcsec']
        if o['final_psf_minor_arcsec'] is not None:
            if o['final_psf_pa_deg'] is None:
                die('If final minor axis is supplied, position angle must be supplied too')
            ddf_kw['beamsize_minor']=o['final_psf_minor_arcsec']
            ddf_kw['beamsize_pa']=o['final_psf_pa_deg']

    if not o['skip_di']:
        ImageName = 'image_full_ampphase_di'
    else:
        ImageName = 'image_full_ampphase'
    
    ddf_image(ImageName,o['full_mslist'],
              cleanmask=CurrentMaskName,
              cleanmode='SSD',ddsols=CurrentDDkMSSolName,
              applysols=o['apply_sols'][5],
              majorcycles=0,
              robust=o['final_robust'],
              colname=colname,use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              AllowNegativeInitHMP=True,
              peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][2],
              normalization=o['normalize'][1],uvrange=uvrange,smooth=True,
              apply_weights=o['apply_weights'][2],use_weightspectrum=o['use_weightspectrum'],catcher=catcher,**ddf_kw)

    if o['exitafter'] == 'fullampphase':
        warn('User specified exit after image_ampphase.')
        stop(2)
        
    separator("MakeMask")
    CurrentMaskName=make_mask(ImageName+'.app.restored.fits',10,external_mask=external_mask,catcher=catcher)

    separator("Finish Deconvolution AP (full mslist)")
    if not o['skip_di']:
        ImageName = 'image_full_ampphase_di_m'
    else:
        ImageName = 'image_full_ampphase_m'
    CurrentBaseDicoModelName=ddf_image(ImageName,o['full_mslist'],
                                       cleanmask=CurrentMaskName,
                                       reuse_psf=True,
                                       reuse_dirty=True,
                                       robust=o['final_robust'],
                                       cleanmode='SSD',ddsols=CurrentDDkMSSolName,
                                       applysols=o['apply_sols'][5],majorcycles=1,
                                       colname=colname,use_dicomodel=True,
                                       dicomodel_base=CurrentBaseDicoModelName,
                                       peakfactor=0.001,automask=True,
                                       automask_threshold=o['thresholds'][2],
                                       normalization=o['normalize'][1],uvrange=uvrange,
                                       apply_weights=o['apply_weights'][2],use_weightspectrum=o['use_weightspectrum'],catcher=catcher,
                                       AllowNegativeInitHMP=True,
                                       RMSFactorInitHMP=.5,
                                       MaxMinorIterInitHMP=10000,smooth=True,**ddf_kw)

    separator("MakeMask")
    CurrentMaskName=make_mask(ImageName+'.app.restored.fits',o['thresholds'][2],external_mask=external_mask,catcher=catcher)
    CurrentBaseDicoModelName=mask_dicomodel(ImageName+'.DicoModel',CurrentMaskName,ImageName+'_masked.DicoModel',catcher=catcher)
            
    separator("DD Calibration (full mslist)")
    CurrentDDkMSSolName=killms_data(ImageName,
                                    o['full_mslist'],'DDS3_full',
                                    colname=colname,
                                    clusterfile=ClusterFile,
                                    dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                    niterkf=o['NIterKF'][6],
                                    CovQ=0.1,
                                    uvrange=killms_uvrange,
                                    wtuv=o['wtuv'],
                                    robust=o['solutions_robust'],
                                    MergeSmooth=o['smoothing'],
                                    dt=o['dt_fast'],catcher=catcher)#,EvolutionSolFile=CurrentDDkMSSolName)

    if o['do_very_slow']:
        separator("Very slow amplitude smooth (full mslist)")
        CurrentDDkMSSolName_FastSmoothed=CurrentDDkMSSolName

        CurrentDDkMSSolName=killms_data(ImageName,
                                        o['full_mslist'],'DDS3_full_slow',
                                        colname=colname,
                                        SolverType="KAFCA",
                                        clusterfile=ClusterFile,
                                        dicomodel='%s.DicoModel'%CurrentBaseDicoModelName,
                                        uvrange=[o['uvmin_very_slow'],1000.],
                                        wtuv=o['wtuv'],
                                        robust=o['solutions_robust'],
                                        SkipSmooth=True,MergeSmooth=True,
                                        SigmaFilterOutliers=o['sigma_clip'],
                                        dt=o['dt_very_slow'],catcher=catcher,
                                        PreApplySols=CurrentDDkMSSolName_FastSmoothed)#,EvolutionSolFile=CurrentDDkMSSolName)

        CurrentDDkMSSolName="[%s,%s]"%(CurrentDDkMSSolName_FastSmoothed,CurrentDDkMSSolName)
    
    if o['low_psf_arcsec'] is not None:
        separator("Low-resolution image")
        # low-res image requested
        low_uvrange=[o['image_uvmin'],2.5*206.0/o['low_psf_arcsec']]
        if o['low_imsize'] is not None:
            low_imsize=o['low_imsize'] # allow over-ride
        else:
            low_imsize=o['imsize']*o['cellsize']/o['low_cell']
            # if mask-low exists then use it
        if os.path.isfile('bootstrap-mask-low.fits'):
            extmask='bootstrap-mask-low.fits'
            # can be empty, in which case recent versions of DDF throw
            # an error, so check and drop it if it is
            hdu=fits.open(extmask)
            if not np.any(hdu[0].data>0):
                warn('Bootstrap external mask is blank, using only internal masking')
                extmask=None
            hdu.close()
        else:
            extmask=None

        ddf_image('image_full_low',o['full_mslist'],
                  cleanmask=extmask,
                  cleanmode='SSD',ddsols=CurrentDDkMSSolName,
                  applysols=o['apply_sols'][6],
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
              applysols=o['apply_sols'][6],
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
            if os.path.isfile('image_dirin_SSD.app.restored.fits'):
                # Normal pipeline run.
                make_extended_mask('image_full_low_im.app.restored.fits','image_dirin_SSD.app.restored.fits',rmsthresh=o['extended_rms'],sizethresh=1500,rootname='full',rmsfacet=o['rmsfacet'])
            elif (not os.path.isfile('image_dirin_SSD.app.restored.fits')) and os.path.isfile('image_full_ampphase_di.app.restored.fits'):
                # Input model was given.
                make_extended_mask('image_full_low_im.app.restored.fits','image_full_ampphase_di.app.restored.fits',rmsthresh=o['extended_rms'],sizethresh=1500,rootname='full',rmsfacet=o['rmsfacet'],ds9region='image_full_ampphase_di_m.tessel.reg')
            # skip_di - we have this one instead
            elif (not os.path.isfile('image_dirin_SSD.app.restored.fits')) and os.path.isfile('image_full_ampphase.app.restored.fits'):
                # Input model was given.
                make_extended_mask('image_full_low_im.app.restored.fits','image_full_ampphase.app.restored.fits',rmsthresh=o['extended_rms'],sizethresh=1500,rootname='full',rmsfacet=o['rmsfacet'],ds9region='image_full_ampphase_m.tessel.reg')
            else:
                # Something may be wrong.
                die('Could not find the required products for the full-bw extended source mask!')
            report('Make_extended_mask returns')
        extmask='full-mask-low.fits'
        make_mask('image_full_low_im.app.restored.fits',3.0,external_mask=extmask,catcher=catcher)

        ddf_image('image_full_low_m',o['full_mslist'],
              cleanmask='image_full_low_im.app.restored.fits.mask.fits',
              cleanmode='SSD',ddsols=CurrentDDkMSSolName,
              applysols=o['apply_sols'][6],
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
            if os.path.isfile('image_full_ampphase_di.app.restored.fits'):
                make_external_mask(external_mask,'image_full_ampphase_di.app.restored.fits',use_tgss=True,clobber=False,extended_use='full-mask-high.fits')
            elif os.path.isfile('image_full_ampphase.app.restored.fits'):
                make_external_mask(external_mask,'image_full_ampphase.app.restored.fits',use_tgss=True,clobber=False,extended_use='full-mask-high.fits')

    # ##########################################################
    if o['exitafter'] == 'fulllow':
        warn('User specified exit after full low.')
        stop(2)

    # before starting the final image, run the download thread if needed
    if o['method'] is not None:
        separator('Offset image downloads')
        report('Checking if optical catalogue download is required')
        from get_cat import get_cat, download_required
        if download_required(o['method']):
            download_thread = threading.Thread(target=get_cat, args=(o['method'],))
            download_thread.start()
        else:
            warn('All data present, skipping download')
            download_thread = None


    ImageName = 'image_full_ampphase_di_m.NS'

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
    ddf_image(ImageName,o['full_mslist'],
              cleanmask=CurrentMaskName,
              reuse_psf=False,
              cleanmode='SSD',ddsols=CurrentDDkMSSolName,
              applysols=o['apply_sols'][6],majorcycles=1,robust=o['final_robust'],
              colname=colname,use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              AllowNegativeInitHMP=True,
              peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][2],
              normalization=o['normalize'][1],uvrange=uvrange,smooth=True,
              apply_weights=o['apply_weights'][2],use_weightspectrum=o['use_weightspectrum'],catcher=catcher,RMSFactorInitHMP=1.,
              PredictSettings=("Clean","DD_PREDICT"),
              **ddf_kw)

    # check for the offset files
    if o['method'] is not None:
        separator('Offset correction')
        # have we got the catalogue?
        if download_thread is not None and download_thread.is_alive():
            warn('Waiting for background download thread to finish...')
            download_thread.join()
        # maybe the thread died, check the files are there
        if download_required(o['method']):
            warn('Retrying download for some or all of the catalogue')
            try:
                get_cat(o['method'])
            except RuntimeError:
                die('Failed to download catalogue with method '+o['method'])

        # we should now have the catalogue, find the offsets
        facet_offset_file='facet-offset.txt'
        if o['restart'] and os.path.isfile(facet_offset_file):
            warn('Offset file already exists, not running offsets.py')
        else:
            run('offsets.py '+' '.join(sys.argv[1:]),log=None)

        # apply the offsets
        ddf_shift(ImageName,facet_offset_file,options=o,catcher=catcher)
    else:
        facet_offset_file=None
            
    spectral_mslist=None
    if o['spectral_restored']:
        import do_spectral_restored
        separator('Spectral restored images')
        spectral_mslist=do_spectral_restored.do_spectral_restored(colname,
                                                  CurrentMaskName,
                                                  CurrentBaseDicoModelName,
                                                  CurrentDDkMSSolName,
                                                  uvrange,
                                                  ddf_kw,
                                                  facet_offset_file,
                                                  options=o,
                                                  catcher=catcher)

    if o['polcubes']:
        from do_polcubes import do_polcubes
        separator('Stokes Q and U cubes')
        cthreads=[]
        flist=[]

        if o['split_polcubes']:
            cubefiles=['image_full_low_StokesQ.cube.dirty.fits','image_full_low_StokesQ.cube.dirty.corr.fits','image_full_low_StokesU.cube.dirty.fits','image_full_low_StokesU.cube.dirty.corr.fits']
        else:
            cubefiles=['image_full_low_QU.cube.dirty.fits','image_full_low_QU.cube.dirty.corr.fits']
        if o['restart'] and os.path.isfile(cubefiles[0]+'.fz') and os.path.isfile(cubefiles[1]+'.fz'):
            warn('Compressed low QU cube product exists, not making new images')
        else:
            do_polcubes(colname,CurrentDDkMSSolName,low_uvrange,'image_full_low',ddf_kw,beamsize=o['low_psf_arcsec'],imsize=low_imsize,cellsize=o['low_cell'],robust=o['low_robust'],options=o,catcher=catcher)
            if o['compress_polcubes']:
                for cubefile in cubefiles:
                    if o['restart'] and os.path.isfile(cubefile+'.fz'):
                        warn('Compressed cube file '+cubefile+'.fz already exists, not starting compression thread')
                    else:
                        report('Starting compression thread for '+cubefile)
                        thread = threading.Thread(target=compress_fits, args=(cubefile,o['fpack_q']))
                        thread.start()
                        cthreads.append(thread)
                        flist.append(cubefile)
        if o['split_polcubes']:
            cubefiles=['image_full_vlow_StokesQ.cube.dirty.fits','image_full_vlow_StokesQ.cube.dirty.corr.fits','image_full_vlow_StokesU.cube.dirty.fits','image_full_vlow_StokesU.cube.dirty.corr.fits']
        else:
            cubefiles=['image_full_vlow_QU.cube.dirty.fits','image_full_vlow_QU.cube.dirty.corr.fits']
        if o['restart'] and os.path.isfile(cubefiles[0]+'.fz') and os.path.isfile(cubefiles[1]+'.fz'):
            warn('Compressed vlow QU cube product exists, not making new images')
        else:
            vlow_uvrange=[o['image_uvmin'],1.6]
            do_polcubes(colname,CurrentDDkMSSolName,vlow_uvrange,'image_full_vlow',ddf_kw,beamsize=o['vlow_psf_arcsec'],imsize=o['vlow_imsize'],cellsize=o['vlow_cell'],robust=o['vlow_robust'],options=o,catcher=catcher)
            if o['compress_polcubes']:
                for cubefile in cubefiles:
                    if o['restart'] and os.path.isfile(cubefile+'.fz'):
                        warn('Compressed cube file '+cubefile+'.fz already exists, not starting compression thread')
                    else:
                        report('Starting compression thread for '+cubefile)
                        thread = threading.Thread(target=compress_fits, args=(cubefile,o['fpack_q']))
                        thread.start()
                        cthreads.append(thread)
                        flist.append(cubefile)
        
    if o['stokesv']:
        separator('Stokes V image')
        ddf_image('image_full_high_stokesV',o['full_mslist'],
                  cleanmode='SSD',ddsols=CurrentDDkMSSolName,
                  applysols=o['apply_sols'][6],stokes='IV',
                  AllowNegativeInitHMP=True,
                  majorcycles=0,robust=o['final_robust'],
                  colname=colname,use_dicomodel=False,
                  uvrange=uvrange,cellsize=o['cellsize'],
                  peakfactor=0.001,
                  smooth=True,automask=True,automask_threshold=5,normalization=o['normalize'][2],
                  catcher=catcher,**ddf_kw)

    if o['polcubes'] and o['compress_polcubes']:
        # cthreads and flist exist
        for thread in cthreads:
            if thread.is_alive():
                warn('Waiting for a compression thread to finish')
                thread.join()
        if o['delete_compressed']:
            for f in flist:
                if os.path.isfile(f+'.fz'):
                    warn('Deleting compressed file %s' % f)
                    os.remove(f)
                else:
                    die('compressed files do not exist, compression must have failed')

    if o['do_dynspec']:
        separator('Dynamic spectra')

        if o['bright_threshold'] is not None and o['method'] is not None:
            warn('Finding bright sources from offsets list')
            from find_bright_offset_sources import find_bright
            bright_exists=find_bright(cutoff=o['bright_threshold'])
        LastImage="image_full_ampphase_di_m.NS.int.restored.fits"
        m=MSList(o['full_mslist'])
        uobsid = set(m.obsids)
    
        for obsid in uobsid:
            LastImageV="image_full_high_stokesV.dirty.corr.fits"
            warn('Running ms2dynspec for obsid %s' % obsid)
            umslist='mslist-%s.txt' % obsid
            print('Writing temporary ms list',umslist)
            with open(umslist,'w') as file:
                for ms,ob in zip(m.mss,m.obsids):
                    if ob==obsid:
                        file.write(ms+'\n')

            g=glob.glob('DynSpec*'+obsid+'*')
            if len(g)>0:
                warn('DynSpecs results directory %s already exists, skipping DynSpecs' % g[0])
            else:
                DicoFacetName="%s.DicoFacet"%LastImage.split(".int.restored.fits")[0]
                runcommand="ms2dynspec.py --ms %s --data %s --model DD_PREDICT --sols %s --rad 2. --imageI %s --imageV %s --LogBoring %i --SolsDir %s --BeamModel LOFAR --BeamNBand 1 --DicoFacet %s  --noff 100 --nMinOffPerFacet 5 --CutGainsMinMax 0.1,1.5 --SplitNonContiguous 1 --SavePDF 1 --FitsCatalog ${DDF_PIPELINE_CATALOGS}/dyn_spec_catalogue_addedexo_addvlotss.fits"%(umslist,colname,CurrentDDkMSSolName,LastImage,LastImageV,o['nobar'],o["SolsDir"],DicoFacetName)

                
                if o['bright_threshold'] is not None:
                    runcommand+=' --srclist brightlist.csv'
                run(runcommand,dryrun=o['dryrun'],log=logfilename('ms2dynspec.log'),quiet=o['quiet'])
                if use_database():
                    ingest_dynspec(obsid)


    if o['compress_ms']:
        separator('Compressing MS for archive')
        if o['skip_di']:
            os.system('archivems.sh . '+o['colname'])
        else:
            os.system('archivems.sh . DATA_DI_CORRECTED')
                
    separator('Write summary and tidy up')
    summary(o)

    # Clear caches if option set
    if o['clearcache_end']:
        extras=[]
        if spectral_mslist is not None:
            extras+=spectral_mslist
        if o['polcubes']:
            extras+=glob.glob('stokes-mslist*.txt')
        full_clearcache(o,extras=extras)
    
    if use_database():
        update_status(None,'Complete',time='end_date',av=4)
        
    return

if __name__=='__main__':
    # Main loop
    report('Welcome to ddf-pipeline, version '+__version__)
    if len(sys.argv)<2:
        warn('pipeline.py must be called with at least one parameter file or a command-line\noption list.\nE.g "pipeline.py example.cfg second_example.cfg --solutions-robust=0.1"\nSee below for a complete list of possible options with their default values.')
        print_options(option_list)
        sys.exit(1)

    o=options(sys.argv[1:],option_list)
    if MyPickle is not None:
        MyPickle.Save(o, "ddf-pipeline.last")

    main(o)
