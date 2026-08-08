#!/usr/bin/env python
"""
ddf-pipeline, a pipeline for LOFAR data reduction
Copyright (C) 2017-2026 Martin Hardcastle (mjh@extragalactic.info) and others

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
        ts = f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}"
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
        f.writelines(f"{k:<20} : {value}\n" for k, value in o.items())

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
    if not os.path.isfile(sols):
        if o['dryrun']:
            return 0.0, 1.0
        raise FileNotFoundError(f'Solutions file not found: {sols}')
    print(f'Reading {sols}')
    S=np.load(sols)
    t = np.concatenate([S["Sols"]["t0"],S["Sols"]["t1"]])
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
    runcommand = f'DDF.py {imagename}.parset --Misc-ConserveMemory=1  --Output-Name={imagename}_shift --Output-Mode=RestoreAndShift --Output-ShiftFacetsFile={shiftfile} --Predict-InitDicoModel {dicomodel} --Cache-SmoothBeam=force --Log-Memory 1 --Cache-Dir={cache_dir}'
    if 'Misc-IgnoreDeprecationMarking' in keywords:
        runcommand += ' --Misc-IgnoreDeprecationMarking=1'
    
    fname = imagename+'_shift.app.facetRestored.fits'
    if options['restart'] and os.path.isfile(fname):
        warn(f'File {fname} already exists, skipping DDF-shift step')
        if verbose:
            print('would have run',runcommand)
    else:
         run(runcommand,dryrun=options['dryrun'],log=logfilename(f'DDF-{imagename}_shift.log',options=options),quiet=options['quiet'])

def ddf_image(
    imagename,
    mslist,

    # Data
    colname='CORRECTED_DATA',

    # Solutions
    ddsols=None,
    applysols=None,

    # Image
    cleanmode=None,
    robust=0,
    imsize=None,
    cellsize=None,
    uvrange=None,
    beamsize=None,
    beamsize_minor=None,
    beamsize_pa=None,
    phasecenter=None,
    do_decorr=None,
    HMPsize=None,
    apply_weights=True,
    use_weightspectrum=False,
    clusterfile=None,
    normalization=None,
    cubemode=False,
    polcubemode=False,
    channels=None,
    startchan=None,
    endchan=None,
    stokes=None,
    freq_nband=2,
    peakfactor=0.1,
    rms_factor=3.0,
    threshold=None,

    # Masking
    cleanmask=None,
    automask=True,
    automask_threshold=None,
    
    # Control
    majorcycles=3,
    reuse_psf=False,
    reuse_dirty=False,
    dirty_from_resid=False,
    conditional_clearcache=False,
    verbose=False,
    saveimages=None,
    predict_column=None,
    smooth=False,
    noweights=False,
    catcher=None,

    # HMP settings
    RMSFactorInitHMP=1.0,
    MaxMinorIterInitHMP=10000,
    AllowNegativeInitHMP=False,
    OuterSpaceTh=None,

    # WSCMS settings
    # allows the option for per pipeline-step scale adjustments
    wscms_MultiScaleBias=None,
    wscms_Scales=None,
    wscms_MaxScale=None,
    wscms_rms_factor=None,
    wscms_peakfactor=None,
    wscms_automask_rms_factor=None,
    wscms_allownegative=None,
    wscms_flux_threshold=None,
    wscms_neg_max_scale=None,
    wscms_NSubMinorIter=250,
    wscms_SubMinorPeakFact=0.85,
    wscms_remove_1px_scale=True,

    # Input model
    use_dicomodel=False,
    dicomodel_base=None,

    # Other
    options=None,
    PredictSettings=None,
    STEP=0):

    if catcher: catcher.check()

    # saveimages lists _additional_ images to save
    if saveimages is None:
        saveimages = ''
    saveimages += 'onNeds'
    if options is None:
        options = o # attempt to get global if it exists

    keywords = parse_parset([os.environ['DDF_DIR']+'/DDFacet/DDFacet/Parset/DefaultParset.cfg'],use_headings=True)

    for key in ('wscms_rms_factor', 'wscms_peakfactor', 'wscms_automask_rms_factor', 'wscms_allownegative', 'wscms_max_scale', 'wscms_flux_threshold', 'thresholds', 'use_external_mask'):
        if STEP > len(options[key]) - 1:
            raise ValueError(f"STEP={STEP} out of range for {key} (length {len(options[key])})")

    # pull defaults from parset if not specified
    if HMPsize            is None: HMPsize            = options['HMPsize']
    if cleanmode          is None: cleanmode          = options['cleanmode']
    if do_decorr          is None: do_decorr          = options['do_decorr']
    if beamsize           is None: beamsize           = options['psf_arcsec']
    if imsize             is None: imsize             = options['imsize']
    if cellsize           is None: cellsize           = options['cellsize']
    if automask_threshold is None: automask_threshold = options['thresholds'][STEP]

    # pull default WSCMS settings from parset if not specified
    if wscms_MultiScaleBias      is None: wscms_MultiScaleBias      = options['wscms_multiscale_bias']
    if wscms_Scales              is None: wscms_Scales              = options['wscms_scales']
    if wscms_MaxScale            is None: wscms_MaxScale            = options['wscms_max_scale'][STEP]
    if wscms_rms_factor          is None: wscms_rms_factor          = options['wscms_rms_factor'][STEP]
    if wscms_peakfactor          is None: wscms_peakfactor          = options['wscms_peakfactor'][STEP]
    if wscms_automask_rms_factor is None: wscms_automask_rms_factor = options['wscms_automask_rms_factor'][STEP]
    if wscms_allownegative       is None: wscms_allownegative       = options['wscms_allownegative'][STEP]
    if wscms_flux_threshold      is None: wscms_flux_threshold      = options['wscms_flux_threshold'][STEP]
    if wscms_neg_max_scale       is None: wscms_neg_max_scale       = options['wscms_neg_max_scale']

    cache_dir=find_cache_dir(options)

    if majorcycles>0:
        fname = imagename+'.app.restored.fits'
    else:
        fname = imagename+'.dirty.fits'

    if PredictSettings is not None and PredictSettings[0]=="Predict":
        fname = f"_has_predicted_OK.{imagename}.info"

    runcommand = f"DDF.py --Misc-ConserveMemory=1 --Output-Name={imagename} --Data-MS={mslist}  --Data-ColName {colname} --Parallel-NCPU={options['NCPU_DDF']} --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter={majorcycles} --Deconv-Mode {cleanmode} --Beam-Model=LOFAR --Weight-Robust {robust} --Image-NPix={imsize} --CF-wmax 50000 --CF-Nw 100 --Output-Also {saveimages} --Image-Cell {float(cellsize)} --Facets-NFacets=11 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1  --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir={cache_dir} --Cache-DirWisdomFFTW={cache_dir} --Debug-Pdb=never --Log-Memory 1"

    if RMSFactorInitHMP is not None:    runcommand += f" --GAClean-RMSFactorInitHMP {RMSFactorInitHMP}"
    if MaxMinorIterInitHMP is not None: runcommand += f" --GAClean-MaxMinorIterInitHMP {MaxMinorIterInitHMP}"
    if AllowNegativeInitHMP:            runcommand += " --GAClean-AllowNegativeInitHMP True"
    if OuterSpaceTh is not None:        runcommand += f" --HMP-OuterSpaceTh {OuterSpaceTh}"

    runcommand += f' --DDESolutions-SolsDir={options["SolsDir"]}'
    runcommand += ' --Cache-Weight=reset'

    if 'Beam-PhasedArrayMode' in keywords: # incompatible change
        runcommand += ' --Beam-PhasedArrayMode=A'
    else:
        runcommand += ' --Beam-LOFARBeamMode=A'
    
    if 'Misc-IgnoreDeprecationMarking' in keywords:
        runcommand += ' --Misc-IgnoreDeprecationMarking=1'

    if 'Beam-At' in keywords:
        runcommand += f' --Beam-At={options["beam_at"]}'
        
    if PredictSettings is None:
        runcommand += " --Output-Mode=Clean"
    else:
        if len(PredictSettings) == 2:
            runcommand += f" --Output-Mode={PredictSettings[0]} --Predict-ColName {PredictSettings[1]}"
        elif len(PredictSettings) == 3:
            runcommand += f" --Output-Mode={PredictSettings[0]} --Predict-ColName {PredictSettings[1]} --Predict-MaskSquare [0,{PredictSettings[2]}]"
        else:
            raise RuntimeError(f'PredictSettings has the wrong dimensions {PredictSettings}')

    if beamsize_minor is not None:
        runcommand += f' --Output-RestoringBeam {beamsize},{beamsize_minor},{beamsize_pa}'
    elif beamsize is not None:
        runcommand += f' --Output-RestoringBeam {beamsize}'
    
    if apply_weights:
        runcommand += ' --Weight-ColName="IMAGING_WEIGHT"'
    else:
        if not use_weightspectrum:
            runcommand += ' --Weight-ColName="None"'
        else:
            runcommand += ' --Weight-ColName="WEIGHT_SPECTRUM"'

    if cubemode or polcubemode:
        # number of channels equals number of distinct freqs in data
        freqs=[]
        mss=[l.rstrip() for l in open(mslist).readlines()]
        for ms in mss:
            with pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False) as t:
                freq=t[0]['REF_FREQUENCY']
                if freq not in freqs:
                    freqs.append(freq)
        channels=len(freqs)
        
        if cubemode:    runcommand += f' --Output-Cubes I --Freq-NBand={channels}'
        if polcubemode: runcommand += f' --Output-Cubes=dD --RIME-PolMode=IQU --Output-Mode=Dirty  --Freq-NBand={channels} --Selection-ChanStart={startchan} --Selection-ChanEnd={endchan}'
    else:
        #if not cubemode and not polcubemode
        runcommand += f' --Freq-NBand={freq_nband}'
    
    if stokes:    runcommand += f' --RIME-PolMode={stokes} --Output-Mode=Dirty'
    if do_decorr: runcommand += ' --RIME-DecorrMode=FT'

    ###########################
    #### SSD / SSD2 / SSD3 ####
    ###########################
    if cleanmode in ('SSD', 'SSD2', 'SSD3'):
        # parameters shared between SSD versions
        runcommand += f' --Deconv-RMSFactor={rms_factor} --Deconv-PeakFactor={peakfactor} --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --SSDClean-NEnlargeData 0'
        if automask:                   runcommand += f' --Mask-Auto=1 --Mask-SigTh={automask_threshold:.2f}' 
        if cleanmask is not None:      runcommand += f' --Mask-External={cleanmask}'
        if options['use_splitisland']: runcommand += f' --SSDClean-MaxIslandSize={options["splitisland_size"]}'

        # SSD2 specific parameters
        if cleanmode == 'SSD2':
            runcommand += ''

        # SSD3 specific parameters
        if cleanmode == 'SSD3':
             runcommand += ' --SSD3-InitType ["HMP_0-100","MultiSlice:Orieux"]'
             runcommand += ' --SSD3-NLookBackModels 5'
             runcommand += ' --SSD3-NLastCycleDeconvAll 2'
             runcommand += ' --SSD3-AllowFacetOverlap 0'
             runcommand += ' --SSD3-RunSimpleClean 0'
             runcommand += ' --SSD3-PropagatePrevGen 0'
             runcommand += ' --SSD3-AlphaScaleModel 1'
             runcommand += ' --SSDClean-SSDCostFunc ["Chi2"]'
             runcommand += ' --GAClean-NMaxGen 30'
             runcommand += ' --GAClean-NSourceKin 30'
             runcommand += ' --Mask-ThFilterRFI 3'
             runcommand += ' --MultiSliceDeconv-HyperSmooth 3.5'
             
    ########################
    #### WSCMS / WSCMS2 ####
    ########################
    elif cleanmode in ('WSCMS', 'WSCMS2'):
        # parameters shared between WSCMS versions
        runcommand += f' --WSCMS-MultiScale=1 --WSCMS-MultiScaleBias={wscms_MultiScaleBias} --WSCMS-NSubMinorIter={wscms_NSubMinorIter} --WSCMS-SubMinorPeakFact={wscms_SubMinorPeakFact} --Deconv-RMSFactor={wscms_rms_factor} --Deconv-PeakFactor={wscms_peakfactor}'
        
        # wscms max scale is not aware of cellsizes, it only knows pixels
        # to prevent divergence at the wide, low, and vlow steps, scale by the standard cellsize (usually 1.5"/px), and ensure no duplicates are created due to rounding.
        # we also exclude the scale px=1 (default:True), to avoid pointsource-fitting divergence
        if wscms_MaxScale is not None:
            wscms_MaxScale_px = int(float(wscms_MaxScale) / cellsize * options['cellsize'])
            runcommand += f' --WSCMS-MaxScale={wscms_MaxScale_px}'
        if wscms_Scales is not None:
            scaled = {int(float(s) / cellsize * options['cellsize']) for s in wscms_Scales}
            # have to make all scales odd before deduplication to avoid double values when using dense scale bank
            scaled = {s + 1 if s % 2 == 0 and s > 0 else s for s in scaled}
            if wscms_remove_1px_scale: scaled -= {1}
            wscms_Scales_px = str(sorted(scaled)).replace(" ", "")
            runcommand += f' --WSCMS-Scales={wscms_Scales_px}'
        if wscms_neg_max_scale is not None:
            # also scale wscms_neg_max_scale_px by cellsize to keep it consistent
            wscms_neg_max_scale_px = int(float(wscms_neg_max_scale) / cellsize * options['cellsize'])
            runcommand += f' --WSCMS-NegMaxScale={wscms_neg_max_scale_px}'

        # dedicated automasking inputs
        runcommand += f' --WSCMS-AutoMask={automask}'
        if automask:
            runcommand += f' --WSCMS-AutoMaskRMSFactor={wscms_automask_rms_factor}'
            runcommand += ' --WSCMS-AutoMaskForceLast=False'

        runcommand += f' --Deconv-AllowNegative={int(wscms_allownegative)}'
        runcommand += f' --Deconv-FluxThreshold={wscms_flux_threshold}'

        if cleanmask is not None: runcommand += f' --Mask-External={cleanmask}'
        runcommand += ' --Mask-FluxImageType=ModelConv'


    if clusterfile is not None:
        runcommand += f' --Facets-CatNodes={clusterfile}'
    
    if applysols is not None:
        if normalization is not None:
            if normalization[:3]=='Abs':
                normalization='Mean'+normalization # backward compat. hack
            runcommand += ' --DDESolutions-GlobalNorm='+normalization
        runcommand += f' --DDESolutions-DDModeGrid={applysols} --DDESolutions-DDModeDeGrid={applysols} --DDESolutions-DDSols={ddsols}'
    if use_dicomodel:
        if dicomodel_base is not None:
            runcommand += f' --Predict-InitDicoModel={dicomodel_base}.DicoModel'
        else:
            raise RuntimeError('use_dicomodel is set but no dicomodel supplied')
        
    if threshold is not None:
        runcommand += f' --Deconv-FluxThreshold={threshold}'
    if uvrange is not None:
        runcommand += f' --Selection-UVRangeKm=[{uvrange[0]},{uvrange[1]}]'
    if dirty_from_resid and reuse_dirty:
        raise RuntimeError('Cannot combine reuse_dirty and dirty_from_resid')

    # possible that crashes could destroy the cache, so need to check
    if dirty_from_resid and os.path.exists(cache_dir+'/'+mslist+'.ddfcache/LastResidual'):
        runcommand += ' --Cache-Dirty forceresidual'
    if reuse_dirty and os.path.exists(cache_dir+'/'+mslist+'.ddfcache/Dirty'):
        runcommand += ' --Cache-Dirty forcedirty'
    if reuse_psf and os.path.exists(cache_dir+'/'+mslist+'.ddfcache/PSF'):
        runcommand += ' --Cache-PSF force'

    if HMPsize is not None:
        runcommand += f' --GAClean-MinSizeInit={HMPsize}'

    if options['nobar']:
        runcommand += ' --Log-Boring=1'

    if smooth:
        runcommand += ' --Beam-Smooth=1'

    if predict_column is not None:
        runcommand += f' --Predict-ColName={predict_column}'
        
    if phasecenter is not None:
        runcommand += f' --Image-PhaseCenterRADEC=[{phasecenter[0]},{phasecenter[1]}]'

    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping DDF step')
        if verbose:
            print('would have run',runcommand)
    else:
        if conditional_clearcache:
            clearcache(mslist,options)
        run(runcommand,dryrun=options['dryrun'],log=logfilename('DDF-'+imagename+'.log',options=options),quiet=options['quiet'])

        # Ugly way to see if predict has been already done
        if PredictSettings is not None and not options['dryrun']:
            os.system(f"touch {fname}")
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
        warn(f'File {OutClusterCat} already exists, skipping clustering step')
        return

    if ".app.restored.fits" not in imagename:
        raise RuntimeError('Input image should be an apparant restored image')

    if options is None:
        options=o # attempt to get global if it exists

    if use_makemask_products:
        runcommand = f"MakeCatalog.py --RestoredIm {imagename} --rmsmean_map [Noise.mean.fits,Noise.fits]"
    else:
        runcommand = f"MakeCatalog.py --RestoredIm {imagename}"
    run(runcommand,dryrun=options['dryrun'],log=logfilename('MakeCatalog-'+imagename+'.log',options=options),quiet=options['quiet'])

    Name=imagename.split(".app.restored.fits")[0]

    #runcommand="ClusterCat.py --SourceCat %s.app.restored.pybdsm.srl.fits --AvoidPolygons MaskDiffuse.pickle --NGen 100 --FluxMin 0.1"%Name
    filenames=[Name+'.app.restored.pybdsm.srl.fits',Name+'.app.restored.pybdsf.srl.fits']
    for filename in filenames:
        if os.path.isfile(filename):
            break
    else:
        if not options['dryrun']:
            die('Catalogue file does not exist!')
        filename=filenames[0]
    if use_makemask_products:
        runcommand = f"ClusterCat.py --SourceCat {filename} --AvoidPolygons MaskDiffuse.pickle --DoPlot=0 --NGen 100 --NCPU {options['NCPU_DDF']}"
    else:
        runcommand = f"ClusterCat.py --SourceCat {filename} --DoPlot=0 --NGen 100 --NCPU {options['NCPU_DDF']}"
    if OutClusterCat is not None:
        runcommand += f" --OutClusterCat {OutClusterCat}"
    runcommand += f" --NCluster {options['ndir']}"
    run(runcommand,dryrun=options['dryrun'],log=logfilename('MakeCluster-'+imagename+'.log',options=options),quiet=options['quiet'])

def make_mask(imagename,thresh,verbose=False,options=None,external_mask=None,catcher=None,OutMaskExtended=None):
    if catcher: catcher.check()

    # mask_use specifies a mask file to use
    if options is None:
        options=o # attempt to get global

    fname=imagename+'.mask.fits'
    if options['dryrun']: return fname
    runcommand = f"MakeMask.py --RestoredIm={imagename} --Th={thresh} --Box=50,2"
    if OutMaskExtended is not None:
        runcommand += f" --OutMaskExtended {OutMaskExtended} --OutNameNoiseMap Noise"


        
    if options['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeMask step')
        if verbose:
            print('Would have run',runcommand)
    else:
        run(runcommand,dryrun=options['dryrun'],log=logfilename('MM-'+imagename+'.log',options=options),quiet=options['quiet'])
        if external_mask is not None:
            if isinstance(external_mask, (list, tuple)):
                for mask in external_mask:
                    merge_mask(fname,mask,fname)
            else:
                merge_mask(fname,external_mask,fname)
    return fname

def killms_data(imagename,mslist,outsols,clusterfile=None,colname='CORRECTED_DATA',niterkf=6,dicomodel=None,
                uvrange=None,wtuv=None,robust=None,catcher=None,dt=None,options=None,
                SolverType="KAFCA",PolMode="Scalar",MergeSmooth=False,NChanSols=None,
                DISettings=None,EvolutionSolFile=None,CovQ=0.1,InterpToMSListFreqs=None,
                SkipSmooth=False,PreApplySols=None,SigmaFilterOutliers=None,UpdateWeights=None):

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
        clipcal_checkname=checkname.replace('.sols.npz','.clipcal_done')
        rootfilename=outsols.split('/')[-1]
        f_=f.replace("/","_")


        #checkname=f+'/killMS.'+outsols+'.sols.npz'
        if options['restart'] and os.path.isfile(checkname):

            warn('Solutions file '+checkname+' already exists, not running killMS step')
            
        else:
            runcommand = f"kMS.py --MSName {f} --SolverType {SolverType} --PolMode {PolMode} --BaseImageName {imagename} --NIterKF {niterkf} --CovQ {CovQ} --LambdaKF={options['LambdaKF']} --NCPU {options['NCPU_killms']} --OutSolsName {outsols} --InCol {colname}"

            # check for option to stop pdb call and use it if present
            
            if 'DebugPdb' in keywords:
                runcommand += ' --DebugPdb=0'
                
            if robust is None:
                runcommand += ' --Weighting Natural'
            else:
                runcommand += f' --Weighting Briggs --Robust={robust}'
            if UpdateWeights is not None:
                runcommand += f' --UpdateWeights={UpdateWeights}'               
            if uvrange is not None:
                if wtuv is not None:
                    runcommand += f' --WTUV={wtuv} --WeightUVMinMax={uvrange[0]},{uvrange[1]}'
                else:
                    runcommand += f' --UVMinMax={uvrange[0]},{uvrange[1]}'
            if options['nobar']:
                runcommand += ' --DoBar=0'

            runcommand += f' --SolsDir={options["SolsDir"]}'
            
            if PreApplySols:
                runcommand += f' --PreApplySols=[{PreApplySols}]'
                
            if DISettings is None:
                if NChanSols is None:
                    NChanSols=1 # reproduce old behaviour
                runcommand += f' --NChanSols {NChanSols}'
                runcommand += ' --BeamMode LOFAR'
                if 'PhasedArrayMode' in keywords: # incompatible change
                    runcommand += ' --PhasedArrayMode=A'
                else:
                    runcommand += ' --LOFARBeamMode=A'
                runcommand += ' --DDFCacheDir='+cache_dir
                if 'BeamAt' in keywords:
                    runcommand += f' --BeamAt={options["beam_at"]}'

                if clusterfile is not None:
                    runcommand += ' --NodesFile '+clusterfile
                if dicomodel is not None:
                    runcommand += ' --DicoModel '+dicomodel
                if EvolutionSolFile is not None:
                    runcommand += ' --EvolutionSolFile '+EvolutionSolFile
                if dt is not None:
                    runcommand += f' --dt {dt}'
            else:
                runcommand += f" --SolverType {DISettings[0]} --PolMode {DISettings[1]} --SkyModelCol {DISettings[2]} --OutCol {DISettings[3]} --ApplyToDir 0"
                _,_,ModelColName,_=DISettings
                if not options['dryrun'] or os.path.isfile(f):
                    _,dt_give,_,n_df_give=give_dt_dnu(f,
                                            DataCol=colname,
                                            ModelCol=ModelColName,
                                            T=10.)
                    if dt is None:
                        dt=dt_give
                    if NChanSols is None:
                        NChanSols=n_df_give
                else:
                    if dt is None: dt=1.0
                    if NChanSols is None: NChanSols=1
                runcommand += f" --dt {dt+1e-4} --NChanSols {NChanSols}"
                
                
            run(runcommand,dryrun=options['dryrun'],log=logfilename('KillMS-'+f_+'_'+rootfilename+'.log',options=options),quiet=options['quiet'])

        # Clip anyway - on IMAGING_WEIGHT by default
        if DISettings is not None:
            ClipCol=DISettings[-1]
        else:
            ClipCol=colname
        runcommand = f"ClipCal.py --MSName {f} --ColName {ClipCol}"

        if options['restart'] and os.path.isfile(clipcal_checkname):
            warn('ClipCal done file '+clipcal_checkname+' already exists, not running ClipCal step')
        else:
            run(runcommand,dryrun=options['dryrun'],log=logfilename('ClipCal-'+f_+'_'+rootfilename+'.log',options=options),quiet=options['quiet'])
            if not options['dryrun']:
                os.system(f"touch {clipcal_checkname}")

    if MergeSmooth:
        outsols=smooth_solutions(mslist,outsols,catcher=None,dryrun=options['dryrun'],InterpToMSListFreqs=InterpToMSListFreqs,
                                 SkipSmooth=SkipSmooth,SigmaFilterOutliers=SigmaFilterOutliers,options=options)
        



    return outsols

def compress_fits(filename,q):
    command = f'fpack -q {q} {filename}'
    run(command, dryrun=o['dryrun'])
    
def make_model(maskname,imagename,catcher=None):
    # returns True if the step was run, False if skipped
    if catcher: catcher.check()

    fname=imagename+'.npy'
    if o['restart'] and os.path.isfile(fname):
        warn('File '+fname+' already exists, skipping MakeModel step')
        return False
    else:
        runcommand = f"MakeModel.py --MaskName={maskname} --BaseImageName={imagename} --NCluster={o['ndir']} --DoPlot=0"
        run(runcommand,dryrun=o['dryrun'],log=logfilename('MakeModel-'+maskname+'.log'),quiet=o['quiet'])
        return True

def mask_dicomodel(indico,maskname,outdico,catcher=None,filterneg=True):
    if catcher: catcher.check()

    if o['restart'] and os.path.isfile(outdico):
        warn('File '+outdico+' already exists, skipping MaskDicoModel step')
    else:
        runcommand = f"MaskDicoModel.py --InDicoModel={indico} --OutDicoModel={outdico}"
        if maskname:
            runcommand += f" --MaskName={maskname}"
        if filterneg:
            runcommand += " --FilterNegComp=1"
        run(runcommand,dryrun=o['dryrun'],log=logfilename('MaskDicoModel-'+outdico+'.log'),quiet=o['quiet'])
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

def smooth_solutions(mslist,ddsols,catcher=None,dryrun=False,InterpToMSListFreqs=None,SkipSmooth=False,SigmaFilterOutliers=None,options=None):
    if options is None:
        options=o
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    full_sollist = []
    start_times = []
    SolsDir=options["SolsDir"]
    if SolsDir is None or SolsDir=="":
        for fname in filenames:
            solname =fname+'/killMS.'+ddsols+'.sols.npz'
            t0,t1 = get_solutions_timerange(solname)
            start_times.append(t0)
            full_sollist.append(solname)
    else:
        for fname in filenames:
            MSName=os.path.abspath(fname).split("/")[-1]
            solname =os.path.abspath(SolsDir)+"/"+MSName+'/killMS.'+ddsols+'.sols.npz'
            t0,t1 = get_solutions_timerange(solname)
            start_times.append(t0)
            full_sollist.append(solname)

    Ustart_times = np.unique(start_times)

    for start_time in Ustart_times:
        if not dryrun:
            with open(f'solslist_{start_time:.2f}.txt', 'w') as f:
                for i in range(0, len(full_sollist)):
                    if start_times[i] == start_time:
                        solname = full_sollist[i]
                        f.write(f'{solname}\n')
        
        checkname=f'{ddsols}_{start_time:.2f}_merged.npz'
        if options['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running MergeSols step')
        else:
            ss=f'MergeSols.py --SolsFilesIn=solslist_{start_time:.2f}.txt --SolFileOut={checkname}'
            if SigmaFilterOutliers:
                ss+=f" --SigmaFilterOutliers {SigmaFilterOutliers}"
            run(ss,dryrun=dryrun)
            
        checkname=f'{ddsols}_{start_time:.2f}_smoothed.npz'
        if options['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running SmoothSols step')
        elif SkipSmooth:
            warn('Skipping smoothing Solutions file')
        else:
            run(f'SmoothSols.py --SolsFileIn={ddsols}_{start_time:.2f}_merged.npz --SolsFileOut={checkname} --InterpMode={options["smoothingtype"]} --NCPU={options["NCPU_killms"]}',dryrun=dryrun)

        smoothoutname=f'{ddsols}_{start_time:.2f}_smoothed.npz'

        if InterpToMSListFreqs:
            interp_outname=f"{smoothoutname}_{start_time:.2f}_interp.npz"
            checkname=interp_outname
            if options['restart'] and os.path.isfile(checkname):
                warn('Solutions file '+checkname+' already exists, not running InterpSols step')
            else:
                command=f"InterpSols.py --SolsFileIn {smoothoutname} --SolsFileOut {interp_outname} --MSOutFreq {InterpToMSListFreqs} --NCPU={options['NCPU_killms']}"
                run(command,dryrun=dryrun)
        
        for i in range(0,len(full_sollist)):
            if start_times[i] == start_time:
                if not SkipSmooth:
                    symsolname = full_sollist[i].replace(ddsols,ddsols+'_smoothed')
                else:
                    symsolname = full_sollist[i].replace(ddsols,ddsols+'_merged')                 
                if not dryrun:
                    # always overwrite the symlink to allow the dataset to move -- costs nothing
                    if os.path.islink(symsolname):
                        warn('Symlink ' + symsolname + ' already exists, recreating')
                        os.unlink(symsolname)
                    if not SkipSmooth:
                        os.symlink(os.path.abspath(f'{ddsols}_{start_time:.2f}_smoothed.npz'),symsolname)
                    else:
                        os.symlink(os.path.abspath(f'{ddsols}_{start_time:.2f}_merged.npz'),symsolname)
                    
                    
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
    warn(f'nch_step={nch_step}')

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
    warn(f'Using ({nt_step},{nch_step}) for self-cal run of {msname} with (<|model|>,std)=({M:.2f},{S:.2f}) giving SNR={SNR:.2f}')
    
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
    report(f'Ingesting dynamic spectra ({obsid}) into the database')
    with SurveysDB() as sdb:
        sdb.cur.execute('lock table spectra write')
        field=os.path.basename(os.getcwd())
        g=glob.glob('*DynSpecs_'+obsid)
        for f in g:
            if '.tgz' in f:
                continue
            bits=f.split('_')
            obsid=bits[1]
            CatName=f+'/Catalog.npy'
            print(f"Loading {CatName}")
            try:
                catalogue=np.load(CatName)
            except:
                print(f"   {CatName} does not exist")
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
            sdb.cur.execute(f'delete from spectra where obsid="{obsid}"')
            for i,r in enumerate(catalogue):
                name=r['Name']
                if isinstance(name,np.bytes_):
                    name=name.decode('utf-8')
                sExec = 'insert into spectra values ( "%s", "%s", "%s", "%s", "%s", "%s", %.7f, %.7f, %g, %g, %g, %g )' % (field+'_'+obsid+'_'+str(i), name, r['Type'], field, obsid, fd[name], r['ra']*180.0/np.pi, r['dec']*180.0/np.pi, r['FluxI'], r['FluxV'], r['sigFluxI'], r['sigFluxV'])
                print(sExec)
                sdb.cur.execute(sExec)

def subtract_vis(mslist=None,colname_a="CORRECTED_DATA",colname_b="DATA_SUB",out_colname="DATA_SUB"):
    from pyrap.tables import table

    if mslist is not None:
        with open(mslist) as f:
            mslist = [msname.rstrip() for msname in f]
    else:
        raise ValueError("mslist is None")

    for msname in mslist:
        report(f'Subtracting: {out_colname} = {colname_a} - {colname_b}')
        t=table(msname,readonly=False)
        d=t.getcol(colname_a)
        p=t.getcol(colname_b)
        d-=p
        if out_colname not in t.colnames():
            report(f'Adding column {out_colname} in {msname}')
            desc=t.getcoldesc(colname_a)
            desc["name"]=out_colname
            desc['comment']=desc['comment'].replace(" ","_")
            t.addcols(desc)
        t.putcol(out_colname,d)
        t.close()

def subtractOuterSquare(o):
    #wide_imsize=o['wide_imsize']
    #NPixSmall=o['imsize'] #int(NPixLarge/float(o['fact_reduce_field']))
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

    ddf_image('image_full_wide', o['mslist'],
            cleanmode=o['cleanmode'],
            colname=colname,
            imsize=o['wide_imsize'],
            cellsize=o['wide_cell'],
            robust=o['wide_robust'],
            beamsize=o['wide_psf_arcsec'],
            phasecenter=o['phasecenter'],
            majorcycles=2,
            peakfactor=0.001,
            automask=True,
            automask_threshold=o['thresholds'][0],
            cleanmask=extmask,
            use_dicomodel=False,
            normalization=o['normalize'][2],
            uvrange=wide_uvrange,
            apply_weights=False,
            use_weightspectrum=o['use_weightspectrum'],
            AllowNegativeInitHMP=True,
            smooth=True,
            catcher=catcher,
            STEP=0,
    )


    external_mask='wide_external_mask.fits'
    make_external_mask(external_mask,'image_full_wide.dirty.fits',use_tgss=True,clobber=False)
    make_mask('image_full_wide.app.restored.fits',o['wide_threshold'],external_mask=external_mask,catcher=catcher)
    
    ddf_image('image_full_wide_im', o['mslist'],
            cleanmode=o['cleanmode'],
            colname=colname,
            imsize=o['wide_imsize'],
            cellsize=o['wide_cell'],
            robust=o['wide_robust'],
            beamsize=o['wide_psf_arcsec'],
            phasecenter=o['phasecenter'],
            majorcycles=1,
            peakfactor=0.001,
            automask=True,
            automask_threshold=o['thresholds'][0],
            cleanmask='image_full_wide.app.restored.fits.mask.fits',
            use_dicomodel=True,
            dicomodel_base='image_full_wide',
            reuse_psf=True,
            dirty_from_resid=True,
            normalization=o['normalize'][2],
            uvrange=wide_uvrange,
            apply_weights=False,
            use_weightspectrum=o['use_weightspectrum'],
            AllowNegativeInitHMP=True,
            smooth=True,
            catcher=catcher,
            STEP=0,
    )

    # predict outside the central rectangle
    NpixMaskSquare = np.floor(0.95*o['imsize']*o['cellsize']/o['wide_cell'])
    
    FileHasPredicted='image_full_wide_predict.HasPredicted'
    if o['restart'] and os.path.isfile(FileHasPredicted):
        warn(f'File {FileHasPredicted} already exists, skipping Predict step')
    else:
        ddf_image('image_full_wide_predict', o['full_mslist'],
                  cleanmode=o['cleanmode'],
                  colname=colname,
                  PredictSettings=('Predict', 'DATA_SUB', NpixMaskSquare),
                  imsize=o['wide_imsize'],
                  cellsize=o['wide_cell'],
                  robust=o['wide_robust'],
                  beamsize=o['wide_psf_arcsec'],
                  phasecenter=o['phasecenter'],
                  #majorcycles=1,
                  peakfactor=0.001,
                  #automask=True,
                  #automask_threshold=o['thresholds'][1],
                  #ddsols='wide_killms_p1',
                  #applysols='AP'
                  cleanmask='image_full_wide.app.restored.fits.mask.fits',
                  use_dicomodel=True,
                  dicomodel_base='image_full_wide_im',
                  #normalization=o['normalize'][0],
                  #uvrange=wide_uvrange,
                  apply_weights=False,
                  use_weightspectrum=o['use_weightspectrum'],
                  catcher=catcher,
                  STEP=0,
        )
        
        if not o['dryrun']:
            os.system(f"touch {FileHasPredicted}")


    # subtract predicted visibilities
    FileHasSubtracted='image_full_wide_predict.HasSubtracted'
    if o['restart'] and os.path.isfile(FileHasSubtracted):
        warn(f'File {FileHasSubtracted} already exists, skipping subtract vis step')
    elif not o['dryrun']:
        subtract_vis(mslist=o['full_mslist'],colname_a=colname,colname_b="DATA_SUB",out_colname="DATA_SUB")
        os.system(f"touch {FileHasSubtracted}")

    ## test subtracted...
    ## sanity check
    ddf_image('image_full_wide_im_sub', o['mslist'],
              cleanmode=o['cleanmode'],
              colname='DATA_SUB',
              imsize=o['wide_imsize'],
              cellsize=o['wide_cell'],
              robust=o['wide_robust'],
              beamsize=o['wide_psf_arcsec'],
              phasecenter=o['phasecenter'],
              majorcycles=1,
              peakfactor=0.001,
              automask=True,
              automask_threshold=o['thresholds'][0],
              cleanmask='image_full_wide.app.restored.fits.mask.fits',
              use_dicomodel=False,
              reuse_psf=True,
              dirty_from_resid=False,
              normalization=o['normalize'][2],
              uvrange=wide_uvrange,
              apply_weights=False,
              use_weightspectrum=o['use_weightspectrum'],
              AllowNegativeInitHMP=True,
              smooth=True,
              catcher=catcher,
              STEP=0,
    )


def main(o=None):
    if o is None and MyPickle is not None:
        o=MyPickle.Load("ddf-pipeline.last")

    lCat=[]
    if ((o['tgss'] is not None) and ('$$' in o['tgss'])):
        if "DDF_PIPELINE_CATALOGS" not in list(os.environ.keys()):
            die("You need to define the environment variable DDF_PIPELINE_CATALOGS where your catalogs are located")
        o["tgss"]=o["tgss"].replace("$$",os.environ["DDF_PIPELINE_CATALOGS"])
        
    if (o['catalogues'] is not None) and np.any(['$$' in l for l in o['catalogues']]):
        if "DDF_PIPELINE_CATALOGS" not in list(os.environ.keys()):
            die("You need to define the environment variable DDF_PIPELINE_CATALOGS where your catalogs are located")
        o["catalogues"]=[l.replace("$$",os.environ["DDF_PIPELINE_CATALOGS"]) for l in o["catalogues"]]
        
    lCat=[]
    if o['catalogues'] is not None:
        lCat+=o["catalogues"]
    if o['tgss'] is not None:
        lCat+=[o["tgss"]]
    if not o['dryrun']:
        for fCat in lCat:
            if not os.path.isfile(fCat):
                die(f"Catalog {fCat} does not exist")

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

    # In dryrun mode, create placeholder mslist files so the pipeline doesn't crash
    # They contain a single dummy MS name for command printing.
    _dryrun_tmpfiles = []
    if o['dryrun'] and o['dummy_msfile']:
        for key in ('mslist', 'full_mslist'):
            if o[key] is not None and not os.path.isfile(o[key]):
                with open(o[key], 'w') as f:
                    f.write('dummy.ms\n')
                _dryrun_tmpfiles.append(o[key])

    # Set column name for first steps
    colname=o['colname']

    # Check if the column exists in one MS. Important to do this
    # before we check imaging weights, because that will create empty
    # versions of e.g. CORRECTED_DATA
    if not o['dryrun']:
        mslist=[s.strip() for s in open(o['mslist']).readlines()]
        t = pt.table(mslist[0])
        try:
            dummy=t.getcoldesc(colname)
        except RuntimeError:
            dummy=None
        t.close()
        if dummy is None:
            die(f'Dataset does not contain the column {colname}')
    
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
    if not o['dryrun']:
        if o['full_mslist'] is not None:
            new=check_imaging_weight(o['full_mslist'])
        else:
            new=check_imaging_weight(o['mslist'])
    else:
        new=False

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
        if o['exitafter'] == 'wide':
            warn('User specified exit after wide field source subtraction')
            stop(2)
    
    # start of 'Big If' for reducing multiple datasets with a pre-made sky model
    if o['basedicomodel'] is None:
        # ##########################################################
        # Initial dirty image to allow an external (TGSS) mask to be made
        separator("Initial dirty")
        ddf_image('image_dirin_SSD_init', o['mslist'],
                  cleanmode=o['cleanmode'],
                  colname=colname,
                  robust=o['image_robust'],
                  phasecenter=o['phasecenter'],
                  majorcycles=0,
                  peakfactor=0.05,
                  automask_threshold=10.0,
                  cleanmask=None,
                  reuse_psf=False,
                  reuse_dirty=False,
                  uvrange=uvrange,
                  apply_weights=o['apply_weights'][0],
                  use_weightspectrum=o['use_weightspectrum'],
                  clusterfile=None,
                  catcher=catcher,
                  STEP=0,
        )


        separator("External mask")
        external_mask='external_mask.fits'
        make_external_mask(external_mask,'image_dirin_SSD_init.dirty.fits',use_tgss=True,clobber=False)
        if o['external_fits_mask'] is not None and not o['dryrun']:
            merge_mask(external_mask,o['external_fits_mask'],external_mask)

        # Deep SSD clean with this external mask and automasking
        separator("DI Deconv (externally defined sources)")
        CurrentBaseDicoModelName=ddf_image('image_dirin_SSD', o['mslist'],
                  cleanmode=o['cleanmode'],
                  colname=colname,
                  robust=o['image_robust'],
                  phasecenter=o['phasecenter'],
                  majorcycles=1,
                  peakfactor=0.01,
                  rms_factor=3,
                  automask=True,
                  automask_threshold=o['thresholds'][0],
                  cleanmask=external_mask,
                  reuse_psf=True,
                  reuse_dirty=True,
                  uvrange=uvrange,
                  apply_weights=o['apply_weights'][0],
                  use_weightspectrum=o['use_weightspectrum'],
                  clusterfile=None,
                  catcher=catcher,
                  STEP=0,
        )

    
        separator("Make the diffuse emission mask")
        # Make the diffuse emission mask
        _=make_mask('image_dirin_SSD.residual01.fits',
                    o['thres_outmaskextended'],
                    external_mask=external_mask,
                    catcher=catcher,
                    OutMaskExtended="MaskDiffuse")

        if o['use_maskdiffuse'] and not o['dryrun'] and external_mask is not None:
            separator("Merge diffuse emission mask into external mask")
            merge_mask(external_mask,"MaskDiffuse.fits",external_mask)
        
        # make a mask from the final image
        separator("Make mask for next iteration")
        CurrentMaskName=make_mask('image_dirin_SSD.app.restored.fits',
                                  o['thresholds'][0],
                                  external_mask=external_mask,
                                  catcher=catcher) if o['use_external_mask'][0] else None
        
        separator("Continue deconvolution")
        CurrentBaseDicoModelName=ddf_image('image_dirin_SSD_m', o['mslist'],
                  cleanmode=o['cleanmode'],
                  colname=colname,
                  PredictSettings=('Clean', 'DD_PREDICT'),
                  robust=o['image_robust'],
                  phasecenter=o['phasecenter'],
                  majorcycles=2,
                  peakfactor=0.001,
                  rms_factor=0,
                  automask=True,
                  automask_threshold=o['thresholds'][0],
                  cleanmask=CurrentMaskName,
                  use_dicomodel=True,
                  dicomodel_base=CurrentBaseDicoModelName,
                  reuse_psf=True,
                  dirty_from_resid=True,
                  uvrange=uvrange,
                  apply_weights=o['apply_weights'][0],
                  use_weightspectrum=o['use_weightspectrum'],
                  RMSFactorInitHMP=1.0,
                  MaxMinorIterInitHMP=10000,
                  clusterfile=None,
                  catcher=catcher,
                  STEP=0,
        )

        if o['exitafter'] == 'initial':
            warn('User specified exit after initial image')
            stop(2)


        #########################
        if o['clusterfile'] is None:
            separator("Cluster the sky model")
            ClusterFile='image_dirin_SSD_m.npy.ClusterCat.npy'
            if o['use_maskdiffuse']:
                clusterGA(imagename="image_dirin_SSD_m.app.restored.fits",
                      OutClusterCat=ClusterFile,
                      use_makemask_products=True)
            else:
                clusterGA(imagename="image_dirin_SSD_m.app.restored.fits",
                      OutClusterCat=ClusterFile,
                      use_makemask_products=False)
        else:
            ClusterFile=o['clusterfile']
            warn('Using user-specifed cluster file '+ClusterFile)

        #########################
        if o['clearcache'] or new or o['redofrom']:
            clearcache(o['mslist'],o)

        separator("Deconv clustered DI image")
        CurrentBaseDicoModelName=ddf_image('image_dirin_SSD_m_c', o['mslist'],
                  cleanmode=o['cleanmode'],
                  colname=colname,
                  PredictSettings=('Clean', 'DD_PREDICT'),
                  robust=o['image_robust'],
                  #reuse_psf=True,
                  #reuse_dirty=True,
                  phasecenter=o['phasecenter'],
                  majorcycles=1,
                  peakfactor=0.001,
                  rms_factor=0,
                  automask=True,
                  automask_threshold=o['thresholds'][0],
                  cleanmask=CurrentMaskName,
                  use_dicomodel=True,
                  #dirty_from_resid=True,
                  dicomodel_base=CurrentBaseDicoModelName,
                  uvrange=uvrange,
                  apply_weights=o['apply_weights'][0],
                  use_weightspectrum=o['use_weightspectrum'],
                  RMSFactorInitHMP=1.0,
                  MaxMinorIterInitHMP=10000,
                  clusterfile=ClusterFile,
                  catcher=catcher,
                  STEP=0,
        )

        if o['exitafter'] == 'dirin':
            warn('User specified exit after image_dirin.')
            stop(2)

        if not o['skip_di']:
            separator("DI CAL")
            ########################
            killms_data('PredictDI_0',o['mslist'],'DIS0',colname=colname,
                        dicomodel=f'{CurrentBaseDicoModelName}.DicoModel',
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
            _=ddf_image('image_dirin_SSD_m_c_di', o['mslist'],
                    cleanmode=o['cleanmode'],
                    colname=colname,
                    PredictSettings=('Clean', 'DD_PREDICT'),
                    robust=o['image_robust'],
                    #reuse_psf=True,
                    phasecenter=o['phasecenter'],
                    majorcycles=0,
                    peakfactor=0.001,
                    rms_factor=0,
                    automask=True,
                    automask_threshold=o['thresholds'][0],
                    cleanmask=CurrentMaskName,
                    use_dicomodel=True,
                    dicomodel_base=CurrentBaseDicoModelName,
                    #dirty_from_resid=True,
                    uvrange=uvrange,
                    #o['apply_weights'][0],
                    apply_weights=True,
                    RMSFactorInitHMP=1.0,
                    MaxMinorIterInitHMP=10000,
                    clusterfile=ClusterFile,
                    catcher=catcher,
                    STEP=0,
            )

            CurrentBaseDicoModelName=ddf_image('image_dirin_SSD_m_c_di_m', o['mslist'],
                      cleanmode=o['cleanmode'],
                      colname=colname,
                      PredictSettings=('Clean', 'DD_PREDICT'),
                      robust=o['image_robust'],
                      phasecenter=o['phasecenter'],
                      majorcycles=1,
                      peakfactor=0.001,
                      rms_factor=0,
                      automask=True,
                      automask_threshold=o['thresholds'][0],
                      cleanmask=CurrentMaskName,
                      use_dicomodel=True,
                      dicomodel_base=CurrentBaseDicoModelName,
                      #dirty_from_resid=True,
                      reuse_psf=True,
                      reuse_dirty=True,
                      uvrange=uvrange,
                      #o['apply_weights'][0],
                      apply_weights=True,
                      RMSFactorInitHMP=1.0,
                      MaxMinorIterInitHMP=10000,
                      clusterfile=ClusterFile,
                      catcher=catcher,
                      STEP=0,
            )

            # make a mask from the full-res image
            separator("Make mask for next iteration")
            CurrentMaskName=make_mask('image_dirin_SSD_m_c_di_m.app.restored.fits',
                                    o['thresholds'][1],
                                    external_mask=external_mask,
                                    catcher=catcher) if o['use_external_mask'][1] else None

            if o['exitafter'] == 'dirin_di':
                warn('User specified exit after image_dirin with DI calibration.')
                stop(2)


        separator("DD calibration")
        CurrentDDkMSSolName=killms_data(CurrentBaseDicoModelName,o['mslist'],'DDS0',colname=colname,
                                        dicomodel=f'{CurrentBaseDicoModelName}.DicoModel',
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
        CurrentBaseDicoModelName=ddf_image('image_phase1', o['mslist'],
                  cleanmode=o['cleanmode'],
                  colname=colname,
                  PredictSettings=('Clean', 'DD_PREDICT'),
                  robust=o['image_robust'],
                  phasecenter=o['phasecenter'],
                  majorcycles=2,
                  peakfactor=0.001,
                  automask=True,
                  automask_threshold=o['thresholds'][1],
                  cleanmask=CurrentMaskName,
                  ddsols=CurrentDDkMSSolName,
                  applysols=o['apply_sols'][0],
                  use_dicomodel=True,
                  dicomodel_base=CurrentBaseDicoModelName,
                  normalization=o['normalize'][0],
                  uvrange=uvrange,
                  apply_weights=o['apply_weights'][1],
                  use_weightspectrum=o['use_weightspectrum'],
                  RMSFactorInitHMP=1.0,
                  MaxMinorIterInitHMP=10000,
                  catcher=catcher,
                  STEP=1,
        )


        if o['exitafter'] == 'phase':
            warn('User specified exit after phase-only deconvolution.')
            stop(2)

        separator("Mask for deeper deconv")
        CurrentMaskName=make_mask('image_phase1.app.restored.fits',o['thresholds'][1],external_mask=external_mask,catcher=catcher) if o['use_external_mask'][1] else None
        CurrentBaseDicoModelName=mask_dicomodel('image_phase1.DicoModel',CurrentMaskName if o['use_external_mask'][1] else None,'image_phase1_masked.DicoModel',catcher=catcher)

        separator("DD calibration")
        CurrentDDkMSSolName=killms_data('image_phase1',o['mslist'],'DDS1',colname=colname,
                                        dicomodel=f'{CurrentBaseDicoModelName}.DicoModel',
                                        CovQ=0.02,
                                        clusterfile=ClusterFile,
                                        niterkf=o['NIterKF'][2],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                                        dt=o['dt_slow'],
                                        catcher=catcher,NChanSols=o['NChanSols'],
                                        EvolutionSolFile=CurrentDDkMSSolName,
                                        MergeSmooth=o['smoothing'])
        ##############################################

        separator("AmpPhase deconv")
        CurrentBaseDicoModelName=ddf_image('image_ampphase1', o['mslist'],
                  cleanmode=o['cleanmode'],
                  colname=colname,
                  PredictSettings=('Clean', 'DD_PREDICT'),
                  robust=o['image_robust'],
                  phasecenter=o['phasecenter'],
                  majorcycles=1,
                  peakfactor=0.001,
                  automask=True,
                  automask_threshold=o['thresholds'][1],
                  cleanmask=CurrentMaskName,
                  ddsols=CurrentDDkMSSolName,
                  applysols=o['apply_sols'][1],
                  use_dicomodel=True,
                  dicomodel_base=CurrentBaseDicoModelName,
                  normalization=o['normalize'][0],
                  uvrange=uvrange,
                  apply_weights=o['apply_weights'][1],
                  use_weightspectrum=o['use_weightspectrum'],
                  #AllowNegativeInitHMP=True,
                  RMSFactorInitHMP=1.0,
                  MaxMinorIterInitHMP=10000,
                  catcher=catcher,
                  STEP=1,
        )

        if o['exitafter'] == 'ampphase':
            warn('User specified exit after amp-phase deconvolution.')
            stop(2)

        separator("Make Mask")
        CurrentMaskName=make_mask('image_ampphase1.app.restored.fits',o['thresholds'][1],external_mask=external_mask,catcher=catcher) if o['use_external_mask'][1] else None
        CurrentBaseDicoModelName=mask_dicomodel('image_ampphase1.DicoModel',CurrentMaskName if o['use_external_mask'][1] else None,'image_ampphase1m_masked.DicoModel',catcher=catcher)

        if not o['skip_di']:
            separator("Second DI calibration")
            ddf_image('Predict_DI1', o['mslist'],
                      cleanmode=o['cleanmode'],
                      colname=colname,
                      PredictSettings=('Predict', 'DD_PREDICT'),
                      robust=o['image_robust'],
                      phasecenter=o['phasecenter'],
                      majorcycles=1,
                      peakfactor=0.001,
                      automask=True,
                      automask_threshold=o['thresholds'][1],
                      cleanmask=CurrentMaskName,
                      ddsols=CurrentDDkMSSolName,
                      applysols=o['apply_sols'][2],
                      use_dicomodel=True,
                      dicomodel_base=CurrentBaseDicoModelName,
                      normalization=o['normalize'][0],
                      uvrange=uvrange,
                      apply_weights=o['apply_weights'][1],
                      use_weightspectrum=o['use_weightspectrum'],
                      RMSFactorInitHMP=1.0,
                      MaxMinorIterInitHMP=10000,
                      catcher=catcher,
                      STEP=1,
            )
            
            separator("Another DI step")
            if o['bootstrap']:
                colname='SCALED_DATA'
            elif o['do_wide']:
                colname='DATA_SUB'
            else:
                colname=o['colname']

            killms_data('PredictDI_1',o['mslist'],'DIS1',colname=colname,
                        dicomodel=f'{CurrentBaseDicoModelName}.DicoModel',
                        #clusterfile=ClusterFile,
                        niterkf=o['NIterKF'][3],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                        catcher=catcher,
                        dt=o['dt_di'],
                        DISettings=("CohJones","IFull","DD_PREDICT","DATA_DI_CORRECTED"),UpdateWeights=0)
            # cubical_data(o['mslist'],
            #              NameSol="DIS1",
            #              n_dt=1,
            #              n_df=2,
            #              n_DT=None,
            #              DataColName=o['colname'],
            #              ModelColName="DD_PREDICT",
            #              OutColName="DATA_DI_CORRECTED")

            colname='DATA_DI_CORRECTED' # again
            CurrentBaseDicoModelName=ddf_image('image_ampphase1_di', o['mslist'],
                      cleanmode=o['cleanmode'],
                      colname=colname,
                      PredictSettings=('Clean', 'DD_PREDICT'),
                      robust=o['image_robust'],
                      phasecenter=o['phasecenter'],
                      majorcycles=1,
                      peakfactor=0.001,
                      automask=True,
                      automask_threshold=o['thresholds'][1],
                      cleanmask=CurrentMaskName,
                      ddsols=CurrentDDkMSSolName,
                      applysols=o['apply_sols'][3],
                      use_dicomodel=True,
                      dicomodel_base=CurrentBaseDicoModelName,
                      normalization=o['normalize'][0],
                      uvrange=uvrange,
                      apply_weights=o['apply_weights'][1],
                      use_weightspectrum=o['use_weightspectrum'],
                      #AllowNegativeInitHMP=True,
                      RMSFactorInitHMP=1.0,
                      MaxMinorIterInitHMP=10000,
                      catcher=catcher,
                      STEP=1,
            )

            if o['exitafter'] == 'ampphase_di':
                warn('User specified exit after amp-phase plus DI deconvolution.')
                stop(2)


        # small mslist cache not needed from this point so clear it to
        # save disk space
        if o['clearcache_end']:
            clearcache(o['mslist'],o)

        if o['full_mslist'] is None:
            warn('No full mslist provided, stopping here')
            summary(o)
            stop(3)

        # #########################################################################
        # ###############                  BIG MSLIST               ###############
        # #########################################################################

        if o['bootstrap']:
            colname='SCALED_DATA'
        else:
            if o['do_wide']:
                colname='DATA_SUB'
            else:
                colname=o['colname']

        if not o['skip_di']:
            separator("Make Mask")
            CurrentMaskName=make_mask('image_ampphase1_di.app.restored.fits',o['thresholds'][1],external_mask=external_mask,catcher=catcher) if o['use_external_mask'][1] else None
            CurrentBaseDicoModelName=mask_dicomodel('image_ampphase1_di.DicoModel',CurrentMaskName if o['use_external_mask'][1] else None,'image_ampphase1_di_masked.DicoModel',catcher=catcher)
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
        CurrentMaskName = o['basemaskname'] if o['use_external_mask'][2] else None
        CurrentBaseDicoModelName = o['basedicomodel']
        CurrentImageName = o['baseimagename']
        external_mask = CurrentMaskName

    separator("DD calibration of full mslist")
    CurrentDDkMSSolName=killms_data(CurrentImageName,o['full_mslist'],'DDS2_full',
                                    colname=colname,
                                    dicomodel=f'{CurrentBaseDicoModelName}.DicoModel',
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
        separator("Compute DD Predict (full mslist)")
        ddf_image('Predict_DDS2', o['full_mslist'],
                  cleanmode=o['cleanmode'],
                  colname=colname,
                  PredictSettings=('Predict', 'DD_PREDICT'),
                  robust=o['image_robust'],
                  phasecenter=o['phasecenter'],
                  majorcycles=1,
                  peakfactor=0.01,
                  automask=True,
                  automask_threshold=o['thresholds'][1],
                  ddsols=CurrentDDkMSSolName,
                  applysols=o['apply_sols'][4],
                  use_dicomodel=True,
                  dicomodel_base=CurrentBaseDicoModelName,
                  normalization=o['normalize'][0],
                  uvrange=uvrange,
                  apply_weights=o['apply_weights'][0],
                  use_weightspectrum=o['use_weightspectrum'],
                  catcher=catcher,
                  STEP=1,
        )

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
                    dicomodel=f'{CurrentBaseDicoModelName}.DicoModel',
                    clusterfile=ClusterFile,
                    niterkf=o['NIterKF'][5],uvrange=killms_uvrange,wtuv=o['wtuv'],robust=o['solutions_robust'],
                    catcher=catcher,
                    dt=o['dt_di'],
                    DISettings=("CohJones","IFull","DD_PREDICT","DATA_DI_CORRECTED"),UpdateWeights=0)
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
    
    ddf_image(ImageName, o['full_mslist'],
            cleanmode=o['cleanmode'],
            colname=colname,
            robust=o['final_robust'],
            phasecenter=o['phasecenter'],
            majorcycles=0,
            peakfactor=0.001,
            automask=True,
            automask_threshold=o['thresholds'][2],
            cleanmask=CurrentMaskName,
            ddsols=CurrentDDkMSSolName,
            applysols=o['apply_sols'][5],
            use_dicomodel=True,
            dicomodel_base=CurrentBaseDicoModelName,
            normalization=o['normalize'][1],
            uvrange=uvrange,
            apply_weights=o['apply_weights'][2],
            use_weightspectrum=o['use_weightspectrum'],
            AllowNegativeInitHMP=True,
            smooth=True,
            catcher=catcher,
            STEP=2,
            **ddf_kw
    )

    if o['exitafter'] == 'fullampphase':
        warn('User specified exit after image_ampphase.')
        stop(2)
        
    separator("MakeMask")
    CurrentMaskName=make_mask(ImageName+'.app.restored.fits',10,external_mask=external_mask,catcher=catcher) if o['use_external_mask'][2] else None

    separator("Finish Deconvolution AP (full mslist)")
    if not o['skip_di']:
        ImageName = 'image_full_ampphase_di_m'
    else:
        ImageName = 'image_full_ampphase_m'

    CurrentBaseDicoModelName=ddf_image(ImageName, o['full_mslist'],
              cleanmode=o['cleanmode'],
              colname=colname,
              robust=o['final_robust'],
              phasecenter=o['phasecenter'],
              majorcycles=1,
              peakfactor=0.001,
              automask=True,
              automask_threshold=o['thresholds'][2],
              cleanmask=CurrentMaskName,
              ddsols=CurrentDDkMSSolName,
              applysols=o['apply_sols'][5],
              use_dicomodel=True,
              dicomodel_base=CurrentBaseDicoModelName,
              reuse_psf=True,
              reuse_dirty=True,
              normalization=o['normalize'][1],
              uvrange=uvrange,
              apply_weights=o['apply_weights'][2],
              use_weightspectrum=o['use_weightspectrum'],
              AllowNegativeInitHMP=True,
              RMSFactorInitHMP=0.5,
              MaxMinorIterInitHMP=10000,
              smooth=True,
              catcher=catcher,
              STEP=2,
              **ddf_kw
    )

    separator("MakeMask")
    CurrentMaskName=make_mask(ImageName+'.app.restored.fits',o['thresholds'][2],external_mask=external_mask,catcher=catcher) if o['use_external_mask'][2] else None
    CurrentBaseDicoModelName=mask_dicomodel(ImageName+'.DicoModel',CurrentMaskName if o['use_external_mask'][2] else None,ImageName+'_masked.DicoModel',catcher=catcher)
            
    separator("DD Calibration (full mslist)")
    CurrentDDkMSSolName=killms_data(ImageName,
                                    o['full_mslist'],'DDS3_full',
                                    colname=colname,
                                    clusterfile=ClusterFile,
                                    dicomodel=f'{CurrentBaseDicoModelName}.DicoModel',
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
                                        dicomodel=f'{CurrentBaseDicoModelName}.DicoModel',
                                        uvrange=[o['uvmin_very_slow'],1000.],
                                        wtuv=o['wtuv'],
                                        robust=o['solutions_robust'],
                                        SkipSmooth=True,MergeSmooth=True,
                                        SigmaFilterOutliers=o['sigma_clip'],
                                        dt=o['dt_very_slow'],catcher=catcher,
                                        PreApplySols=CurrentDDkMSSolName_FastSmoothed)#,EvolutionSolFile=CurrentDDkMSSolName)

        CurrentDDkMSSolName=f"[{CurrentDDkMSSolName_FastSmoothed},{CurrentDDkMSSolName}]"

    if o['low_psf_arcsec'] is not None: # low-res image requested
        separator("Low-resolution image")

        ddf_kw={}
        ddf_kw['beamsize']=o['low_psf_arcsec']
        if o['low_psf_minor_arcsec'] is not None:
            if o['low_psf_pa_deg'] is None:
                die('If low minor axis is supplied, position angle must be supplied too')
            ddf_kw['beamsize_minor']=o['low_psf_minor_arcsec']
            ddf_kw['beamsize_pa']=o['low_psf_pa_deg']

        low_uvrange=[o['image_uvmin'],2.5*206.0/o['low_psf_arcsec']]
        if o['low_imsize'] is not None:
            low_imsize=o['low_imsize'] # allow over-ride
        else:
            low_imsize=o['imsize']*o['cellsize']/o['low_cell']
            # if mask-low exists then use it
        
        extmask=None
        if os.path.isfile('bootstrap-mask-low.fits'):
            extmask='bootstrap-mask-low.fits'
            # can be empty, in which case recent versions of DDF throw
            # an error, so check and drop it if it is
            hdu=fits.open(extmask)
            if not np.any(hdu[0].data>0):
                warn('Bootstrap external mask is blank, using only internal masking')
                extmask=None
            hdu.close()

        ddf_image('image_full_low', o['full_mslist'],
                cleanmode=o['cleanmode'],
                colname=colname,
                imsize=low_imsize,
                cellsize=o['low_cell'],
                robust=o['low_robust'],
                phasecenter=o['phasecenter'],
                majorcycles=2,
                peakfactor=0.001,
                automask=True,
                automask_threshold=5,
                cleanmask=extmask,
                ddsols=CurrentDDkMSSolName,
                applysols=o['apply_sols'][6],
                use_dicomodel=False,
                normalization=o['normalize'][2],
                uvrange=low_uvrange,
                AllowNegativeInitHMP=True,
                smooth=True,
                catcher=catcher,
                STEP=3,
                **ddf_kw,
        )

        make_mask('image_full_low.app.restored.fits',o['low_threshold'],external_mask=extmask,catcher=catcher)

        ddf_image('image_full_low_im', o['full_mslist'],
                  cleanmode=o['cleanmode'],
                  colname=colname,
                  imsize=low_imsize,
                  cellsize=o['low_cell'],
                  robust=o['low_robust'],
                  phasecenter=o['phasecenter'],
                  majorcycles=1,
                  peakfactor=0.001,
                  automask=True,
                  automask_threshold=5,
                  cleanmask='image_full_low.app.restored.fits.mask.fits',
                  ddsols=CurrentDDkMSSolName,
                  applysols=o['apply_sols'][6],
                  use_dicomodel=True,
                  dicomodel_base='image_full_low',
                  reuse_psf=True,
                  dirty_from_resid=True,
                  normalization=o['normalize'][2],
                  uvrange=low_uvrange,
                  AllowNegativeInitHMP=True,
                  smooth=True,
                  catcher=catcher,
                  STEP=3,
                  **ddf_kw,
        )

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
                if not o['dryrun']:
                    die('Could not find the required products for the full-bw extended source mask!')
            report('Make_extended_mask returns')

        extmask='full-mask-low.fits'
        make_mask('image_full_low_im.app.restored.fits',o['low_threshold'],external_mask=extmask,catcher=catcher)

        ddf_image('image_full_low_m', o['full_mslist'],
                cleanmode=o['cleanmode'],
                colname=colname,
                imsize=low_imsize,
                cellsize=o['low_cell'],
                robust=o['low_robust'],
                beamsize=o['low_psf_arcsec'],
                phasecenter=o['phasecenter'],
                majorcycles=1,
                peakfactor=0.001,
                rms_factor=o['final_rmsfactor'],
                automask=True,
                automask_threshold=4,
                cleanmask='image_full_low_im.app.restored.fits.mask.fits',
                ddsols=CurrentDDkMSSolName,
                applysols=o['apply_sols'][6],
                use_dicomodel=True,
                dicomodel_base='image_full_low_im',
                reuse_psf=True,
                dirty_from_resid=True,
                normalization=o['normalize'][2],
                uvrange=low_uvrange,
                AllowNegativeInitHMP=True,
                smooth=True,
                catcher=catcher,
                STEP=3,
        )

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
        if not o['dryrun'] and download_required(o['method']):
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

    ddf_image(ImageName, o['full_mslist'],
            cleanmode=o['cleanmode'],
            colname=colname,
            PredictSettings=('Clean', 'DD_PREDICT'),
            robust=o['final_robust'],
            phasecenter=o['phasecenter'],
            majorcycles=1,
            peakfactor=0.001,
            automask=True,
            automask_threshold=o['thresholds'][2],
            cleanmask=CurrentMaskName,
            ddsols=CurrentDDkMSSolName,
            applysols=o['apply_sols'][6],
            use_dicomodel=True,
            dicomodel_base=CurrentBaseDicoModelName,
            reuse_psf=False,
            normalization=o['normalize'][1],
            uvrange=uvrange,
            apply_weights=o['apply_weights'][2],
            use_weightspectrum=o['use_weightspectrum'],
            AllowNegativeInitHMP=True,
            RMSFactorInitHMP=1.0,
            smooth=True,
            catcher=catcher,
            STEP=3,
            **ddf_kw
    )


    # check for the offset files
    if o['method'] is not None:
        separator('Offset correction')
        # have we got the catalogue?
        if not o['dryrun']:
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
            run('offsets.py '+' '.join(sys.argv[1:]),log=None,dryrun=o['dryrun'])

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
        if o['clearcache_end']:
            full_clearcache(o)
        from do_polcubes import do_polcubes
        separator('Stokes Q and U cubes')
        cthreads=[]
        flist=[]
        pol_mslists=[]

        if o['split_polcubes']:
            cubefiles=['image_full_low_StokesQ.cube.dirty.fits','image_full_low_StokesQ.cube.dirty.corr.fits','image_full_low_StokesU.cube.dirty.fits','image_full_low_StokesU.cube.dirty.corr.fits']
        else:
            cubefiles=['image_full_low_QU.cube.dirty.fits','image_full_low_QU.cube.dirty.corr.fits']
        if o['restart'] and os.path.isfile(cubefiles[0]+'.fz') and os.path.isfile(cubefiles[1]+'.fz'):
            warn('Compressed low QU cube product exists, not making new images')
        else:
            pol_mslists=do_polcubes(colname,CurrentDDkMSSolName,low_uvrange,'image_full_low',o['full_mslist'],ddf_kw,beamsize=o['low_psf_arcsec'],imsize=low_imsize,cellsize=o['low_cell'],robust=o['low_robust'],options=o,catcher=catcher)
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
            do_polcubes(colname,CurrentDDkMSSolName,vlow_uvrange,'image_full_vlow',o['full_mslist'],ddf_kw,beamsize=o['vlow_psf_arcsec'],imsize=o['vlow_imsize'],cellsize=o['vlow_cell'],robust=o['vlow_robust'],options=o,catcher=catcher)
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
        
    mslist_file = o['full_mslist'] or o['mslist']
    if mslist_file is not None and os.path.isfile(mslist_file):
        try:
            m=MSList(mslist_file)
            uobsid = set(m.obsids)
        except RuntimeError:
            if not o['dryrun']:
                raise
            m=None
            uobsid = set()
    else:
        m=None
        uobsid = set()
    stokesv_mslists=[]
    for obsid in uobsid:
        umslist = f'mslist-{obsid}.txt'
        stokesv_mslists.append(umslist)
        if not o['dryrun']:
            print('Writing ms list for obsids',umslist)
            with open(umslist,'w') as file:
                for ms,ob in zip(m.mss,m.obsids):
                    if ob==obsid:
                        file.write(ms+'\n')
    if o['stokesv']:
        for obsid in uobsid:
            separator(f'Stokes V image for {obsid}')
            ddf_image(f'image_full_high_stokesV_{obsid}', f'mslist-{obsid}.txt',
                      cleanmode=o['cleanmode'],
                      colname=colname,
                      cellsize=o['cellsize'],
                      robust=o['final_robust'],
                      phasecenter=o['phasecenter'],
                      stokes='IV',
                      majorcycles=0,
                      peakfactor=0.001,
                      automask=True,
                      automask_threshold=5,
                      ddsols=CurrentDDkMSSolName,
                      applysols=o['apply_sols'][6],
                      use_dicomodel=False,
                      normalization=o['normalize'][2],
                      uvrange=uvrange,
                      AllowNegativeInitHMP=True,
                      smooth=True,
                      catcher=catcher,
                      STEP=3,
                      **ddf_kw
            )

    if o['polcubes'] and o['compress_polcubes']:
        # cthreads and flist exist
        for thread in cthreads:
            if thread.is_alive():
                warn('Waiting for a compression thread to finish')
                thread.join()
        if o['delete_compressed']:
            for f in flist:
                if os.path.isfile(f+'.fz'):
                    warn(f'Deleting compressed file {f}')
                    os.remove(f)
                else:
                    if not o['dryrun']:
                        die('compressed files do not exist, compression must have failed')

    if o['do_dynspec']:
        separator('Dynamic spectra')

        if o['bright_threshold'] is not None and o['method'] is not None:
            warn('Finding bright sources from offsets list')
            from find_bright_offset_sources import find_bright
            bright_exists=find_bright(cutoff=o['bright_threshold'])
        LastImage="image_full_ampphase_di_m.NS.int.restored.fits"

        for obsid in uobsid:
            LastImageV=f"image_full_high_stokesV_{obsid}.dirty.corr.fits"
            warn(f'Running ms2dynspec for obsid {obsid}')
            umslist=f'mslist-{obsid}.txt'
            g=glob.glob(f'DynSpec*{obsid}*')
            if len(g)>0:
                warn(f'DynSpecs results directory {g[0]} already exists, skipping DynSpecs')
            else:
                DicoFacetName=f"{LastImage.split('.int.restored.fits')[0]}.DicoFacet"
                runcommand = f"ms2dynspec.py --ms {umslist} --data {colname} --model DD_PREDICT --sols {CurrentDDkMSSolName} --rad 2. --imageI {LastImage} --imageV {LastImageV} --LogBoring {o['nobar']} --SolsDir {o['SolsDir']} --BeamModel LOFAR --BeamNBand 1 --DicoFacet {DicoFacetName} --noff 100 --nMinOffPerFacet 5 --CutGainsMinMax 0.1,1.5 --SplitNonContiguous 1 --SavePDF 1 --FitsCatalog ${{DDF_PIPELINE_CATALOGS}}/dyn_spec_catalogue_addedexo_addvlotss.fits"
                
                if o['bright_threshold'] is not None and os.path.isfile('brightlist.csv'):
                    runcommand+=' --srclist brightlist.csv'
                run(runcommand,dryrun=o['dryrun'],log=logfilename('ms2dynspec.log'),quiet=o['quiet'])
                if use_database():
                    ingest_dynspec(obsid)


    if o['compress_ms'] and not o['dryrun']:
        separator(f'Compressing MS for archive -- column {colname}')
        os.system(f'archivems.sh . {colname}')
                
    separator('Write summary and tidy up')
    summary(o)

    # Clear caches if option set
    if o['clearcache_end']:
        extras=[]
        if spectral_mslist is not None:
            extras+=spectral_mslist
        if o['polcubes']:
            extras+=pol_mslists
        if o['stokesv']:
            extras+=stokesv_mslists
        full_clearcache(o,extras=extras)
    
    if use_database():
        update_status(None,'Complete',time='end_date',av=4)

    # Remove temporary msfiles made for dryrun mode
    for f in _dryrun_tmpfiles:
        os.remove(f)

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
