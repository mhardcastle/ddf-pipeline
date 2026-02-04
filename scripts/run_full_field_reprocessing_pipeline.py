#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from surveys_db import SurveysDB, tag_field, use_database, get_cluster
from parset import option_list
from options import options,print_options
import sys
import os
import glob
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import time
from subprocess import call
import subprocess
from reprocessing_utils import prepare_field,do_rclone_reproc_tape_upload,convert_summary_cfg,do_rclone_reproc_tape_download
import argparse
import threading
from auxcodes import run,warn,report, MSList
import numpy as np
import datetime
from pipeline import ingest_dynspec,ddf_image
from rclone import RClone
from getcpus import getcpus
import pyrap.tables as pt
import shutil
    
def check_cube_format(header):
    try:
        assert header["CTYPE1"].startswith("RA")
        assert header["CTYPE2"].startswith("DEC")
        assert header["CTYPE3"].startswith("STOKES")
        assert header["CTYPE4"].startswith("FREQ")
    except AssertionError:
        raise ValueError("Input cube must be in order: RA,DEC,STOKES,FREQ")

def redo_cube_headers(incubename,outcubename,stokesparam):
    hdu = fits.open(incubename)[0]
    header = hdu.header
    data = hdu.data
    check_cube_format(header)
    print('data shape:', hdu.data.shape)
    old_crval4 = header["CRVAL4"]
    # Survey individual Q or U cubes will be ordered: [Freqs, Dec, RA]
    data = data[:,:,:]
    header_axes_attributes = ["NAXIS", "CTYPE",  "CRVAL", "CRPIX", "CDELT", "CUNIT", "CROTA"]
    for attr in header_axes_attributes:
        attr += "4"
        if attr in header:
            del header[attr]
    header["CTYPE3"] = "freq"
    header["CUNIT3"] = "Hz"
    header["CRPIX3"] = 1
    header["CRVAL3"] = old_crval4
    header["CDELT3"] = 97656.25
    header["NAXIS"] = 4
    header["CTYPE4"] = "STOKES"
    header["CUNIT4"] = " "
    header["CRPIX4"] = 1.
    if stokesparam == 'Q':
        header["CRVAL4"] = 2
    if stokesparam == 'U':
        header["CRVAL4"] = 3
    header["CDELT4"] = 1.
    cube_shape = (1,hdu.data.shape[0], hdu.data.shape[1], hdu.data.shape[2])
    fits.writeto(outcubename, data.reshape(cube_shape), header, overwrite=True)
    del data # to clear from memory

def do_run_subtract(field):

    # Run subtract code
    executionstr = 'sub-sources-outside-region.py --timeavg=1 --boxfile=fullfield'
    print(executionstr)
    result=os.system(executionstr)
    if result!=0:
        raise RuntimeError('sub-sources-outside-region.py failed with error code %i' % result)

def do_run_dynspec(field):
    # Run subtract code
    try:
        DoRunDDF=("DDFacet ended successfully" not in open("image_full_ampphase_di_m.NS_SUB.log","r").readlines()[-1])
    except:
        DoRunDDF=True
    if DoRunDDF:
        executionstr = 'sub-sources-outside-region.py --timeavg=1 --freqavg=1 --boxfile=fullfield --onlyPredict --AlsoMakeResidualImage'
        print(executionstr)
        result=os.system(executionstr)
        if result!=0:
            raise RuntimeError('sub-sources-outside-region.py failed with error code %i' % result)
    else:
        print("DDFacet has already been successfully run, skipping")

    # executionstr = 'ms2dynspec.py --ms=big-mslist.txt --data DATA_SUB --model '

    ListMSName=[l.strip() for l in open("big-mslist.txt","r").readlines()]
    ListObsName=sorted(list(set([MSName.split("/")[-1].split("_")[0] for MSName in ListMSName])))
    AllOutputExist=True
    for ObsID in ListObsName:
        tgzName="DynSpecs_%s.tgz"%ObsID
        if not os.path.isfile(tgzName):
            AllOutputExist=False
            print("DynSpecMS output %s does not exist"%tgzName)
        else:
            print("DynSpecMS output %s exists"%tgzName)
            
    if not AllOutputExist:
        executionstr = 'ms2dynspec.py --tolerance 0.3 --ms big-mslist.txt --data DATA --model PREDICT_SUB --sols [DDS3_full_smoothed,DDS3_full_slow] --rad 2. --SolsDir SOLSDIR --BeamModel LOFAR --BeamNBand 1 --DicoFacet image_full_ampphase_di_m.NS_SUB.DicoFacet --noff 100 --nMinOffPerFacet 5 --CutGainsMinMax 0.1,1.5 --SplitNonContiguous 1 --imageI image_full_ampphase_di_m.NS.int.restored.fits --imageV image_full_high_stokesV.dirty.corr.fits --SavePDF 1 --FitsCatalog ${DDF_PIPELINE_CATALOGS}/dyn_spec_catalogue_addedexo_addvlotss.fits'
        print(executionstr)
        result=os.system(executionstr)
        if result!=0:
            raise RuntimeError('ms2dynspec.py failed with error code %i' % result)
    else:
        print("All DynSpecMS output exists, skipping... ")    

def transient_image(msfilename,imagename,galactic=False,options=None):
    t = pt.table(msfilename,readonly=True)
    numtimes = len(np.unique(t.getcol('TIME')))
    print(numtimes,msfilename)

    # Check if (compressed!) output files already exist
    if galactic:
        chanout = 16
        template='%s*-MFS-dirty.fits.fz'
    else:
        chanout = 1
        template='%s*-dirty.fits.fz'

    existingfiles = glob.glob(template%imagename)
    if len(existingfiles) == numtimes:
        print('All output files already exist, skipping making images again')
    else:
        # Create 8s images and 2 min images (assuming time resolution is 8s...). 
        if chanout > 1:
            os.system('wsclean -make-psf -intervals-out %s -padding 1.6 -interval 0 %s -auto-threshold 5 -channels-out %s -deconvolution-channels 3 -fit-spectral-pol 3 -scale 6asec -size 2200 2200 -join-channels -minuv-l 50 -maxuv-l 5000 -mgain 0.8 -multiscale -niter 0 -nmiter 0 -no-update-model-required -parallel-deconvolution 1500 -no-reorder -taper-gaussian 40asec -use-wgridder -weight briggs -0.25 -name %s_8sec %s'%(numtimes,numtimes,chanout,imagename,msfilename))
            os.system('wsclean -make-psf -intervals-out %s -padding 1.6 -interval 0 %s -auto-threshold 5 -channels-out %s -deconvolution-channels 3 -fit-spectral-pol 3 -scale 6asec -size 2200 2200 -join-channels -minuv-l 50 -maxuv-l 5000 -mgain 0.8 -multiscale -niter 0 -nmiter 0 -no-update-model-required -parallel-deconvolution 1500 -no-reorder -taper-gaussian 40asec -use-wgridder -weight briggs -0.25 -name %s_2min %s'%(int(numtimes/15),numtimes,chanout,imagename,msfilename))
        else:
            os.system('wsclean -make-psf -intervals-out %s -padding 1.6 -interval 0 %s -auto-threshold 5 -channels-out %s -scale 6asec -size 2200 2200 -minuv-l 50 -maxuv-l 5000 -mgain 0.8 -multiscale -niter 0 -nmiter 0 -no-update-model-required -parallel-deconvolution 1500 -no-reorder -taper-gaussian 40asec -use-wgridder -weight briggs -0.25 -name %s_8sec %s'%(numtimes,numtimes,chanout,imagename,msfilename))
            os.system('wsclean -make-psf -intervals-out %s -padding 1.6 -interval 0 %s -auto-threshold 5 -channels-out %s -scale 6asec -size 2200 2200 -minuv-l 50 -maxuv-l 5000 -mgain 0.8 -multiscale -niter 0 -nmiter 0 -no-update-model-required -parallel-deconvolution 1500 -no-reorder -taper-gaussian 40asec -use-wgridder -weight briggs -0.25 -name %s_2min %s'%(int(numtimes/15),numtimes,chanout,imagename,msfilename))

    outimages = glob.glob('%s*dirty.fits'%imagename)
    for outimage in outimages:
        compress_fits(outimage,o['fpack_q'])
        os.remove(outimage)
    unwanted=glob.glob(imagename+'*.fits')
    for outimage in unwanted:
        os.remove(outimage)

def image_vlow(ncpu,wd=None,cache='.'):
    if wd is not None:
        os.chdir(wd)
    run('CleanSHM.py')
    if not os.path.exists('image_full_vlow_nocut.app.restored.fits'):
        run('DDF.py --Output-Name=image_full_vlow_nocut --Data-MS=big-mslist.txt --Deconv-PeakFactor 0.001000 --Data-ColName DATA --Parallel-NCPU=%i --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=2 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust -0.20000 --Image-NPix=2000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 15.00000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=%s --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --GAClean-AllowNegativeInitHMP True --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 60.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=4.00 --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=[DDS3_full_smoothed,DDS3_full_slow] --Selection-UVRangeKm=[0.000000,7.0] --GAClean-MinSizeInit=10 --Beam-Smooth=1 --Debug-Pdb=never --Cache-DirWisdomFFTW=%s' % (ncpu,cache,cache))
    vlowmask='image_full_vlow_nocut.app.restored.fits.mask.fits'
    if not os.path.exists(vlowmask):
        vlowmask = make_mask('image_full_vlow_nocut.app.restored.fits',3.0)
        
    if not os.path.exists('image_full_vlow_nocut_m.app.restored.fits'):
        run('DDF.py --Output-Name=image_full_vlow_nocut_m --Data-MS=big-mslist.txt --Deconv-PeakFactor 0.001000 --Data-ColName DATA --Parallel-NCPU=%i --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=2 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust -0.20000 --Image-NPix=2000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 15.00000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=%s --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --GAClean-AllowNegativeInitHMP True --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 60.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=3.00 --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=[DDS3_full_smoothed,DDS3_full_slow] --Selection-UVRangeKm=[0.000000,7.0] --GAClean-MinSizeInit=10 --Beam-Smooth=1 --Debug-Pdb=never --Cache-DirWisdomFFTW=%s --Predict-InitDicoModel=image_full_vlow_nocut.DicoModel --Mask-External=%s'% (ncpu,cache,cache,vlowmask))


def image_vlow_sub(useIDG=False,stringMultiscaleScales = "0,4,8,16,32,64,128,256", stringMultiscaleScalesIDG = "0,4,8,16,32,64",
                sizePixels=2000,taperGaussian=60.0,beamsize=60.0,scale=15.0,numberOfSubbands=6,uvDistanceMinLambda=0,IDGMode="cpu",name='WSCLEAN_low'):

    os.system('rm -rf *.archive?')
    print('Subtracting at low resolution')
    run('sub-sources-outside-region.py -b fullfield -p SUB -c DATA --uselowres -t 8 -f 5')
    g=glob.glob('*.archive?')
    stringMSs=' '.join(g)
    print("Imaging with multiscale CLEAN (at low-resolution)...")
    if (useIDG):
        # Generate IDG configuration file.
        print("Generating IDG configuration file...")
        with open("aconfig.txt", "w") as configurationFile:
            configurationFile.write("aterms=[beam]\nbeam.differential = true\nbeam.update_interval = 600\nbeam.usechannelfreq = true")

        command = "wsclean -no-update-model-required -size " + str(sizePixels) + " " + str(sizePixels) + " -reorder -weight briggs -0.5 -weighting-rank-filter 3 -clean-border 1 -mgain 0.8 -no-fit-beam -data-column DATA -join-channels -channels-out " + str(numberOfSubbands) + " -padding 1.2 -multiscale -multiscale-scales " + stringMultiscaleScalesIDG + " -auto-mask 3.0 -auto-threshold 2.5 -taper-gaussian " + str(taperGaussian) + "arcsec -circular-beam -beam-size " + str(beamsize) + "arcsec -pol i -name " + name + " -scale " + str(scale) + "arcsec -niter 100000 -minuv-l " + str(uvDistanceMinLambda) + " -use-idg -idg-mode " + IDGMode + " -aterm-kernel-size 16 -aterm-config aconfig.txt " + stringMSs
    else:
        command = "wsclean -no-update-model-required -size " + str(sizePixels) + " " + str(sizePixels) + " -reorder -weight briggs -0.5 -weighting-rank-filter 3 -clean-border 1 -mgain 0.8 -no-fit-beam -data-column DATA -join-channels -channels-out " + str(numberOfSubbands) + " -padding 1.2 -multiscale -multiscale-scales " + stringMultiscaleScales + " -auto-mask 3.0 -auto-threshold 2.5 -taper-gaussian " + str(taperGaussian) + "arcsec -circular-beam -beam-size " + str(beamsize) + "arcsec -pol i -name " + name + " -scale " + str(scale) + "arcsec -niter 100000 -minuv-l " + str(uvDistanceMinLambda) + " -baseline-averaging 10.0 " + stringMSs

    print(command)
    run(command)
    print('Now correcting for beam...')
    
    f='image_full_vlow_nocut_m.app.restored.fits'
    dapp=fits.open(f)
    dint=fits.open(f.replace('app','int'))                     
    beam=dint[0].data/dapp[0].data                    
    wsuc=fits.open('WSCLEAN_low-MFS-image.fits')
    wsuc[0].data*=beam[0,0,12:-13,12:-13]                 
    wsuc.writeto('WSCLEAN_low-MFS-image-int.fits')  

def compress_fits(filename,q):
    if not os.path.exists(filename+'.fz'):
        command='fpack -q %i %s' % (q,filename)
        run(command)

def do_highres_pol(outname,field,options=None):
    # Makes a single 6" QU cube for entire field using all epochs
    uvrange=[o['image_uvmin'],o['uvmax']]
    ddf_kw = {}
    cubefiles = ['image_full_polhigh_StokesQ.cube.dirty.fits','image_full_polhigh_StokesQ.cube.dirty.corr.fits','image_full_polhigh_StokesU.cube.dirty.fits','image_full_polhigh_StokesU.cube.dirty.corr.fits']
    cthreads=[]
    flist=[]
    ddf_kw = {}
    mslistname='big-mslist.txt'
    do_polcubes('DATA','[DDS3_full_smoothed,DDS3_full_slow]',uvrange,outname,mslistname,ddf_kw,beamsize=o['final_psf_arcsec'],imsize=o['imsize'],cellsize=o['cellsize'],robust=o['final_robust'],options=o,catcher=None)
    #,[0.1,1000.0],'image_full_polhigh',mslistname,ddf_kw,beamsize=6.0,imsize=o['imsize'],cellsize=1.5,robust=-0.5,options=o,catcher=None)

    # Redo the headers
    for cubefile in cubefiles:
        stokesparam = cubefile.split('_')[3].split('.')[0].replace('Stokes','')
        redo_cube_headers(cubefile,cubefile,stokesparam)
    
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
    if o['compress_polcubes']:
        # cthreads and flist exist
        for thread in cthreads:
            if thread.is_alive():
                warn('Waiting for a compression thread to finish')
                thread.join()
        #if o['delete_compressed']:
        #    for f in flist:
        #        warn('Deleting compressed file %s' % f)
        #        #os.remove(f)
    #update_status(None,'Complete')

def do_epoch_pol(outname,mslistname,options=None):

    low_uvrange=[o['image_uvmin'],2.5*206.0/o['low_psf_arcsec']]
    if o['low_imsize'] is not None:
        low_imsize=o['low_imsize'] # allow over-ride
    else:
        low_imsize=o['imsize']*o['cellsize']/o['low_cell']
    cubefiles=['%s_QU.cube.dirty.fits'%outname,'%s_QU.cube.dirty.corr.fits'%outname]
    cthreads=[]
    flist=[]
    ddf_kw = {}
    do_polcubes('DATA','[DDS3_full_smoothed,DDS3_full_slow]',low_uvrange,outname,mslistname,ddf_kw,beamsize=o['low_psf_arcsec'],imsize=low_imsize,cellsize=o['low_cell'],robust=o['low_robust'],options=o,catcher=None)

    # Redo the headers
    #for cubefile in cubefiles:
    #    stokesparam = cubefile.split('_')[3].split('.')[0].replace('Stokes','')
    #    redo_cube_headers(cubefile,cubefile,stokesparam)
    
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
    if o['compress_polcubes']:
        # cthreads and flist exist
        for thread in cthreads:
            if thread.is_alive():
                warn('Waiting for a compression thread to finish')
                thread.join()
        #if o['delete_compressed']:
        #    for f in flist:
        #        warn('Deleting compressed file %s' % f)
        #        #os.remove(f)
    #update_status(None,'Complete')

    
def do_run_high_v(outname,mslist,options=None):
    uvrange=[o['image_uvmin'],o['uvmax']]
    ddf_kw = {}
    ddf_image(outname,mslist,
                  cleanmode='SSD',ddsols='[DDS3_full_smoothed,DDS3_full_slow]',
                  applysols=o['apply_sols'][6],stokes='IV',
                  AllowNegativeInitHMP=True,
                  majorcycles=0,robust=o['final_robust'],
                  colname='DATA',use_dicomodel=False,
                  uvrange=uvrange,cellsize=o['cellsize'],
                  peakfactor=0.001,
                  smooth=True,automask=True,automask_threshold=5,normalization=o['normalize'][2],
                  catcher=None,options=o,**ddf_kw)

def make_mask(imagename,thresh):

    fname=imagename+'.mask.fits'
    runcommand = "MakeMask.py --RestoredIm=%s --Th=%s --Box=50,2"%(imagename,thresh)
    run(runcommand)
    return fname

def update_status(name,operation,status,time=None,workdir=None,av=None,survey=None):
    # modified from surveys_db.update_status
    # utility function to just update the status of a field
    # name can be None (work it out from cwd), or string (field name)

    with SurveysDB(survey=survey) as sdb:
        idd=sdb.get_ffr(name,operation)
        if idd is None:
            raise RuntimeError('Unable to find database entry for field "%s".' % id)
        idd['status']=status
        tag_field(sdb,idd,workdir=workdir)
        if time is not None and idd[time] is None:
            idd[time]=datetime.datetime.now()
        sdb.set_ffr(idd)

def get_next(reverse=False):
    with SurveysDB(readonly=True) as sdb:
        sdb.cur.execute('select distinct ffr.id,ffr.priority from full_field_reprocessing ffr left join fields f on ffr.id=f.id where ffr.status="Not started" and ffr.clustername is NULL order by priority desc,ra'+' desc' if reverse)
        results=sdb.cur.fetchall()
    if len(results)==0:
        return None
    else:
        return results[0]['id']
        
def get_galactic_coords(msfile):
    """Return galactic longitude and latitude (degrees) from an MS file."""
    import pyrap.tables as pt
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    with pt.table(msfile + '/FIELD', readonly=True) as field_table:
        phase_dir = field_table.getcol('PHASE_DIR')[0, 0]  # shape: (n_fields, 1, 2)
        ra_rad, dec_rad = phase_dir
        ra_deg = ra_rad * 180.0 / np.pi
        dec_deg = dec_rad * 180.0 / np.pi

    skycoord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    gal = skycoord.galactic
    return gal.l.deg, gal.b.deg

def upload_results(operation, field, resultfilestar=None, nodb=False):
    # Lookup table of where the operations go. Next time, Tim, let's
    # call the directories the same thing as the operations (-:
    # If resultfilestar is None then try to find the files to upload
    # This allows recovery from a failed upload
    op_locations={'FullSub':'Subtracted_data',
                  'TransientImage':'Subtracted_snapshot_images',
                  'DynSpecMS':'DynSpecMS',
                  'StokesV':'StokesV_imaging',
                  'HighPol':'Pol_highres',
                  'EpochPol':'Pol_Epoch',
                  'VLow_image':'VLow_imaging',
                  'VLow_sub_image':'VLow_sub_imaging'
                  }
    glob_patterns={'DynSpecMS':'DynSpecs_*.tar',
                   'FullSub':'*_object*archive*.tar',
                   'TransientImage':'*_snapshot_images.tar',
                   'StokesV':'V_high_maps.tar',
                   'HighPol':'stokes_highres.tar',
                   'EpochPol':'stokes_epochpol.tar',
                   'VLow_image':'*_low_images.tar',
                   'VLow_sub_image':'*_low_sub_images.tar'}
    
    if os.path.isdir(field): os.chdir(field)
    if resultfilestar is None:
        resultfilestar=glob.glob(glob_patterns[operation])
    if operation not in op_locations:
        raise RuntimeError('Upload location undefined')
    target=op_locations[operation]+'/'
    if not nodb: update_status(field,operation,'Uploading')
    result=do_rclone_reproc_tape_upload(field,os.getcwd(),resultfilestar,target)
    if result['code']==0:
        if not nodb: update_status(field,operation,'Verified',time='end_date')
    else:
        if not nodb: update_status(field,operation,'Upload failed')
        raise RuntimeError(operation+' upload failed')


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Reprocessing LoTSS fields')
    parser.add_argument('--StokesV', help='Include Stokes V reprocessing', action='store_true')
    parser.add_argument('--FullSub', help='Include full field subtraction', action='store_true')
    parser.add_argument('--HighPol', help='Include full field 6asec QU cube', action='store_true')
    parser.add_argument('--EpochPol', help='Create individual epoch QU cube at 20" and 45"', action='store_true')    
    parser.add_argument('--Dynspec', help='Process with DynSpecMS', action='store_true')
    parser.add_argument('--Field',help='LoTSS fieldname',type=str,default="")
    parser.add_argument('--Cache',help='DDF cache dir',type=str,default=".")
    parser.add_argument('--NoDBSync',help='Do not update the reprocessing database',action='store_true')
    parser.add_argument('--NoUpload',help='Do not upload the processed files',action='store_true')
    parser.add_argument('--NCPU',help='Number of CPU (0 for all)',type=int,default=0)
    parser.add_argument('--Force', help='Process anyway disregarding status in database',action='store_true')
    parser.add_argument('--TransientImage',help='Only possible if doing FullSub and then images output data at 1 image per time slot',action='store_true')
    parser.add_argument('--IgnorePrerequisites',help='Ignore sequence prerequisites, e.g. if the data already exist',action='store_true')
    parser.add_argument('--VLow_image',help='Image the data at very low resolution with DDFacet',action='store_true')
    parser.add_argument('--VLow_sub_image',help='Image the source subtracted data at very low resolution with WSClean',action='store_true')
    parser.add_argument('--Download_FullSub',help='Download the full subtracted data from tape',action='store_true')
    parser.add_argument('--Skip_DataPrepare',help='Do not download the uv-data (if e.g. working on downloaded subtracted data)',action='store_true')
    parser.add_argument('--DoDatabaseOps',help='Do the operations specified in the database (in this case operations to do should have status other than Verified)',action='store_true')
    parser.add_argument('--Upload',help='Retry failed uploads for this field, then stop',action='store_true')
    parser.add_argument('--GetNext',help='Get the next field to do in the database, and do that (combine with DoDatabaseOps)',action='store_true')
    parser.add_argument('--Reverse',help='Do GetNext in reverse RA order',action='store_true')

    args = vars(parser.parse_args())

    args['DynSpecMS']=args['Dynspec'] ## because option doesn't match database value
    if args['NCPU']==0: args['NCPU']=getcpus()
    print('Input arguments: ',args)
    if args['GetNext']:
        field=get_next(reverse=args.Reverse)
        print('Selected field is',field)
        if not field:
            raise RuntimeError('Cannot find a field to do!')
    else:
        field = args['Field']

    if not args['NoDBSync']:
        with SurveysDB(readonly=True) as sdb:
            sdb.cur.execute('select * from full_field_reprocessing where id="%s"'%field)
            # sdb.cur.execute('select * from fields where status=Verified and id="%s"'%field)
            results=sdb.cur.fetchall()
            print('Requested field database:',results)
        if len(results)==0:
            raise RuntimeError('Requested field is not in database')
        
        with SurveysDB(readonly=True) as sdb:
            fieldinfo=sdb.get_field(field)
            print('Field info',fieldinfo)
            if fieldinfo is None:
                raise RuntimeError('Field',field,'does not exist in the database')

        if args['Upload']:
            # Special mode to retry a failed upload
            field=args['Field']
            for r in results:
                if r['status']=='Upload failed':
                    print('Retrying upload for operation',r['operation'])
                    upload_results(r['operation'],field)
            sys.exit(0)
    
        if args['DoDatabaseOps']:
            for r in results:
                if r['status']=='Verified': continue
                if r['clustername']!=get_cluster() and r['clustername'] is not None:
                    raise RuntimeError('Trying to run non-local fields??')
                print('Enabling operation',r['operation'])
                args[r['operation']]=True
                
    startdir = os.getcwd()

    for option in ['StokesV','FullSub','HighPol','DynSpecMS','EpochPol','TransientImage','VLow_image','VLow_sub_image']:
        if args[option] and not args['NoDBSync']:
            with SurveysDB(readonly=False) as sdb:
                tmp = sdb.get_ffr(field,option)
                if tmp['status'] not in ['Not started','Staging','Staged','Downloaded','Unpacked','Queued'] or (tmp['clustername'] is not None and tmp['clustername']!=get_cluster()):
                    print('Status of',option,tmp['status'])
                    if not args["NoDBSync"] and not args['Force']:
                        raise RuntimeError('Field already processing')

    if not os.path.exists(startdir+'/'+field):
        os.system('mkdir %s'%field)
        os.chdir(field)
        print('Downloading field',field)
        for option in ['StokesV','FullSub','HighPol','DynSpecMS','EpochPol','TransientImage','VLow_image','VLow_sub_image']:
            if args[option] and not args["NoDBSync"]:
                update_status(field,option,'Downloading')
        if not args['Skip_DataPrepare']:
            prepare_field(field,startdir +'/'+field)
        else:
            prepare_field(field,startdir + '/' + field,Mode='Misc',operations=['download','untar'])

    else:
        os.chdir(field)

    # should be downloaded and ready to go, set status to Queued now,
    # then each operation will be set to started as we go
    for option in ['StokesV','FullSub','HighPol','DynSpecMS','EpochPol','TransientImage','VLow_image','VLow_sub_image']:
        if args[option] and not args["NoDBSync"]:
            update_status(field,option,'Queued')


    print('Reading summary file (summary.txt)')
    if not os.path.exists('summary.txt'):
        print('Cannot read the summary.txt file - failing')
        sys.exit(0)
    o = convert_summary_cfg(option_list)
    o['NCPU_DDF'] = args['NCPU']
    o['NCPU_killms'] = args['NCPU']
    o['colname'] = 'DATA'
    print(o)
        
    if not args['Skip_DataPrepare']:
        print('Checking big-mslist')
        m=MSList('big-mslist.txt')
        print('Looking for different epochs')
        uobsid = set(m.obsids)
        epoch_mslists=[]
        for obsid in uobsid:
            umslist='mslist-%s.txt' % obsid
            epoch_mslists.append(umslist)
            print('Writing ms list for obsids',umslist)
            with open(umslist,'w') as file:
                for ms,ob in zip(m.mss,m.obsids):
                    if ob==obsid:
                        file.write(ms+'\n')
    
    #for option in ['StokesV','FullSub','HighPol','DynSpecMS','EpochPol','TransientImage','VLow_image','VLow_sub_image']:
    #    if args[option] and not args["NoDBSync"]:
    #        print('Changing',option,'status to Started','for',field)
    #        update_status(field,option,'Started',time='start_date')

    if args['Download_FullSub']:
        do_rclone_reproc_tape_download(field,'./','Subtracted_data',verbose=False)

    if args['FullSub']:
        if args['Download_FullSub']:
            print('Full subtracted data already downloaded, skipping source subtraction')
        else:
            if not args["NoDBSync"]: update_status(field,'FullSub','Started',time='start_date')
            do_run_subtract(field)
            resultfiles = glob.glob('*object*sub*archive?')
            if not args["NoDBSync"]:
                resultfilestar = []
                for resultfile in resultfiles:

                     p = subprocess.run(
                       [
                         "tar",
                         "--exclude=table.lock",
                         "--exclude=*.lock",
                         "-cvf", '%s.tar'%resultfile,
                         resultfile,
                      ],
                      text=True,
                      capture_output=True,
                      )

                     # GNU tar: 0=ok, 1=warnings (e.g. "file changed as we read it"), 2=fatal
                     if p.returncode == 2:
                         raise RuntimeError(
                         f"Tar of {resultfile} failed (fatal).\n"
                         f"STDERR:\n{p.stderr}\n"
                         f"STDOUT:\n{p.stdout}"
                         )
                     elif p.returncode == 1:
                         # Warning only    ^`^t keep going, but print so it lands in Slurm logs
                         print(f"Tar warning for {resultfile}:\n{p.stderr}")
                     resultfilestar.append('%s.tar'%resultfile)

                if not args["NoUpload"]:
                     update_status(field,'FullSub','Uploading')
                     result = do_rclone_reproc_tape_upload(field,os.getcwd(),resultfilestar,'Subtracted_data/')
                     if result['code']==0:
                        update_status(field,'FullSub','Verified',time='end_date')
                     else:
                        update_status(field,'FullSub','Upload failed')
   
    if args['TransientImage']:
        if not args["NoDBSync"]: update_status(field,'TransientImage','Started',time='start_date')
        resultfiles = glob.glob('*object*sub*archive?')
        if len(resultfiles)==0:
            print('Source subtraction (--FullSub) is needed to do TransientImage')
            sys.exit(0)
        resultdir=field+'_snapshot_images'
        try:
            os.mkdir(resultdir)
        except FileExistsError:
            pass # should only happen on a restart, and then files moved in there should be kept
        for resultfile in resultfiles:
            imagefile = resultfile.split('_')[0] + '_epoch_' + resultfile.split('archive')[-1]
            print(resultfile)
            print(imagefile)
            gal_l, gal_b = get_galactic_coords(resultfiles[0])
            if abs(gal_b) < 10.0:
                transient_image(resultfile,imagefile,galactic=True,options=o)
            else:
                transient_image(resultfile,imagefile,galactic=False,options=o)
            g=glob.glob(imagefile+'*dirty.fits.fz')
            assert(g)
            for f in g:
                shutil.move(f,resultdir)
        if not args['NoUpload']:
            d=os.system('tar -cvf %s_snapshot_images.tar %s'%(field,resultdir))
            if d!=0:
                raise RuntimeError('Tar of %s_snapshot_images failed'%field)	
            upload_results('TransientImage',field,['%s_snapshot_images.tar'%field],nodb=args['NoDBSync'])

    if args['Dynspec'] or args['DynSpecMS']:
        if not args["NoDBSync"]: update_status(field,'DynSpecMS','Started',time='start_date')
        do_run_dynspec(field)
        ingest_dynspec()
        OutDir="DynSpecs_%s"%field
        os.system("mkdir -p %s"%OutDir)
        resultfiles = glob.glob('*DynSpecs_*.tgz')
        if not resultfiles:
            update_status(field,'DynSpecMS','Failed')
            raise RuntimeError('Failed to find any results file from DynSpecMS!')
        for resultfile in resultfiles:
            os.system('cp %s %s'%(resultfile,OutDir))
        if not args['NoUpload']:
            os.system('tar -cvf %s.tar %s'%(OutDir,OutDir))
            upload_results('DynSpecMS',field,['%s.tar' % OutDir],nodb=args['NoDBSync'])
            
    if args['StokesV']:
        if not args["NoDBSync"]: update_status(field,'StokesV','Started',time='start_date')
        for obsid in uobsid:
            print('Stokes V image for %s'%obsid)
            do_run_high_v('image_full_high_stokesV_%s'%obsid,'mslist-%s.txt'%obsid,options=o)

            resultfiles = glob.glob('image_full_high_stokesV_%s*dirty*.fits'%obsid)
            os.system('mkdir V_high_maps')
            for resultfile in resultfiles:
                os.system('cp %s V_high_maps'%(resultfile))
        if not args['NoUpload']:
            os.system('tar -cvf V_high_maps.tar V_high_maps')
            resultfilestar = ['V_high_maps.tar']            
            upload_results('StokesV',field,resultfilestar,nodb=args['NoDBSync'])

    if args['HighPol']:
        if not args["NoDBSync"]: update_status(field,'HighPol','Started',time='start_date')
        from do_polcubes import do_polcubes
        os.system('mkdir logs')
        do_highres_pol('image_full_polhigh',field,options=o)
        resultfiles = glob.glob('*fz')
        print('Compressed pol cubes',resultfiles)
        os.system('mkdir stokes_highres')
        for resultfile in resultfiles:
            os.system('mv %s stokes_highres/'%(resultfile))
        if not args['NoUpload']:
            os.system('tar -cvf stokes_highres.tar stokes_highres')
            upload_results('HighPol',field,['stokes_highres.tar'],nodb=args['NoDBSync'])

    if args['EpochPol']:
        if not args["NoDBSync"]: update_status(field,'EpochPol','Started',time='start_date')

        from do_polcubes import do_polcubes
        os.system('mkdir logs')
        for obsid in uobsid:
            do_epoch_pol('image_full_low_QU_%s'%obsid,'mslist-%s.txt'%obsid,options=o)
            resultfiles = glob.glob('*%s*fz'%obsid)
            print('Compressed pol cubes',resultfiles)
            os.system('mkdir stokes_epochpol')
            for resultfile in resultfiles:
                os.system('mv %s stokes_epochpol/'%(resultfile))
        if not args['NoUpload']:
            os.system('tar -cvf stokes_epochpol.tar stokes_epochpol')
            upload_results('EpochPol',field,['stokes_epochpol.tar'],nodb=args['NoDBSync'])

    if args['VLow_image']:
        if not args["NoDBSync"]: update_status(field,'VLow_image','Started',time='start_date')
        image_vlow(args['NCPU'],cache=args['Cache'])
        OutDir = '%s_low_images'%field
        os.system('mkdir %s'%OutDir)
        resultfiles = glob.glob('image_full_vlow_nocut*')
        for resultfile in resultfiles:
            os.system('cp %s %s'%(resultfile,OutDir))
        if not args['NoUpload']:
            os.system('tar -cvf %s.tar %s'%(OutDir,OutDir))
            resultfilestar = ['%s.tar'%OutDir]
            upload_results('VLow_image',field,resultfilestar,nodb=args['NoDBSync'])

    if args['VLow_sub_image']:
        if not args["NoDBSync"]: update_status(field,'VLow_sub_image','Started',time='start_date')
        image_vlow_sub() # does a new subtraction so must be run last
        OutDir = '%s_low_sub_images'%field
        os.system('mkdir %s'%OutDir)
        resultfiles = glob.glob('WSCLEAN_low*')
        for resultfile in resultfiles:
            os.system('cp %s %s'%(resultfile,OutDir))
        if not args['NoUpload']:
            os.system('tar -cvf %s.tar %s'%(OutDir,OutDir))
            resultfilestar = ['%s.tar'%OutDir]
            upload_results('VLow_sub_image',field,resultfilestar,nodb=args['NoDBSync'])

