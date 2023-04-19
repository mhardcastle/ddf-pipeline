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
from reprocessing_utils import prepare_field,do_rclone_disk_upload,do_rclone_tape_pol_upload
import argparse
import threading
from auxcodes import run,warn,report
import numpy as np
import pipeline
import datetime

from rclone import RClone
from sdr_wrapper import SDR
from time import sleep

def stage_field(cname,f,verbose=False,Mode='Imaging+Misc'):
    # stage a dataset from SDR or rclone repos 
    # this should really live in reprocessing_utils and share code with do_sdr_and_rclone_download -- low-level stage functionality should be moved to sdr_wrapper
    # as with that function cname is the field name, f is the processing directory
    if not os.path.isdir(f):
        os.makedirs(f)
    s=SDR(target=f)
    try:
        files=s.get_status(cname)
    except RuntimeError:
        files=None
    if files:
        if verbose: print('Initiating SDR stage for field',cname)
        if Mode=="Imaging":
            tarfiles=['images.tar','uv.tar']
        elif Mode=="Misc":
            tarfiles=['misc.tar']
        elif Mode=="Imaging+Misc":
            tarfiles=['images.tar','uv.tar','misc.tar',"stokes_small.tar"]

        # code adapted from download_and_stage
        for f in tarfiles:
            if f not in files:
                raise RuntimeError('File not found!')
            else:
                if files[f]=='OFL':
                    s.stage(cname,f)

        if verbose:
            print('Waiting for files to be online:')
        
        while True:
            files=s.get_status(cname)
            count=0
            for f in tarfiles:
                if files[f]=='DUL':
                    count+=1
            if verbose:
                print('%i/%i... ' % (count,len(tarfiles)),end='')
                sys.stdout.flush()
            if count==len(tarfiles):
                if verbose: print()
                break
            else:
                sleep(30)

    else:
        # staging for rclone goes here.
        pass
    

def update_status(name,operation,status,time=None,workdir=None,survey=None):
    # modified from surveys_db update_status
    # utility function to just update the status of a ffr entry
    # name can be None (work it out from cwd), or string (field name)

    with SurveysDB(survey=survey) as sdb:
        idd=sdb.get_ffr(id,operation)
        if idd is None:
            raise RuntimeError('Unable to find database entry for field "%s".' % id)
        idd['status']=status
        tag_field(sdb,idd,workdir=workdir)
        if time is not None and idd[time] is None:
            idd[time]=datetime.datetime.now()
        sdb.set_ffr(idd)

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
        executionstr = 'ms2dynspec.py --ms big-mslist.txt --data DATA --model PREDICT_SUB --sols [DDS3_full_smoothed,DDS3_full_slow] --rad 2. --SolsDir SOLSDIR --BeamModel LOFAR --BeamNBand 1 --DicoFacet image_full_ampphase_di_m.NS_SUB.DicoFacet --noff 100 --nMinOffPerFacet 5 --CutGainsMinMax 0.1,1.5 --SplitNonContiguous 1 --imageI image_full_ampphase_di_m.NS.int.restored.fits --imageV image_full_high_stokesV.dirty.corr.fits --SavePDF 1 --FitsCatalog ${DDF_PIPELINE_CATALOGS}/dyn_spec_catalogue_addedexo_addvlotss.fits'
        print(executionstr)
        result=os.system(executionstr)
        if result!=0:
            raise RuntimeError('ms2dynspec.py failed with error code %i' % result)
    else:
        print("All DynSpecMS output exists, skipping... ")    



def compress_fits(filename,q):
    command='fpack -q %i %s' % (q,filename)
    run(command)

def do_highres_pol(field):

    cubefiles = ['image_full_polhigh_StokesQ.cube.dirty.fits','image_full_polhigh_StokesQ.cube.dirty.corr.fits','image_full_polhigh_StokesU.cube.dirty.fits','image_full_polhigh_StokesU.cube.dirty.corr.fits']
    cthreads=[]
    flist=[]
    ddf_kw = {}
    do_polcubes('DATA','[DDS3_full_smoothed,DDS3_full_slow]',[0.1,1000.0],'image_full_polhigh',ddf_kw,beamsize=6.0,imsize=o['imsize'],cellsize=1.5,robust=-0.5,options=o,catcher=None)

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

    

    
def do_run_high_v(field):

    executionstr = 'DDF.py --Parallel-NCPU=12 --Output-Name=image_full_high_stokesV --Data-MS=big-mslist.txt --Deconv-PeakFactor 0.001000 --Data-ColName DATA --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=0 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-PhasedArrayMode=A --Weight-Robust -0.50000 --Image-NPix=20000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 1.500000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --GAClean-AllowNegativeInitHMP True --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 6.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-PolMode=IV --Output-Mode=Dirty --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=5.00 --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=[DDS3_full_smoothed,DDS3_full_slow] --Selection-UVRangeKm=[0.100000,1000.0000] --GAClean-MinSizeInit=10 --Beam-Smooth=1'

    print(executionstr)
    result=os.system(executionstr)
    if result!=0:
        raise RuntimeError('sub-sources-outside-region.py failed with error code %i' % result)

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

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Reprocessing LoTSS fields')
    parser.add_argument('--StokesV', help='Include Stokes V reprocessing', action='store_true')
    parser.add_argument('--FullSub', help='Include full field subtraction', action='store_true')
    parser.add_argument('--HighPol', help='Include full field 6asec QU cube', action='store_true')
    parser.add_argument('--Dynspec', help='Process with DynSpecMS', action='store_true')
    parser.add_argument('--Field',help='LoTSS fieldname',type=str,default="")
    parser.add_argument('--Parset',help='DDF pipeline parset',type=str)
    parser.add_argument('--NoDBSync',help='DDF pipeline parset',type=int,default=0)
    args = vars(parser.parse_args())
    args['DynSpecMS']=args['Dynspec'] ## because option doesn't match database value
    
    field = args['Field']
    if args['Parset']:
        o = options(args['Parset'],option_list)
        print(o)
    print('Input arguments: ',args)


    with SurveysDB(readonly=True) as sdb:
        sdb.cur.execute('select * from full_field_reprocessing where id="%s"'%field)
        # sdb.cur.execute('select * from fields where status=Verified and id="%s"'%field)
        results=sdb.cur.fetchall()
        print('Requested field database:',results)
    if len(results)==0:
        raise RuntimeError('Requested field is not in database')

    startdir = os.getcwd()

    for option in ['StokesV','FullSub','HighPol','DynSpecMS']:
        if args[option]:
            with SurveysDB(readonly=False) as sdb:
                tmp = sdb.get_ffr(field,option)
                if tmp['status'] not in ['Not started','Staged','Downloaded','Unpacked','Queued'] or (tmp['clustername'] is not None and tmp['clustername']!=get_cluster()):
                    print('Status of',option,tmp['status'])
                    if not args["NoDBSync"]:
                        raise RuntimeError('Field already processing')

    if not os.path.exists(startdir+'/'+field):
        os.system('mkdir %s'%field)
        os.chdir(field)
        print('Downloading field',field)
        for option in ['StokesV','FullSub','HighPol','DynSpecMS']:
            if args[option] and not args["NoDBSync"]:
                update_status(field,option,'Downloading')
        prepare_field(field,startdir +'/'+field)
    else:
        os.chdir(field)

    for option in ['StokesV','FullSub','HighPol','DynSpecMS']:
        if args[option] and not args["NoDBSync"]:
            print('Changing',option,'status to Started')
            update_status(field,option,'Started',time='start_date')

    if args['FullSub']:
        do_run_subtract(field)

        # resultfiles = glob.glob('*sub*archive*')
        # resultfilestar = []
        # for resultfile in resultfiles:
        #     d=os.system('tar -cvf %s.tar %s'%(resultfile,resultfile))
        #     if d!=0:
        #         raise RuntimeError('Tar of %s failed'%resultfile)	
        #     resultfilestar.append('%s.tar'%resultfile)

        # do_rclone_disk_upload(field,os.getcwd(),resultfilestar,'subtract_pipeline/')

        # with SurveysDB(readonly=False) as sdb:
        #     tmp = sdb.get_ffr(field,'FullSub')
        #     tmp['status'] == 'Verified'
        #     sdb.set_ffr(tmp)

    if args['Dynspec']:
        
        do_run_dynspec(field)
        
        pipeline.ingest_dynspec()

        OutDir="DynSpecs_%s"%field
        os.system("mkdir -p %s"%OutDir)
        resultfiles = glob.glob('DynSpecs_*.tgz')
        for resultfile in resultfiles:
            os.system('cp %s %s'%(resultfile,OutDir))
        os.system('tar -cvf %s.tar %s'%(OutDir,OutDir))
        resultfilestar = ['%s.tar'%OutDir]

        if not args["NoDBSync"]:
            update_status(field,'DynSpecMS','Uploading')
            result=do_rclone_disk_upload(field,os.getcwd(),resultfilestar,'DynSpecMS_reprocessing')
            if result['code']==0:
                update_status(field,'DynSpecMS','Verified',time='end_date')
            else:
                update_status(field,'DynSpecMS','Upload failed')

            
    if args['StokesV']:
        do_run_high_v(field)

        resultfiles = glob.glob('image_full_high_stokesV*dirty*.fits')
        os.system('mkdir V_high_maps')
        for resultfile in resultfiles:
            os.system('cp %s V_high_maps'%(resultfile))
        os.system('tar -cvf V_high_maps.tar V_high_maps')
        resultfilestar = ['V_high_maps.tar']

            
        if not args["NoDBSync"]:
            update_status(field,'StokesV','Uploading')
            do_rclone_disk_upload(field,os.getcwd(),resultfilestar,'Stokes_V_imaging')
            update_status(field,'StokesV','Verified',time='end_date')

    if args['HighPol']:
        from do_polcubes import do_polcubes
        os.system('mkdir logs')
        do_highres_pol(field)
        resultfiles = glob.glob('*fz')
        print('Compressed pol cubes',resultfiles)
        os.system('mkdir stokes_highres')
        for resultfile in resultfiles:
            os.system('mv %s stokes_highres/'%(resultfile))
        os.system('tar -cvf stokes_highres.tar stokes_highres')

        if not args["NoDBSync"]:
            resultfilestar = ['stokes_highres.tar']
            print('Starting upload of',resultfilestar)
            update_status(field,'HighPol','Uploading')
            do_rclone_tape_pol_upload(field,os.getcwd(),resultfilestar,'')
            update_status(field,'HighPol','Verified',time='end_date')
