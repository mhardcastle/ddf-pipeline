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
from reprocessing_utils import prepare_field,do_rclone_disk_upload,do_rclone_tape_pol_upload,convert_summary_cfg
import argparse
import threading
from auxcodes import run,warn,report, MSList
import numpy as np
import pipeline
import datetime
from pipeline import ddf_image
from rclone import RClone
from sdr_wrapper import SDR
from time import sleep
import pyrap.tables as pt


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

def transient_image(msfilename,imagename,galactic=False,options=None):
    t = pt.table(msfilename,readonly=True)
    numtimes = len(np.unique(t.getcol('TIME')))
    print(numtimes,msfilename)

    # Create 8s images and 2 min images (assuming time resolution is 8s...). 
    if galactic:
        chanout = 16
    else:
        chanout = 1
    if chanout > 1:
        os.system('wsclean -make-psf -intervals-out %s -padding 1.6 -interval 0 %s -auto-threshold 5 -channels-out %s -deconvolution-channels 3 -fit-spectral-pol 3 -scale 6asec -size 2200 2200 -join-channels -minuv-l 50 -maxuv-l 5000 -mgain 0.8 -multiscale -niter 0 -nmiter 0 -no-update-model-required -parallel-deconvolution 1500 -no-reorder -taper-gaussian 40asec -use-wgridder -weight briggs -0.25 -name %s_8sec %s'%(numtimes,numtimes,chanout,imagename,msfilename))
        os.system('wsclean -make-psf -intervals-out %s -padding 1.6 -interval 0 %s -auto-threshold 5 -channels-out %s -deconvolution-channels 3 -fit-spectral-pol 3 -scale 6asec -size 2200 2200 -join-channels -minuv-l 50 -maxuv-l 5000 -mgain 0.8 -multiscale -niter 0 -nmiter 0 -no-update-model-required -parallel-deconvolution 1500 -no-reorder -taper-gaussian 40asec -use-wgridder -weight briggs -0.25 -name %s_2min %s'%(int(numtimes/15),int(numtimes/15),chanout,imagename,msfilename))
    else:
        os.system('wsclean -make-psf -intervals-out %s -padding 1.6 -interval 0 %s -auto-threshold 5 -channels-out %s -scale 6asec -size 2200 2200 -minuv-l 50 -maxuv-l 5000 -mgain 0.8 -multiscale -niter 0 -nmiter 0 -no-update-model-required -parallel-deconvolution 1500 -no-reorder -taper-gaussian 40asec -use-wgridder -weight briggs -0.25 -name %s_8sec %s'%(numtimes,numtimes,chanout,imagename,msfilename))
        os.system('wsclean -make-psf -intervals-out %s -padding 1.6 -interval 0 %s -auto-threshold 5 -channels-out %s -scale 6asec -size 2200 2200 -minuv-l 50 -maxuv-l 5000 -mgain 0.8 -multiscale -niter 0 -nmiter 0 -no-update-model-required -parallel-deconvolution 1500 -no-reorder -taper-gaussian 40asec -use-wgridder -weight briggs -0.25 -name %s_2min %s'%(int(numtimes/15),int(numtimes/15),chanout,imagename,msfilename))

    outimages = glob.glob('%s*dirty.fits'%imagename)
    for outimage in outimages:
        compress_fits(outimage,o['fpack_q'])


def compress_fits(filename,q):
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
    cubefiles = ['%s_StokesQ.cube.dirty.fits'%outname,'%s_StokesQ.cube.dirty.corr.fits'%outname,'%s_StokesU.cube.dirty.fits'%outname,'%s_StokesU.cube.dirty.corr.fits'%outname]
    cthreads=[]
    flist=[]
    ddf_kw = {}
    do_polcubes('DATA','[DDS3_full_smoothed,DDS3_full_slow]',low_uvrange,outname,mslistname,ddf_kw,beamsize=o['low_psf_arcsec'],imsize=low_imsize,cellsize=o['low_cell'],robust=o['low_robust'],options=o,catcher=None)

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


def update_status(name,operation,status,time=None,workdir=None,av=None,survey=None):
    # modified from surveys_db.update_status
    # utility function to just update the status of a field
    # name can be None (work it out from cwd), or string (field name)

    with SurveysDB(survey=survey) as sdb:
        idd=sdb.get_ffr(name,operation)
        if idd is None:
            raise RuntimeError('Unable to find database entry for field "%s".' % id)
        idd['status']=status
        #tag_field(sdb,idd,workdir=workdir) # Turned this off as columns are not long enough in DB.
        if time is not None and idd[time] is None:
            idd[time]=datetime.datetime.now()
        sdb.set_ffr(idd)

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Reprocessing LoTSS fields')
    parser.add_argument('--StokesV', help='Include Stokes V reprocessing', action='store_true')
    parser.add_argument('--FullSub', help='Include full field subtraction', action='store_true')
    parser.add_argument('--HighPol', help='Include full field 6asec QU cube', action='store_true')
    parser.add_argument('--EpochPol', help='Create individual epoch QU cube at 20" and 45"', action='store_true')    
    parser.add_argument('--Dynspec', help='Process with DynSpecMS', action='store_true')
    parser.add_argument('--Field',help='LoTSS fieldname',type=str,default="")
    parser.add_argument('--NoDBSync',help='Do not Update the reprocessing database',type=int,default=0)
    parser.add_argument('--NCPU',help='Number of CPU',type=int,default=32)
    parser.add_argument('--Force', help='Process anyway disregarding status in database',action='store_true')
    parser.add_argument('--TransientImage',help='Only possible if doing FullSub and then images output data at 1 image per time slot',action='store_true')
    args = vars(parser.parse_args())
    args['DynSpecMS']=args['Dynspec'] ## because option doesn't match database value
    print('Input arguments: ',args)

    field = args['Field']

    with SurveysDB(readonly=True) as sdb:
        sdb.cur.execute('select * from full_field_reprocessing where id="%s"'%field)
        # sdb.cur.execute('select * from fields where status=Verified and id="%s"'%field)
        results=sdb.cur.fetchall()
        print('Requested field database:',results)
    if len(results)==0:
        raise RuntimeError('Requested field is not in database')

    with SurveysDB(readonly=not update_status) as sdb:
        fieldinfo=sdb.get_field(field)
        print('Field info',fieldinfo)
        if fieldinfo is None:
            raise RuntimeError('Field',field,'does not exist in the database')

    startdir = os.getcwd()

    for option in ['StokesV','FullSub','HighPol','DynSpecMS']:
        if args[option]:
            with SurveysDB(readonly=False) as sdb:
                tmp = sdb.get_ffr(field,option)
                if tmp['status'] not in ['Not started','Staged','Downloaded','Unpacked','Queued'] or (tmp['clustername'] is not None and tmp['clustername']!=get_cluster()):
                    print('Status of',option,tmp['status'])
                    if not args["NoDBSync"] and not args['Force']:
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


    print('Reading summary file (summary.txt)')
    if not os.path.exists('summary.txt'):
        print('Cannot read the summary.txt file - failing')
        sys.exit(0)
    o = convert_summary_cfg(option_list)
    o['NCPU_DDF'] = args['NCPU']
    o['NCPU_killms'] = args['NCPU']
    o['colname'] = 'DATA'
    print(o)
        
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
    
    for option in ['StokesV','FullSub','HighPol','DynSpecMS']:
        if args[option] and not args["NoDBSync"]:
            print('Changing',option,'status to Started','for',field)
            update_status(field,option,'Started')#,time='start_date')

    if args['FullSub']:
        #do_run_subtract(field)
        resultfiles = glob.glob('*sub*archive?')
        if not args["NoDBSync"]:
            resultfilestar = []
            for resultfile in resultfiles:
                d=os.system('tar -cvf %s.tar %s'%(resultfile,resultfile))
                if d!=0:
                    raise RuntimeError('Tar of %s failed'%resultfile)	
                resultfilestar.append('%s.tar'%resultfile)
            update_status(field,'FullSub','Uploading')
            result = do_rclone_disk_upload(field,os.getcwd(),resultfilestar,'subtract_pipeline/')
            if result['code']==0:
                update_status(field,'FullSub','Verified')#,time='end_date')
            else:
                update_status(field,'FullSub','Upload failed')
        if args['TransientImage']:
            for resultfile in resultfiles:
                imagefile = resultfile.split('_')[0] + '_epoch_' + resultfile.split('archive')[-1]
                print(resultfile)
                print(imagefile)
                if abs(fieldinfo['gal_b']) < 10.0:
                    transient_image(resultfile,imagefile,galactic=True,options=o)
                else:
                    transient_image(resultfile,imagefile,galactic=False,options=o) 
            os.system('mkdir %s_snapshot_images'%field)
            os.system('mv %s*dirty.fits.fz %s_snapshot_images'%(imagefile,field))
            d=os.system('tar -cvf %s_snapshot_images.tar %s_snapshot_images'%(field,field))
            if d!=0:
                raise RuntimeError('Tar of %s_snapshot_images failed'%field)	
            result = do_rclone_disk_upload(field,os.getcwd(),['%s_snapshot_images.tar'%field],'subtract_snapshot_images/')
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
                update_status(field,'DynSpecMS','Verified')#,time='end_date')
            else:
                update_status(field,'DynSpecMS','Upload failed')

            
    if args['StokesV']:
        for obsid in uobsid:
            print('Stokes V image for %s'%obsid)
            do_run_high_v('image_full_high_stokesV_%s'%obsid,'mslist-%s.txt'%obsid,options=o)

            resultfiles = glob.glob('image_full_high_stokesV_%s*dirty*.fits'%obsis)
            os.system('mkdir V_high_maps')
            for resultfile in resultfiles:
                os.system('cp %s V_high_maps'%(resultfile))
        os.system('tar -cvf V_high_maps.tar V_high_maps')
        resultfilestar = ['V_high_maps.tar']

            
        if not args["NoDBSync"]:
            update_status(field,'StokesV','Uploading')
            do_rclone_disk_upload(field,os.getcwd(),resultfilestar,'Stokes_V_imaging')
            update_status(field,'StokesV','Verified')#,time='end_date')

    if args['HighPol']:
        from do_polcubes import do_polcubes
        os.system('mkdir logs')
        do_highres_pol('image_full_polhigh',field,options=o)
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
            update_status(field,'HighPol','Verified')#,time='end_date')

    if args['EpochPol']:
        from do_polcubes import do_polcubes
        os.system('mkdir logs')
        for obsid in uobsid:
            do_epoch_pol('image_full_low_QU_%s'%obsid,'mslist-%s.txt'%obsid,options=o)
            resultfiles = glob.glob('*fz')
            print('Compressed pol cubes',resultfiles)
            os.system('mkdir stokes_highres')
            for resultfile in resultfiles:
                os.system('mv %s stokes_highres/'%(resultfile))
        os.system('tar -cvf stokes_highres.tar stokes_highres')
