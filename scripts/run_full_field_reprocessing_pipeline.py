#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from surveys_db import update_reprocessing_extract, get_next_extraction, SurveysDB
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
from reprocessing_utils import *
import argparse
import threading
from auxcodes import run,warn,report
import numpy as np

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

    executionstr = 'DDF.py --Parallel-NCPU=12 --Output-Name=image_full_high_stokesV --Data-MS=big-mslist.txt --Deconv-PeakFactor 0.001000 --Data-ColName DATA --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=0 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust -0.50000 --Image-NPix=20000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 1.500000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --GAClean-AllowNegativeInitHMP True --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 6.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-PolMode=IV --Output-Mode=Dirty --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=5.00 --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=[DDS3_full_smoothed,DDS3_full_slow] --Selection-UVRangeKm=[0.100000,1000.0000] --GAClean-MinSizeInit=10 --Beam-Smooth=1'

    print(executionstr)
    result=os.system(executionstr)
    if result!=0:
        raise RuntimeError('sub-sources-outside-region.py failed with error code %i' % result)



if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Reprocessing LoTSS fields')
    parser.add_argument('--StokesV', help='Include Stokes V reprocessing', action='store_true')
    parser.add_argument('--FullSub', help='Include full field subtraction', action='store_true')
    parser.add_argument('--HighPol', help='Include full field 6asec QU cube', action='store_true')
    parser.add_argument('--Field',help='LoTSS fieldname',type=str)
    parser.add_argument('--Parset',help='DDF pipeline parset',type=str)
    args = vars(parser.parse_args())

    field = args['Field']
    o = options(args['Parset'],option_list)
    print(o)
    print('Input arguments: ',args)


    with SurveysDB(readonly=True) as sdb:
        sdb.cur.execute('select * from full_field_reprocessing where id="%s"'%field)
        results=sdb.cur.fetchall()
        print('Requested field database:',results)
    if len(results)==0:
        raise RuntimeError('Requested field is not in database')

    startdir = os.getcwd()

    if not os.path.exists(startdir+'/'+field):
        os.system('mkdir %s'%field)
        os.chdir(field)
        print('Downloading field',field)
        prepare_field(field,startdir +'/'+field)
    else:
        os.chdir(field)

    if args['StokesV']:
        with SurveysDB(readonly=False) as sdb:
            tmp = sdb.get_ffr(field,'StokesV')
            if tmp['status'] != 'Not started':
                print('Status of StokesV',tmp['status'])
                raise RunTimeError('Field already processing')
            tmp['status'] = 'Started'
            print('Changing StokesV status to Started')
            sdb.set_ffr(tmp)
        
    if args['FullSub']:
        with SurveysDB(readonly=False) as sdb:
            tmp = sdb.get_ffr(field,'FullSub')
            if tmp['status'] != 'Not started':
                print('Status of FullSub',tmp['status'])
                raise RunTimeError('Field already processing')
            tmp['status'] = 'Started'
            print('Changing FullSub status to Started')
            sdb.set_ffr(tmp)

    if args['HighPol']:
        with SurveysDB(readonly=False) as sdb:
            tmp = sdb.get_ffr(field,'HighPol')
            if tmp['status'] != 'Not started':
                print('Status of HighPol',tmp['status'])
                raise RunTimeError('Field already processing')
            tmp['status'] = 'Started'
            print('Changing HighPol status to Started')
            sdb.set_ffr(tmp)


    if args['FullSub']:
        do_run_subtract(field)

        resultfiles = glob.glob('*sub*archive*')
        resultfilestar = []
        for resultfile in resultfiles:
            d=os.system('tar -cvf %s.tar %s'%(resultfile,resultfile))
            if d!=0:
                raise RuntimeError('Tar of %s failed'%resultfile)	
            resultfilestar.append('%s.tar'%resultfile)


        do_rclone_disk_upload(field,os.getcwd(),resultfilestar,'subtract_pipeline/')

        with SurveysDB(readonly=False) as sdb:
            tmp = sdb.get_ffr(field,'FullSub')
            tmp['status'] == 'Verified'
            sdb.set_ffr(tmp)

    if args['StokesV']:
        do_run_high_v(field)

        resultfiles = glob.glob('image_full_high_stokesV*dirty*.fits')
        os.system('mkdir V_high_maps')
        for resultfile in resultfiles:
            os.system('cp %s V_high_maps'%(resultfile))
        os.system('tar -cvf V_high_maps.tar V_high_maps')
        resultfilestar = ['V_high_maps.tar']

        do_rclone_disk_upload(field,os.getcwd(),resultfilestar,'Stokes_V_imaging')

        with SurveysDB(readonly=False) as sdb:
            tmp = sdb.get_ffr(field,'StokesV')
            tmp['status'] == 'Verified'
            sdb.set_ffr(tmp)

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

        resultfilestar = ['stokes_highres.tar']
        print('Starting upload of',resultfilestar)
        do_rclone_tape_pol_upload(field,os.getcwd(),resultfilestar,'')

        with SurveysDB(readonly=False) as sdb:
            tmp = sdb.get_ffr(field,'HighPol')
            tmp['status'] == 'Verified'
            sdb.set_ffr(tmp)
        # SEe /net/lofar8/data2/shimwell/testing-highres-pol
