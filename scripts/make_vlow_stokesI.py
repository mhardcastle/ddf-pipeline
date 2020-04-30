#!/usr/bin/env python

from __future__ import print_function
import argparse
import os,sys
from subprocess import call
import pyrap.tables as pt
import numpy as np
import argparse
from make_mslists import make_list
import glob
from astropy.io import ascii
from astropy.io import fits
from scipy.signal import convolve2d
from scipy.ndimage import gaussian_filter
from time import sleep
from auxcodes import run
from getcpus import getcpus
from surveys_db import use_database,SurveysDB,get_id

# Reimage at very low resolution

def die(error,cname=None):
    update_status(cname,'Failed')
    raise RuntimeError(error)

def update_status(name,status):
    if not use_database():
        return
    if name is None:
        # work it out
        id=get_id(workdir=os.getcwd())
    else:
        id=name
        
    with SurveysDB() as sdb:
      idd=sdb.get_field(id)
      if idd is None:
          raise RuntimeError('Unable to find database entry for field "%s".' % id)
      idd['vlow_image']=status
      sdb.set_field(idd)

def make_mask(imagename,thresh):

    fname=imagename+'.mask.fits'
    runcommand = "MakeMask.py --RestoredIm=%s --Th=%s --Box=50,2"%(imagename,thresh)
    run(runcommand)
    return fname            

def get_solutions_timerange(sols):
    t = np.load(sols)['BeamTimes']
    return np.min(t),np.max(t)

def fixsymlinks():
    #dds3smoothed = glob.glob('SOLSDIR/*/*killMS.DDS3_full_smoothed*npz')
    dds3 = glob.glob('SOLSDIR/*/killMS.DDS3_full.sols.npz')
    for i in range(0,len(dds3)):
        symsolname = dds3[i].split('killMS.DDS3_full.sols.npz')[0] + 'killMS.DDS3_full_smoothed.sols.npz' #dds3smoothed[i]
        solname = dds3[i]
        ddsols = 'DDS3_full'
        start_time,t1 = get_solutions_timerange(solname)
        # Rounding different on different computers which is a pain.
        start_time = glob.glob('%s_%s*_smoothed.npz'%(ddsols,int(start_time)))[0].split('_')[2]

        if os.path.islink(symsolname):
            print('Symlink ' + symsolname + ' already exists, recreating')
            os.unlink(symsolname)
            os.symlink(os.path.relpath('../../%s_%s_smoothed.npz'%(ddsols,start_time)),symsolname)
        else:
            print('Symlink ' + symsolname + ' does not yet exist, creating')
            os.symlink(os.path.relpath('../../%s_%s_smoothed.npz'%(ddsols,start_time)),symsolname)
            
    return


def striparchivename():
  mslist = glob.glob('L*_SB*.ms.archive')
  for ms in mslist:
      outname = ms.rstrip('.archive')
      cmd = 'ln -s ' + ms + ' ' + outname
      print (cmd)
      os.system(cmd)

  return


def do_rsync_download(cname,basedir,f):

    workdir=basedir+'/'+cname # directory on remote host to download

    target=os.environ['DDF_PIPELINE_LEIDENUSER']+'@ssh.strw.leidenuniv.nl:'

    while True:
        excludeinclude = ' --include="image_full_ampphase_di_m.NS.mask01.fits" --include="image_full_ampphase_di_m.NS.app.restored.fits" --include="image_full_low_m*fits" --exclude="*QU_*" --exclude="*fits*" --exclude="*.tgz*" --exclude="*QU_*" --exclude="*DDS0*" --exclude="*DDS1*" --exclude="*DDS2*" --exclude="*.corrupted" '
        s= 'rsync -azvh --timeout=20 --progress --perms --chmod=a+rwx'+ excludeinclude + target+workdir + ' ' + f
        #'cd '+workdir+'; rsync -avz --progress --safe-links --inplace --append --partial --timeout=20 '+' '.join(f)+' '+target+'/disks/paradata/shimwell/LoTSS-DR2/archive/'+name
        print('Running command:',s)
        retval=call(s,shell=True)
        if retval==0:
            break
        print('Non-zero return value',retval)
        if retval!=30:
            die('rsync failed unexpectedly',cname)
        sleep(10)


def do_download(cname, basedir='.'):
    # check whether download is needed -- saves rsync
    if os.path.isdir(basedir+'/'+cname):
        os.chdir(basedir+'/'+cname)
        if os.path.isfile('big-mslist.txt'):
            files=[l.rstrip() for l in open('big-mslist.txt').readlines()]
            for f in files:
                if not os.path.exists(f):
                    break
            else:
                return os.getcwd()
    
    update_status(cname,'Downloading')
    os.chdir(basedir)
    do_rsync_download(cname,'/disks/paradata/shimwell/LoTSS-DR2/archive/',os.getcwd())
    os.chdir(cname)
    striparchivename()

    success=make_list(workdir=os.getcwd())
    if not success:
        update_status(cname,'Download failed')
        raise RuntimeError('Failed to make mslist')
    
    #filechecker()
    fixsymlinks()
    update_status(cname,'Downloaded')
    return os.getcwd() # return directory where everything has been done

        
def image_vlow(wd=None):
    if wd is not None:
        os.chdir(wd)
    update_status(None,'Running')
    run('CleanSHM.py')
    run('DDF.py --Output-Name=image_full_vlow_nocut --Data-MS=big-mslist.txt --Deconv-PeakFactor 0.001000 --Data-ColName DATA --Parallel-NCPU=%i --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=2 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust -0.20000 --Image-NPix=2000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 15.00000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --GAClean-AllowNegativeInitHMP True --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 60.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=4.00 --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=[DDS3_full_smoothed,DDS3_full_slow] --Selection-UVRangeKm=[0.000000,7.0] --GAClean-MinSizeInit=10 --Beam-Smooth=1 --Debug-Pdb=never' % getcpus())
    vlowmask = make_mask('image_full_vlow_nocut.app.restored.fits',3.0)
    run('DDF.py --Output-Name=image_full_vlow_nocut_m --Data-MS=big-mslist.txt --Deconv-PeakFactor 0.001000 --Data-ColName DATA --Parallel-NCPU=%i --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=2 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust -0.20000 --Image-NPix=2000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 15.00000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --GAClean-AllowNegativeInitHMP True --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 60.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=3.00 --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=[DDS3_full_smoothed,DDS3_full_slow] --Selection-UVRangeKm=[0.000000,7.0] --GAClean-MinSizeInit=10 --Beam-Smooth=1 --Debug-Pdb=never --Predict-InitDicoModel=image_full_vlow_nocut.DicoModel --Mask-External=%s'% (getcpus(),vlowmask))
    update_status(None,'Complete')

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Reimage DR2 field at very low resolution')
    parser.add_argument('-f','--field', help='DR2 pointing name', required=True, type=str)
    args = vars(parser.parse_args())
    
    field = args['field']
    
    # Download the appropriate data
    directory=do_download(field,basedir='/beegfs/car/mjh/vlow')

    image_vlow(directory)
    
