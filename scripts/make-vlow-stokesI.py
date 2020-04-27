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

# Reimage at very low resolution

def make_mask(imagename,thresh):

    fname=imagename+'.mask.fits'
    runcommand = "MakeMask.py --RestoredIm=%s --Th=%s --Box=50,2"%(imagename,thresh)
    os.system(runcommand)
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
      cmd = 'mv ' + ms + ' ' + outname
      print (cmd)
      os.system(cmd)

  return


def do_rsync_download(cname,basedir,f):
    workdir=basedir+'/'+cname

    #if os.environ['DDF_PIPELINE_CLUSTER']!='paracluster':
    target='lofararchive@ssh.strw.leidenuniv.nl:'
    #else:
    #    target=''

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
            raise RuntimeError('rsync failed unexpectedly')
        time.sleep(10)

def image_vlow():
    os.system('DDF.py --Output-Name=image_full_vlow_nocut --Data-MS=big-mslist.txt --Deconv-PeakFactor 0.001000 --Data-ColName DATA --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=2 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust -0.20000 --Image-NPix=3000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 10.00000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --GAClean-AllowNegativeInitHMP True --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 60.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=5.00 --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=[DDS3_full_smoothed,DDS3_full_slow] --Selection-UVRangeKm=[0.000000,6.0] --GAClean-MinSizeInit=10 --Beam-Smooth=1')
    vlowmask = make_mask('image_full_vlow_nocut.app.restored.fits',5.0)
    os.system('DDF.py --Output-Name=image_full_vlow_nocut_m --Data-MS=big-mslist.txt --Deconv-PeakFactor 0.001000 --Data-ColName DATA --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=2 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust -0.20000 --Image-NPix=3000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 10.00000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --GAClean-AllowNegativeInitHMP True --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Clean --Output-RestoringBeam 60.000000 --Weight-ColName="IMAGING_WEIGHT" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=5.00 --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=[DDS3_full_smoothed,DDS3_full_slow] --Selection-UVRangeKm=[0.000000,6.0] --GAClean-MinSizeInit=10 --Beam-Smooth=1 --Predict-InitDicoModel=image_full_vlow_nocut.DicoModel --Mask-External=%s'%vlowmask)


    return

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Add corrupted model to DR2 data')
    parser.add_argument('-f','--field', help='DR2 pointing name', required=True, type=str)
    args = vars(parser.parse_args())
    
    print args


    field = args['field']
    inarchivedir = '/disks/paradata/shimwell/LoTSS-DR2/archive/'
    workdir = os.getcwd()
    
    # Download the appropriate data
    do_rsync_download(field,inarchivedir,workdir)

    # Prepare data
    os.chdir(field)
    striparchivename()
    success=make_list(workdir=os.getcwd())
    fixsymlinks()
    msfiles   = ascii.read('big-mslist.txt',data_start=0)
    msfiles   = list(msfiles[:][msfiles.colnames[0]])

    image_vlow()
    
