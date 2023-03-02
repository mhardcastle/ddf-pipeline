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

# Simulate a model image from DR2 solutions

def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]   
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]  
    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

def get_solutions_timerange(sols):
    t = np.load(sols)['BeamTimes']
    return np.min(t),np.max(t)

def fixsymlinks():
    # Code from Tim for fixing symbolic links for DDS3_
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


def predict_model(fakemap):
    os.system('DDF.py --Output-Name=Predict_ADD --Data-MS=big-mslist.txt --Deconv-PeakFactor 0.010000 --Data-ColName DATA --Parallel-NCPU=32 --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=1 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-PhasedArrayMode=A --Weight-Robust -0.150000 --Image-NPix=20000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 1.500000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --GAClean-RMSFactorInitHMP 1.000000 --GAClean-MaxMinorIterInitHMP 10000.000000 --DDESolutions-SolsDir=SOLSDIR --Cache-Weight=reset --Output-Mode=Predict --Predict-ColName ADDED_SOURCE --Output-RestoringBeam 12.000000 --Weight-ColName="None" --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=5.00 --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=[DDS3_full_smoothed,DDS3_full_slow] --Predict-InitDicoModel=None --Selection-UVRangeKm=[0.100000,1000.000000] --GAClean-MinSizeInit=10 --Predict-FromImage=%s'%(fakemap))
    return

def add_predicted(msfiles):
    for msfile in msfiles:
        print 'Creating simulated column for ',msfile
        # Add noise to the predict_modeled column
        t = pt.table(msfile,readonly=False)
        d=t.getcol("DATA")
        predict=t.getcol("ADDED_SOURCE")

        # Append new column containing all sources
        desc = t.getcoldesc('DATA')
        desc['name']='SIMULATED'
        t.addcols(desc)
        t.putcol('SIMULATED',predict+d)
        t.close()
    return

def convolve_model(model,template,outname):
    hdu_template = fits.open(template)
    hdu_model = fits.open(model)
    convfwhm = hdu_template[0].header['BMAJ']/hdu_template[0].header['CDELT1']
    minindex,maxindex = int(hdu_template[0].header['NAXIS1']*0.1),int(hdu_template[0].header['NAXIS1']*0.9)
    print 'Starting model convolution'
    kernel = makeGaussian(size=35,fwhm=convfwhm,center=None)
    print 'Made convolution kernel'
    #hdu_model[0].data[0,0,6000:14000,6000:14000] = convolve2d(hdu_model[0].data[0,0,6000:14000,6000:14000],kernel,mode='same',boundary='fill')
    hdu_model[0].data[0,0,minindex:maxindex,minindex:maxindex] = convolve2d(hdu_model[0].data[0,0,minindex:maxindex,minindex:maxindex],kernel,mode='same',boundary='fill')
    #hdu_model[0].data[0,0,:,:] = convolve2d(hdu_model[0].data[0,0,:,:],kernel,mode='same')#,boundary='fill')
    hdu_model.writeto(outname)
    return outname


def add_to_map(inmap1,inmap2,outname):

    map1 = fits.open(inmap1)
    map2 = fits.open(inmap2)
    data1= map1[0].data
    data2 = map2[0].data
    newdata = data1 + data2
    map1[0].data = newdata
    map1.writeto(outname)


    
if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Add corrupted model to DR2 data')
    parser.add_argument('-f','--field', help='DR2 pointing name', required=True, type=str)
    parser.add_argument('-m','--model', help='Model to inject (full path)', required=True, type=str)
    parser.add_argument('-r','--realmap',help='Real map of the data (sizes matching the model',required=True,type=str)

    args = vars(parser.parse_args())
    
    print args

    field = args['field']
    fakemap = args['model']
    realmap = args['realmap']
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


    # Make some comparison images by adding in the image plane
    convolved_model = convolve_model(fakemap,realmap,fakemap.replace('.fits','-convolved.fits'))
    add_to_map(convolved_model,realmap,convolved_model.replace('.fits','-plusimage.fits'))

    # Predict the input model and add it to the visibilities
    predict_model(fakemap)
    add_predicted(msfiles)

    print 'FINISHED '
    print 'Now run LoTSS-DR2 pipeline on the SIMULATED data column'
