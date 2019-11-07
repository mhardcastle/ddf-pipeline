from do_polcubes import do_polcubes
from parset import option_list
from options import options,print_options
import casacore.tables as pt
import os,sys
import numpy as np
import argparse
import pyregion
from astropy.io import fits
from astropy.wcs import WCS
from astropy.io import ascii
import glob
from subprocess import call
import argparse

def do_rsync_download(cname,basedir,f):
    workdir=basedir+'/'+cname

    #if os.environ['DDF_PIPELINE_CLUSTER']!='paracluster':
    target='lofararchive@ssh.strw.leidenuniv.nl:'
    #else:
    #    target=''

    while True:
	excludeinclude = ' --include="image_full_ampphase_di_m.NS.mask01.fits" --include="image_full_ampphase_di_m.NS.app.restored.fits" --exclude="*QU_*" --exclude="*fits*" --exclude="*.tgz*" --exclude="*QU_*" --exclude="*DDS0*" --exclude="*DDS1*" --exclude="*DDS2*" --exclude="*.corrupted" '
        s= 'rsync -azvh --timeout=20 --progress --perms --chmod=a+rwx'+ excludeinclude + target+workdir + ' ' + f
        #'cd '+workdir+'; rsync -avz --progress --safe-links --inplace --append --partial --timeout=20 '+' '.join(f)+' '+target+'/disks/paradata/shimwell/LoTSS-DR2/archive/'+name
        print 'Running command:',s
        retval=call(s,shell=True)
        if retval==0:
            break
        print 'Non-zero return value',retval
        if retval!=30:
            raise RuntimeError('rsync failed unexpectedly')
        sleep(10)

def striparchivename():
  mslist = glob.glob('L*_SB*.ms.archive')
  for ms in mslist:
      outname = ms.rstrip('.archive')
      cmd = 'mv ' + ms + ' ' + outname
      print (cmd)
      os.system(cmd)

  return

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

def filechecker():
  '''
  Check if files are present to avoid errors to avoid crashes
  '''
  if not os.path.isfile('image_dirin_SSD_m.npy.ClusterCat.npy'):
    raise IOError('image_dirin_SSD_m.npy.ClusterCat.npy does not exist')   
  if not os.path.isdir('SOLSDIR'):
   raise IOError('SOLSDIR directory does not exist')

  solsfiletmp = glob.glob('DDS3_full*smoothed.npz')
  if len(solsfiletmp) < 1:
     raise IOError('Cannot find the DDS3_full*smoothed.npz file(s)')

  solsfiletmp = glob.glob('DDS3_full_slow*.npz')
  if len(solsfiletmp) < 1:
     raise IOError('Cannot find the DDS3_full_slow*.npz file(s)')

  return

parser = argparse.ArgumentParser(description='Redo vlow resolution cubes')
parser.add_argument('--field', help='field name', default=None, type=str)
args = vars(parser.parse_args())
cname = args['field']
os.system('cp /net/lofar9/data2/shimwell/smoothbeam/tier1-jul2018.cfg .')
if cname == None:
  print 'No field specified -- exiting'
  sys.exit(0)

do_rsync_download(cname,'/disks/paradata/shimwell/LoTSS-DR2/archive/',os.getcwd())
os.chdir(cname)
ddf_kw = {}

striparchivename()

from make_mslists import make_list
success=make_list(workdir=os.getcwd())

filechecker()
fixsymlinks()

solsfile = glob.glob('DDS3_full*smoothed.npz')
if len(solsfile) < 1:
     print 'Cannot find the correct solution file -- exiting'
     sys.exit()
solsfile = str(solsfile[0])

o = options('reprocess-vlow.cfg',option_list)
do_polcubes('DATA','[DDS3_full_smoothed,DDS3_full_slow]',[o['image_uvmin'],1.600000],'image_full_vlow',ddf_kw,beamsize=o['vlow_psf_arcsec'],imsize=o['vlow_imsize'],cellsize=o['vlow_cell'],robust=o['vlow_robust'],options=o,catcher=None)

