#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
from do_polcubes import do_polcubes
from parset import option_list
from options import options,print_options
from make_mslists import make_list
from auxcodes import run,warn,report
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
import threading
from surveys_db import use_database,SurveysDB,get_id
from upload import do_upload_vlow

def die(error,cname=None):
    update_status(cname,'Failed')
    raise RuntimeError(error)

def do_rsync_download(cname,basedir,f):
    workdir=basedir+'/'+cname

    #if os.environ['DDF_PIPELINE_CLUSTER']!='paracluster':
    target=os.environ['DDF_PIPELINE_LEIDENUSER']+'@ssh.strw.leidenuniv.nl:'
    #else:
    #    target=''

    while True:
	excludeinclude = ' --include="image_full_ampphase_di_m.NS.mask01.fits" --include="image_full_ampphase_di_m.NS.app.restored.fits" --exclude="*QU_*" --exclude="*fits*" --exclude="*.tgz*" --exclude="*QU_*" --exclude="*DDS0*" --exclude="*DDS1*" --exclude="*DDS2*" --exclude="*.corrupted" --exclude="*.pickle" --exclude="*.DicoModel" --exclude="*.txt" --exclude="*.png" '
        s= 'rsync -azvh --timeout=20 --progress --perms --chmod=a+rwx'+ excludeinclude + target+workdir + ' ' + f
        print('Running command:',s)
        retval=call(s,shell=True)
        if retval==0:
            break
        print('Non-zero return value',retval)
        if retval!=30:
            die('rsync failed unexpectedly',cname)
        sleep(10)
        
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
      idd['vlow_reprocess']=status
      sdb.set_field(idd)
        
def striparchivename():
  mslist = glob.glob('L*_SB*.ms.archive')
  for ms in mslist:
      outname = ms.rstrip('.archive')
      cmd = 'ln -s ' + ms + ' ' + outname
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
      die('image_dirin_SSD_m.npy.ClusterCat.npy does not exist')   
  if not os.path.isdir('SOLSDIR'):
      die('SOLSDIR directory does not exist')

  solsfiletmp = glob.glob('DDS3_full*smoothed.npz')
  if len(solsfiletmp) < 1:
      die('Cannot find the DDS3_full*smoothed.npz file(s)')

  solsfiletmp = glob.glob('DDS3_full_slow*.npz')
  if len(solsfiletmp) < 1:
      die('Cannot find the DDS3_full_slow*.npz file(s)')

  return

def compress_fits(filename,q):
    command='fpack -q %i %s' % (q,filename)
    run(command)

def do_download(cname, basedir='.'):
    update_status(cname,'Downloading')
    os.chdir(basedir)
    do_rsync_download(cname,'/disks/paradata/shimwell/LoTSS-DR2/archive/',os.getcwd())
    os.chdir(cname)
    striparchivename()

    success=make_list(workdir=os.getcwd())
    if not success:
        update_status(cname,'Download failed')
        raise RuntimeError('Failed to make mslist')
    
    filechecker()
    fixsymlinks()
    os.system('cp '+os.environ['DDF_DIR']+'/ddf-pipeline/examples/tier1-jul2018.cfg reprocess-vlow.cfg')
    update_status(cname,'Downloaded')
    return os.getcwd() # return directory where everything has been done

def run_reprocess(wd=None):
    # by default assume we're in the working directory at this point
    if wd is not None:
        os.chdir(wd)
    update_status(None,'Running')
    solsfile = glob.glob('DDS3_full*smoothed.npz')
    if len(solsfile) < 1:
        die('Cannot find the correct solution file -- exiting')
    solsfile = str(solsfile[0])
    o = options('reprocess-vlow.cfg',option_list)
    cubefiles=['image_full_vlow_QU.cube.dirty.fits','image_full_vlow_QU.cube.dirty.corr.fits']
    cthreads=[]
    flist=[]
    ddf_kw = {}
    do_polcubes('DATA','[DDS3_full_smoothed,DDS3_full_slow]',[o['image_uvmin'],1.600000],'image_full_vlow',ddf_kw,beamsize=o['vlow_psf_arcsec'],imsize=o['vlow_imsize'],cellsize=o['vlow_cell'],robust=o['vlow_robust'],options=o,catcher=None)
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
            if thread.isAlive():
                warn('Waiting for a compression thread to finish')
                thread.join()
        if o['delete_compressed']:
            for f in flist:
                warn('Deleting compressed file %s' % f)
                os.remove(f)
    update_status(None,'Complete')

    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='Redo vlow resolution cubes')
    parser.add_argument('--field', help='field name', default=None, type=str)
    args = vars(parser.parse_args())
    cname = args['field']

    if cname == None:
      print('No field specified -- exiting')
      sys.exit(0)

    do_download(cname)
    run_reprocess()
    os.chdir('..')
    do_upload_vlow(cname,'.')
