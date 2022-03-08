#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from surveys_db import update_reprocessing_extract, get_next_extraction, SurveysDB
import sys
import os
import glob
from astropy.coordinates import SkyCoord
from astropy import units as u
import time
from subprocess import call
from rclone import RClone
from sdr_wrapper import SDR
from fixsymlinks import fixsymlinks
from make_mslists import make_list

def do_rclone_extract_upload(cname,basedir,f,directory):
    '''
    Upload extract results
    '''
    rc=RClone('maca_sksp_disk_extract.conf',debug=True)
    rc.get_remote()
    print(rc.remote,'maca_sksp_disk_extract.conf')
    rc.multicopy(basedir,f,rc.remote+directory+'/'+cname)


def do_rclone_disk_upload(cname,basedir,f,directory):
    '''
    Upload results to surf disk
    '''
    rc=RClone('maca_sksp_disk_subtract.conf',debug=True)
    rc.get_remote()
    print(rc.remote,'maca_sksp_disk_subtract.conf')
    rc.multicopy(basedir,f,rc.remote+directory+'/'+cname)

def do_rclone_tape_pol_upload(cname,basedir,f,directory):
    '''
    Upload results to surf tape
    '''
    rc=RClone('maca_sksp_tape_pol.conf',debug=True)
    rc.get_remote()
    print(rc.remote,'maca_sksp_tape_pol.conf')
    rc.multicopy(basedir,f,rc.remote+directory+'/'+cname)
    
def untar(f,tarfiles,verbose=False):
    print('Untarring files')
    for t in tarfiles:
        if verbose:
            print(t)
        d=os.system('cd %s; tar xf %s; rm %s' % (f,t,t))
        if d!=0:
            raise RuntimeError('untar %s failed!' % t)

def do_sdr_and_rclone_download(cname,f,verbose=False):
    if not os.path.isdir(f):
        os.makedirs(f)
    s=SDR(target=f)
    try:
        status=s.get_status(cname)
    except RuntimeError:
        status=None
    if status:
        if verbose: print('Initiating SDR download for field',cname)
        tarfiles=['images.tar','uv.tar']
        s.download_and_stage(cname,tarfiles)
        untar(f,tarfiles,verbose=verbose)
    else:
        if verbose: print('Trying rclone download for field',cname)
        do_rclone_download(cname,f,verbose=verbose)
    
def do_rclone_download(cname,f,verbose=False):
    '''
    Download required data from field cname into location f
    '''
    tarfiles=['images.tar','uv.tar']
    for macaroon, directory in [('maca_sksp_tape_DR2_readonly.conf',''),('maca_sksp_tape_DDF.conf','archive/'),('maca_sksp_tape_DDF_readonly.conf','other/')]:
        try:
            rc=RClone(macaroon,debug=True)
        except RuntimeError:
            print('Macaroon',macaroon,'does not exist!')
            continue
        rc.get_remote()
        d=rc.multicopy(rc.remote+directory+cname,tarfiles,f)
        if d['err'] or d['code']!=0:
            continue
        break
        
    else:
        raise RuntimeError('Failed to download from any source')
    untar(f,tarfiles,verbose=verbose)
    

def striparchivename():
  mslist = glob.glob('L*_SB*.ms.archive')
  for ms in mslist:
      outname = ms.rstrip('.archive')
      if os.path.exists(outname):
          if os.path.islink(outname):
              print ('Link to',outname,'already exists')
              continue
          else:
              raise RuntimeError(ms+' and '+outname+' both exist in the directory!')
      cmd = 'ln -s ' + ms + ' ' + outname
      print (cmd)
      os.system(cmd)

  return

def prepare_field(field,processingdir):

    cdir = os.getcwd()
    os.chdir(processingdir)

    do_sdr_and_rclone_download(field,processingdir)

    striparchivename()
    fixsymlinks('DDS3_full')
    success=make_list(workdir=os.getcwd())
  
    os.chdir(cdir)

    return
