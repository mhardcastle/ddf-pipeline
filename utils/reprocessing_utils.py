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
    return rc.multicopy(basedir,f,rc.remote+directory+'/'+cname)

def do_rclone_disk_upload(cname,basedir,f,directory):
    '''
    Upload results to surf disk
    '''
    rc=RClone('maca_sksp_disk_subtract.conf',debug=True)
    rc.get_remote()
    print(rc.remote,'maca_sksp_disk_subtract.conf')
    return rc.multicopy(basedir,f,rc.remote+directory+'/'+cname)
    
def do_rclone_tape_pol_upload(cname,basedir,f,directory):
    '''
    Upload results to surf tape
    '''
    rc=RClone('maca_sksp_tape_pol.conf',debug=True)
    rc.get_remote()
    print(rc.remote,'maca_sksp_tape_pol.conf')
    return rc.multicopy(basedir,f,rc.remote+directory+'/'+cname)
    
def untar(f,tarfiles,verbose=False):
    print('Untarring files')
    for t in tarfiles:
        if verbose:
            print(t)
        d=os.system('cd %s; tar xf %s; rm %s' % (f,t,t))
        if d!=0:
            raise RuntimeError('untar %s failed!' % t)

def do_sdr_and_rclone_download(cname,f,verbose=False,Mode="Imaging+Misc",operations=['download','untar']):
    ''' download tar files from field cname to location f. Try SDR and if that doesn't work try rclone. '''
    if not os.path.isdir(f):
        os.makedirs(f)
    s=SDR(target=f)
    try:
        status=s.get_status(cname)
    except RuntimeError:
        status=None
    if status:
        if Mode=="Imaging":
            tarfiles=['images.tar','uv.tar']
        elif Mode=="Misc":
            tarfiles=['misc.tar']
        elif Mode=="Imaging+Misc":
            tarfiles=['images.tar','uv.tar','misc.tar',"stokes_small.tar"]

        if 'download' in operations:
            if verbose: print('Initiating SDR download for field',cname)
            s.download_and_stage(cname,tarfiles,progress_bar=verbose)
        if 'untar' in operations:
            #tarfiles = glob.glob('*tar')
            untar(f,tarfiles,verbose=verbose)
    else:
        if verbose: print('Trying rclone download for field',cname)
        do_rclone_download(cname,f,verbose=verbose,Mode=Mode,operations=operations)

def do_rclone_download(cname,f,verbose=False,Mode="Imaging+Misc",operations=['download','untar']):
    '''
    Download required data from field cname into location f
    '''
    #tarfiles=['images.tar','uv.tar']
    for macaroon, directory in [('maca_sksp_tape_DR2_readonly.conf',''),('maca_sksp_tape_DDF.conf','archive/'),('maca_sksp_tape_DDF_readonly.conf','other/')]:
        try:
            rc=RClone(macaroon,debug=True)
        except RuntimeError:
            print('Macaroon',macaroon,'does not exist!')
            continue
        rc.get_remote()
        files=rc.get_files(directory+cname)
        print(files)
        tarfiles=None
        if files:
            if Mode=="Imaging":
                tarfiles=[fl for fl in files if 'images' in fl or 'uv' in fl]
            elif Mode=="Misc":
                tarfiles=[fl for fl in files if 'misc.tar'==fl]
            elif Mode=="Imaging+Misc":
                tarfiles=[fl for fl in files if 'images' in fl or 'uv' in fl or 'misc.tar'==fl or "stokes_small.tar"==fl]
            
        if 'download' in operations and tarfiles is not None:
            d=rc.multicopy(rc.remote+directory+cname,tarfiles,f)
            if d['err'] or d['code']!=0:
                continue # try next source

        if tarfiles is not None:
            break # out of loop to unpack
        
    else:
        raise RuntimeError('Failed to download from any source')
    #tarfiles = glob.glob('*tar')
    if 'untar' in operations and tarfiles is not None:
        untar(f,tarfiles,verbose=verbose)
    

def striparchivename(workdir='.'):
  mslist = glob.glob(workdir+'/L*.ms.archive')
  for ms in mslist:
      outname = ms.rstrip('.archive')
      if os.path.exists(outname):
          if os.path.islink(outname):
              print ('Link to',outname,'already exists')
              continue
          else:
              raise RuntimeError(ms+' and '+outname+' both exist in the directory!')
      cmd = 'ln -s ' + ms + ' ' + outname
      print(cmd)
      os.system(cmd)

  return

def prepare_field(field,processingdir,verbose=False,Mode="Imaging+Misc",operations=['download','untar','fixsymlinks','makelist']):

    ''' General function to prepare a field for reprocessing.
    field: the field to download
    processingdir: where we will download and unpack the data (i.e. normally in a directory with the same name as the field
    verbose: set True to have more interactive information about what's being done
    Mode: controls what tarfiles are downloaded. Default is 'Imaging+Misc' which gets all files including misc.tar, needed for some reprocessing because it includes summary.txt'
    operations: list of operations to carry out. By default all of them will be done but a subset can be passed to break up the prepare_field run.
    '''
    
    if not os.path.isdir(processingdir):
        if verbose:
            print('Creating directory',processingdir)
        os.mkdir(processingdir)

    if 'download' in operations or 'untar' in operations:
        if verbose:
            print('Calling sdr/rclone code to download/untar')
        do_sdr_and_rclone_download(field,processingdir,verbose=verbose,Mode=Mode,operations=operations)

    if 'fixsymlinks' in operations:
        if verbose: print('Fixing symlinks')
        striparchivename(workdir=processingdir)
        fixsymlinks('DDS3_full',workdir=processingdir)
    if 'makelist' in operations:
        if verbose: print('Making list')
        success=make_list(workdir=processingdir)
    else:
        success=True
        
    return success
