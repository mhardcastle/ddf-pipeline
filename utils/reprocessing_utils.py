#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from surveys_db import SurveysDB
from db_utils import update_reprocessing_extract, get_next_extraction
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
from parset import option_list
from options import options,print_options
from time import sleep

def convert_summary_cfg(option_list,summaryname='summary.txt',newconfigname='tier1-reprocessing.cfg'):
    summary = open(summaryname,'r')
    summary_dict = {}
    print('Summary of the reprocessing pipeline inputs from summary.txt')
    for line in summary:
        line = line.replace(" ", "")
        line = line.strip().split(':')
        if len(line) != 2:
            continue
        summary_dict[line[0]] = line[1]
        
    newconfig = open(newconfigname,'w')
    option_types = []
    for option in option_list:
        if option[0] not in option_types:
            option_types.append(option[0])


    for option_type in option_types:
        newconfig.write('[%s]\n'%option_type)
        for option in option_list:
            if option[0] != option_type:
                continue
            try:
                summary_value = summary_dict[option[1]]
            except KeyError:
                try:
                    summary_value = summary_dict[option[0]+'_'+option[1]]
                except KeyError:
                    print('Cannot find option for parameter',option[0],option[1])
                    continue
            if summary_value == 'None':
                continue
            print(option,summary_value)
            newconfig.write('%s=%s\n'%(option[1],summary_value))
    newconfig.close()

    o = options(newconfigname,option_list)
    return o

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
    rc=RClone('maca_sksp_disk_reproc.conf',debug=True)
    rc.get_remote()
    print(rc.remote,'maca_sksp_disk_reproc.conf')
    return rc.multicopy(basedir,f,rc.remote+directory+'/'+cname)

def do_rclone_reproc_tape_upload(cname,basedir,f,directory):
    '''
    Upload results to surf disk
    '''
    rc=RClone('maca_sksp_tape_reproc.conf',debug=True)
    rc.get_remote()
    print(rc.remote,'maca_sksp_tape_reproc.conf')
    return rc.multicopy(basedir,f,rc.remote+directory+'/'+cname)


def do_rclone_reproc_tape_download(cname,f,directory,verbose=False):
    '''
    Download required data from field cname into location f
    '''
    rc=RClone('maca_sksp_tape_reproc.conf',debug=True)
    rc.get_remote()
    files=rc.get_files(directory+'/' + cname)
    if verbose: print(files)
    tarfiles=None
    if files:
        tarfiles = [fl for fl in files]
        
    if tarfiles is not None:
        d=rc.multicopy(rc.remote+directory+'/' + cname,tarfiles,f)

    untar(f,tarfiles,verbose=verbose)

    
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

def stage_field_SDR(cname,f,verbose=False,Mode='Imaging+Misc',timeout=None,sleeptime=30):
    timer=0
    s=SDR(target=f)
    # Next line raises RuntimeError if files don't exist
    files=s.get_status(cname)
    if verbose: print('Initiating SDR stage for field',cname)
    if Mode=="Imaging":
        tarfiles=['images.tar','uv.tar']
    elif Mode=="Misc":
        tarfiles=['misc.tar']
    elif Mode=="Imaging+Misc":
        tarfiles=['images.tar','uv.tar','misc.tar',"stokes_small.tar"]

    for f in tarfiles:
        if f not in files:
            raise RuntimeError(f'File {f} not found!')
        else:
            if files[f]=='OFL':
                s.stage(cname,f)

    if verbose:
        print('Waiting for files to be online:')

    while True:
        if timeout is not None and timer>timeout:
            raise RuntimeError('Timer exceeded')
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
            sleep(sleeptime)
            timer+=sleeptime

def stage_field_rclone(cname,f,verbose=False,Mode='Imaging+Misc',timeout=None,sleeptime=30,
                       macaroon='maca_sksp_tape_DDF.conf', directory='archive/'):
    timer=0
    rc=RClone(macaroon,debug=True)
    # Next line raises RuntimeError if directory does not exist
    files=rc.get_files(directory+cname)
    if verbose: print('Initiating rclone/ADA stage for field',cname)
    if Mode=="Imaging":
        tarfiles=['images.tar']+[f for f in files if 'uv' in f]
    elif Mode=="Misc":
        tarfiles=['misc.tar']
    elif Mode=="Imaging+Misc":
        tarfiles=['images.tar','misc.tar',"stokes_small.tar"]+[f for f in files if 'uv' in f]
    to_stage=[]
    for f in tarfiles:
        staged=rc.check_stage(directory+cname+'/'+f)
        if 'ONLINE' not in staged:
            to_stage.append(f)
            rc.stage(directory+cname+'/'+f)

    if verbose:
        print('Waiting for files to be online:')
        while True:
            if timeout is not None and timer>timeout:
                raise RuntimeError('Timer exceeded')
            count=0
            for f in to_stage:
                staged=rc.check_stage(directory+cname+'/'+f)
                if 'ONLINE' in staged:
                    count+=1
            if verbose:
                print('%i/%i... ' % (count,len(to_stage)),end='')
                sys.stdout.flush()
            if count==len(to_stage):
                if verbose: print()
                break
            else:
                sleep(sleeptime)
                timer+=sleeptime
        
def stage_field(cname,f,verbose=False,Mode='Imaging+Misc',timeout=None,sleeptime=30,order=['rclone','SDR']):
    ''' Wrapper around staging for rclone and SDR. Try them in the specified order '''
    for step in order:
        if verbose:
            print('Stage: Trying step',step)
        stage={'rclone':stage_field_rclone,'SDR':stage_field_SDR}[step]
        try:
            stage(cname,f,verbose=verbose,Mode=Mode,timeout=timeout,sleeptime=sleeptime)
            success=True
        except RuntimeError as e:
            print('Error',e,'caught')
            success=False
        if success:
            break
    else:
        raise RuntimeError('Unable to stage from any source')

def do_sdr_and_rclone_download(cname,f,verbose=False,Mode="Imaging+Misc",operations=['download','untar'],order=['rclone','SDR']):
    ''' download tar files from field cname to location f and optionally untar. Use rclone and SDR in an order specified by the user '''
    if not os.path.isdir(f):
        os.makedirs(f)

    for step in order:
        if verbose:
            print('Download: Trying step',step)
        download={'rclone':do_rclone_download,'SDR':do_sdr_download}[step]
        try:
            download(cname,f,verbose=verbose,Mode=Mode,operations=operations)
            success=True
        except RuntimeError as e:
            print('Error',e,'caught')
            success=False
        if success:
            break
    else:
        raise RuntimeError('Unable to download from any source')

def do_sdr_download(cname,f,verbose=False,Mode="Imaging+Misc",operations=['download','untar']):
    s=SDR(target=f)
    # next line will raise a RuntimeError if directory is not present
    try:
        status=s.get_status(cname)
    except:
        raise RuntimeError('Failed to find SDR data')
    if Mode=="Imaging":
        tarfiles=['images.tar','uv.tar']
    elif Mode=="Misc":
        tarfiles=['misc.tar']
    elif Mode=="Imaging+Misc":
        tarfiles=['images.tar','uv.tar','misc.tar',"stokes_small.tar"]
    elif Mode=="ImageOnly":
        tarfiles=['images.tar']
    else:
        raise NotImplementedError('Unknown mode '+Mode+' requested')

    if 'download' in operations:
        if verbose: print('Initiating SDR download for field',cname)
        s.download_and_stage(cname,tarfiles,progress_bar=verbose)
    if 'untar' in operations:
        #tarfiles = glob.glob('*tar')
        untar(f,tarfiles,verbose=verbose)
        
def do_rclone_download(cname,f,verbose=False,Mode="Imaging+Misc",operations=['download','untar']):
    '''
    Download required data from field cname into location f
    '''
    #tarfiles=['images.tar','uv.tar']
    for macaroon, directory in [('maca_sksp_tape_DDF.conf','archive/'),('maca_sksp_tape_DDF_readonly.conf','other/'),('maca_sksp_tape_DDF_readonly.conf','archive/')]:
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
            else:
                raise NotImplementedError('Unknown mode requested')
            
        if 'download' in operations and tarfiles is not None:
            d=rc.multicopy(rc.remote+directory+cname,tarfiles,f)
            if d['err'] or d['code']!=0:
                continue # try next source

        if tarfiles is not None:
            break # out of loop to unpack
        
    else:
        raise RuntimeError('Failed to find rclone data')
    
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

    from make_mslists import make_list

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
