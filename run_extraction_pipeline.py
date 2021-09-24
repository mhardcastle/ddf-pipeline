#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from surveys_db import *
import sys
import os
import glob
from astropy.coordinates import SkyCoord
from astropy import units as u
import time
from subprocess import call
from rclone import RClone

macaroonpath = '/project/lotss/Software/ddf-operations/testing/'

def do_rclone_upload(cname,basedir,f,directory):
    '''
    Untested upload code
    '''
    rc=RClone('%s/maca_sksp_disk_extract.conf'%macaroonpath,debug=True)
    rc.get_remote()
    print(rc.remote,'%s/maca_sksp_disk_extract.conf'%macaroonpath)
    rc.multicopy(basedir,f,rc.remote+directory+'/'+cname)

def do_rclone_download(cname,f):
    '''
    Download required data from field cname into location f
    '''
    tarfiles=['images.tar','uv.tar']
    for macaroon, directory in [('%s/maca_sksp_tape_DR2_readonly.conf'%macaroonpath,''),('%s/maca_sksp_tape_DDF.conf'%macaroonpath,'archive/'),('%s/maca_sksp_tape_DDF.conf'%macaroonpath,'other/')]:
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
        
        '''
        for t in tarfiles:
            print('Downloading',t)
            d=rc.copy(rc.remote+directory+cname+'/'+t,f)
            if d['err'] or d['code']!=0:
                break # failed so break out of inner loop
        else:
            break # succeeded so break out of outer loop
        '''
    else:
        raise RuntimeError('Failed to download from any source')
    print('Untarring files')
    for t in tarfiles:
        d=os.system('cd %s; tar xf %s; rm %s' % (f,t,t))
        if d!=0:
            raise RuntimeError('untar %s failed!' % t)

def create_ds9_region(filename,ra,dec,size):

    sc=SkyCoord('%sdeg'%ra,'%sdeg'%dec,frame='icrs')
    scstring = sc.to_string(style='hmsdms',sep=':',precision=2)
    openfile = open(filename,'w')
    openfile.write('# Region file format: DS9 version 4.1\n')
    openfile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    openfile.write('fk5\n')
    openfile.write('box(%s,%s",%s",0.000136627)\n'%(scstring.replace(' ',','),size*60.0*60.0,size*60.0*60.0))

    openfile.close()
    return(filename)


def stage_field_data(field):

    for macaroon, directory in [('%s/maca_sksp_tape_DR2_readonly.conf'%macaroonpath,''),('%s/maca_sksp_tape_DDF.conf'%macaroonpath,'archive/'),('%s/maca_sksp_tape_DDF.conf'%macaroonpath,'other/')]:
         print('ada --tokenfile %s --stage %s/%s'%(macaroon,directory,field))
         os.system('ada --tokenfile %s --stage %s/%s'%(macaroon,directory,field))
    return

def unstage_field_data(field):

    for macaroon, directory in [('%s/maca_sksp_tape_DR2_readonly.conf'%macaroonpath,''),('%s/maca_sksp_tape_DDF.conf'%macaroonpath,'archive/'),('%s/maca_sksp_tape_DDF.conf'%macaroonpath,'other/')]:
         print('ada --tokenfile %s --unstage %s/%s'%(macaroon,directory,field))
         os.system('ada --tokenfile %s --unstage %s/%s'%(macaroon,directory,field))
    return

def check_staged(field):

   notstaged = True
   while notstaged:
       for macaroon, directory in [('%s/maca_sksp_tape_DR2_readonly.conf'%macaroonpath,''),('%s/maca_sksp_tape_DDF.conf'%macaroonpath,'archive/'),('%s/maca_sksp_tape_DDF.conf'%macaroonpath,'other/')]:
           print('ada --tokenfile %s --longlist %s/%s'%(macaroon,directory,field))
           stageinfo = os.popen('ada --tokenfile %s --longlist %s/%s'%(macaroon,directory,field)).read().split()
           numofnotstaged = stageinfo.count('NEARLINE')
           if numofnotstaged == 0:
              notstage = False
              return
           else:
              print(stageinfo)
              print('%s not staged'%numofnotstaged)
              time.sleep(600)


def do_run_extract(field,name):

        # Run subtract code
        os.chdir(startdir)
        os.chdir(name)
        os.chdir(field)
        print ('sub-sources-outside-region.py -b %s/%s.ds9.reg -p %s'%(startdir,name,name,name))
        result=os.system('sub-sources-outside-region.py -b %s.ds9.reg -p %s'%(workdir,name,name))
        if result!=0:
            raise RuntimeError('sub-sources-outside-region.py failed with error code %i' % result)

if __name__=='__main__':

    if len(sys.argv)==1:
        target,field,ra,dec,size = get_next_extraction()
        force = False
    else:
        force = True
        target = sys.argv[1]

    startdir = os.getcwd()
    os.system('mkdir %s'%target)
    os.chdir(target)  

    # Stage and wait for coming online
    stage_field_data(field)
    check_staged(field)

    update_reprocessing_extract(target,field,'STARTED')

    do_rclone_download(field,startdir+'/'+target + '/'+field)
    unstage_field_data(field)
    os.chdir(field)
    create_ds9_region('%s.ds9.reg'%(target),ra,dec,size)

    singularityfile = '/project/lotss/Software/lofar_sksp_fedora31_ddf.sif'

    executionstr = 'singularity exec -B  /scratch/,/home/lotss-tshimwell/,/project/lotss/ %s sub-sources-outside-region.py -b %s.ds9.reg -p %s'%(singularityfile,target,target)
    print(executionstr)
    os.system(executionstr)

    resultfiles = glob.glob('*archive*')
    resultfilestar = []
    for resultfile in resultfiles:
        d=os.system('tar -cvf %s.tar %s'%(resultfile,resultfile))
        if d!=0:
            raise RuntimeError('Tar of %s failed'%resultfile)	
        resultfilestar.append('%s.tar'%resultfile)

    do_rclone_upload(field,os.getcwd(),resultfilestar,target)

    update_reprocessing_extract(target,field,'EDONE')

