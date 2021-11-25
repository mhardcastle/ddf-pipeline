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

def do_rclone_upload(cname,basedir,f,directory):
    '''
    Upload extract results
    '''
    rc=RClone('maca_sksp_disk_extract.conf',debug=True)
    rc.get_remote()
    print(rc.remote,'maca_sksp_disk_extract.conf')
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
    for macaroon, directory in [('maca_sksp_tape_DR2_readonly.conf',''),('maca_sksp_tape_DDF.conf','archive/'),('maca_sksp_tape_DDF.conf','other/')]:
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
        field = sys.argv[2]
        with SurveysDB(readonly=True) as sdb:
            sdb.cur.execute('select * from reprocessing where id="%s"' % target)
            results=sdb.cur.fetchall()
        if len(results)==0:
            raise RuntimeError('Requested target is not in database')

        fields = results[0]['fields'].split(',')
        if field not in fields:
            raise RuntimeError('Requested field is not in target list')
        bad_pointings = results[0]['bad_pointings']
        if bad_pointings is None:
            bad_pointings = ['']
        else:
            bad_pointings = bad_pointings.split(',')
        if field in bad_pointings:
            raise RuntimeError('Field is in bad pointing list')
        ra=results[0]['ra']
        dec=results[0]['decl']
        size=results[0]['size']
            

    startdir = os.getcwd()
    os.system('mkdir %s'%target)
    os.chdir(target)  

    update_reprocessing_extract(target,field,'STARTED')

    do_sdr_and_rclone_download(field,startdir+'/'+target + '/'+field)

    os.chdir(field)
    create_ds9_region('%s.ds9.reg'%(target),ra,dec,size)

    executionstr = 'sub-sources-outside-region.py -b %s.ds9.reg -p %s'%(target,target)
    print(executionstr)
    result=os.system(executionstr)
    if result!=0:
        raise RuntimeError('Failed to run sub-sources')

    resultfiles = glob.glob('*archive*')
    resultfilestar = []
    for resultfile in resultfiles:
        d=os.system('tar -cvf %s.tar %s'%(resultfile,resultfile))
        if d!=0:
            raise RuntimeError('Tar of %s failed'%resultfile)	
        resultfilestar.append('%s.tar'%resultfile)

    do_rclone_upload(field,os.getcwd(),resultfilestar,target)

    update_reprocessing_extract(target,field,'EDONE')

