#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from auxcodes import report,warn,die
from surveys_db import *
from download import download_dataset
from download_field import download_field

from auxcodes import MSList
import sys
import os
import glob
from astropy.coordinates import SkyCoord
from astropy import units as u
import time
from subprocess import call
from rclone import RClone

def do_rclone_upload(cname,basedir,f):
    '''
    Untested upload code
    '''
    rc=Rclone('add-your-macaroon-here.conf',debug=True)
    rc.multicopy(basedir,f,rc.remote+'destination-directory/'+cname)

def do_rsync_upload(cname,basedir,f):
    workdir=basedir+'/'+cname

    #if os.environ['DDF_PIPELINE_CLUSTER']!='paracluster':
    target='lofararchive@ssh.strw.leidenuniv.nl:'
    #else:
    #    target=''

    while True:
        s= 'rsync -avz --relative --progress --perms --chmod=ugo+rX --safe-links --partial --timeout=20 '+' '.join(f)+' '+target+'/disks/paradata/shimwell/LoTSS-DR2/archive_extract/'+cname 
        print('Running command:',s)
        retval=call(s,shell=True)
        if retval==0:
            break
        print('Non-zero return value',retval)
        if retval!=30:
            raise RuntimeError('rsync failed unexpectedly')
        time.sleep(10)

def do_rclone_download(cname,f):
    '''
    Download required data from field cname into location f
    '''
    tarfiles=['images.tar','uv.tar']
    for macaroon, directory in [('maca_sksp_tape_DR2_readonly.conf',''),('maca_sksp_tape_DDF.conf','archive/')]:
        try:
            rc=RClone(macaroon,debug=True)
        except RuntimeError:
            print('Macaroon',m,'does not exist!')
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
    print('Untarring images')
    for t in tarfiles:
        d=os.system('cd %s; tar xf %s' % (f,t))
        if d!=0:
            raise RuntimeError('untar %s failed!' % t)

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
        print('Running command:',s)
        retval=call(s,shell=True)
        if retval==0:
            break
        print('Non-zero return value',retval)
        if retval!=30:
            raise RuntimeError('rsync failed unexpectedly')
        time.sleep(10)

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


def do_run_subtract(name,basedir,inarchivedir,outarchivedir,force=False):
    startdir = os.getcwd()
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    sdb.close()
    fields = extractdict['fields'].split(',')
    extract_status = extractdict['extract_status'].split(',')
    try:
        bad_pointings = extractdict['bad_pointings'].split(',')
    except AttributeError:
        bad_pointings = ['']
    print('Working on ',name, 'in fields', fields,'which have status',extract_status)
    
    for i in range(0,len(fields)):
        os.chdir(startdir)
        if not(extract_status[i] == 'EREADY' or (force and extract_status[i] == 'STARTED')):
            continue
        field = fields[i]
        if field in bad_pointings:
            print('Field',field,'in bad pointings -- skipping and setting to BADP')
            sdb=SurveysDB()
            extractdict = sdb.get_reprocessing(name)
            extract_status[i] = 'BADP'
            extractdict['extract_status'] = ','.join(extract_status)
            sdb.db_set('reprocessing',extractdict)
            sdb.close()
            continue
        workdir=basedir+'/'+name
        try:
            os.mkdir(workdir)
        except OSError:
            warn('Working directory already exists')
        print('In directory', os.getcwd())
        os.chdir(workdir)
        # Update status to running here
        extract_status[i] = 'STARTED'
        sdb=SurveysDB()
        extractdict = sdb.get_reprocessing(name)
        extractdict['extract_status'] = ','.join(extract_status)
        sdb.db_set('reprocessing',extractdict)
        sdb.close()
        print('Updated status to STARTED for',field,name)
        time.sleep(2.0)
        report('Copying data from %s'%inarchivedir)
        
        # WANT TO MAKE THIS INTO A RSYNC SO THAT IT CAN BE DONE OUTSIDE LEIDEN
        #os.system('cp -r %s/%s %s'%(inarchivedir,field,workdir))
        do_rsync_download(field,inarchivedir,workdir)

        # Update status to copied here
        extract_status[i] = 'COPIED'
        sdb=SurveysDB()
        extractdict = sdb.get_reprocessing(name)
        extractdict['extract_status'] = ','.join(extract_status)
        sdb.db_set('reprocessing',extractdict)
        sdb.close()
        print('Updated status to COPIED for',field,name)


        # Create boxfile
        create_ds9_region('%s.ds9.reg'%name,extractdict['ra'],extractdict['decl'],extractdict['size'])


        # Run subtract code
        print(os.getcwd(), 'working here')
        os.chdir(field)
        print ('sub-sources-outside-region.py -b %s/%s.ds9.reg -p %s'%(workdir,name,name))
        result=os.system('sub-sources-outside-region.py -b %s/%s.ds9.reg -p %s'%(workdir,name,name))
        if result!=0:
            raise RuntimeError('sub-sources-outside-region.py failed with error code %i' % result)
        
        # Archive the results need an rsync code this is just the *archive file that needs to be archived.
        #os.system('mkdir %s/%s'%(outarchivedir,name))
        #os.system('mkdir %s/%s/%s'%(outarchivedir,name,field))
        os.chdir(workdir)
        f = glob.glob('%s/*.archive*'%(field))
        do_rsync_upload(name,field,f)

        #print  ('cp -r %s_%s.dysco.sub.shift.avg.weights.ms.archive %s/%s/%s'%(field,name,outarchivedir,name,field))
        #os.system('cp -r %s_%s.dysco.sub.shift.avg.weights.ms.archive %s/%s/%s'%(field,name,outarchivedir,name,field))


        # update the database to give success
        extract_status[i] = 'EDONE'
        sdb=SurveysDB()
        extractdict = sdb.get_reprocessing(name)
        extractdict['extract_status'] = ','.join(extract_status)
        sdb.db_set('reprocessing',extractdict)
        sdb.close()
        print('Updated status to EDONE for',field,name)

    # update the database to give selfcal status as SREADY
    selfcal_status = 'SREADY'
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    extractdict['selfcal_status'] = selfcal_status
    sdb.db_set('reprocessing',extractdict)
    sdb.close()
    print('Updated status to SREADY for',name)

if __name__=='__main__':

    if len(sys.argv)==1:
        target = get_next_extraction()['id']
        force = False
    else:
        force = True
        target = sys.argv[1]

    # Takes the targetname, the current directory (the working directory), and the directory that contains the LoTSS-DR2 archive
    do_run_subtract(target,os.getcwd(),'/disks/paradata/shimwell/LoTSS-DR2/archive/','/disks/paradata/shimwell/LoTSS-DR2/archive_extract/',force=force)

