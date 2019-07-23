#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from auxcodes import report,warn,die
from surveys_db import *
from download import download_dataset
from download_field import download_field
from run_job import do_run_job
from unpack import unpack
from make_mslists import make_list,list_db_update
from average import average
from auxcodes import MSList
import sys
import os
import glob
from astropy.coordinates import SkyCoord
from astropy import units as u
import time
from subprocess import call

def do_rsync_upload(cname,basedir,f):
    workdir=basedir+'/'+cname

    #if os.environ['DDF_PIPELINE_CLUSTER']!='paracluster':
    target='lofararchive@ssh.strw.leidenuniv.nl:'
    #else:
    #    target=''

    while True:
        s= 'rsync -avz --relative --progress --safe-links --inplace --append --partial --timeout=20 '+' '.join(f)+' '+target+'/disks/paradata/shimwell/LoTSS-DR2/archive_extract/'+cname +'/selfcal/' 
        print 'Running command:',s
        retval=call(s,shell=True)
        if retval==0:
            break
        print 'Non-zero return value',retval
        if retval!=30:
            raise RuntimeError('rsync failed unexpectedly')
        sleep(10)

def do_rsync_download(cname,basedir,f):
    workdir=basedir+'/'+cname

    #if os.environ['DDF_PIPELINE_CLUSTER']!='paracluster':
    target='lofararchive@ssh.strw.leidenuniv.nl:'
    #else:
    #    target=''

    while True:
        s= 'rsync -azvh --timeout=20 --progress '+target+workdir + ' ' + f
        print 'Running command:',s
        retval=call(s,shell=True)
        if retval==0:
            break
        print 'Non-zero return value',retval
        if retval!=30:
            raise RuntimeError('rsync failed unexpectedly')
        sleep(10)

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


def do_run_selfcal(name,basedir,inarchivedir,outarchivedir):
    startdir = os.getcwd()
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    sdb.close()
    fields = extractdict['fields'].split(',')
    selfcal_status = extractdict['selfcal_status']

    print 'Working on ',name, 'in fields', fields,'current selfcal status',selfcal_status
    
  
    workdir=basedir+'/'+name
    try:
        os.mkdir(workdir)
    except OSError:
        warn('Working directory already exists')
    print 'In directory', os.getcwd()
    os.chdir(workdir)
    # Update status to running here
    selfcal_status = 'STARTED'
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    extractdict['selfcal_status'] = selfcal_status
    sdb.db_set('reprocessing',extractdict)
    sdb.close()
    print 'Updated status to STARTED for',name
    time.sleep(2.0)
    
    
    fieldstring = ''
    for field in fields:
        # WANT TO MAKE THIS INTO A RSYNC SO THAT IT CAN BE DONE OUTSIDE LEIDEN
        #os.system('cp -r %s/%s %s'%(inarchivedir,field,workdir))
        observations = glob.glob('%s/%s/%s/%s_%s*archive*'%(inarchivedir,name,field,field,name))
        for observation in observations:                     
            report('Copying data from %s'%observation)
            #'+ inarchivedir +'/' + name + '/' + field +'/' + field +'_'+name +'.dysco.sub.shift.avg.weights.ms.archive')
            do_rsync_download(observation.split('/')[-1],inarchivedir +'/'+name + '/'+field +'/',workdir)
            #do_rsync_download(field +'_'+name +'.dysco.sub.shift.avg.weights.ms.archive',inarchivedir +'/'+name + '/'+field +'/',workdir)
            fieldstring += observation.split('/')[-1] + ' '
            #'%s_%s.dysco.sub.shift.avg.weights.ms.archive '%(field,name)
    fieldstring = fieldstring[:-1]

    # Update status to copied here
    report('Updating %s status to copied'%name)
    selfcal_status = 'COPIED'
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    extractdict['selfcal_status'] = selfcal_status
    uvminstr =  str(extractdict['uvmin'])
    uvmin =  extractdict['uvmin']
    sdb.db_set('reprocessing',extractdict)
    sdb.close()
    
    
    # Create boxfile
    report('Create ds9 region file for extraction')
    create_ds9_region('%s.ds9.reg'%name,extractdict['ra'],extractdict['decl'],extractdict['size'])
    

    # Run subtract code
    print os.getcwd(), 'working here'
    
    
    if uvmin > 0.0:
       print ('runwsclean.py --uvmin=%s -b  %s.ds9.reg -i %s %s'%(uvminstr,name,name+"_image",fieldstring))
       os.system('runwsclean.py --uvmin=%s -b  %s.ds9.reg -i %s %s'%(uvminstr,name,name+"_image",fieldstring))
    else:    
       print ('runwsclean.py -b  %s.ds9.reg -i %s %s'%(name,name+"_image",fieldstring))
       os.system('runwsclean.py -b  %s.ds9.reg -i %s %s'%(name,name+"_image",fieldstring))

    report('Archiving the results to %s'%outarchivedir)
    os.chdir(workdir)
    f = glob.glob('%s.ds9.tar.gz'%(name))
    do_rsync_upload(name,outarchivedir,f)
    

    # update the database to give success
    selfcal_status = 'SDONE'
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    extractdict['selfcal_status'] = selfcal_status
    sdb.db_set('reprocessing',extractdict)
    sdb.close()
    print 'Updated status to SDONE for',name
    

if __name__=='__main__':
    target = get_next_selfcalibration()['id']
    print target
    # Takes the targetname, the current directory (the working directory), and the directory that contains the LoTSS-DR2 archive
    do_run_selfcal(target,os.getcwd(),'/disks/paradata/shimwell/LoTSS-DR2/archive_extract/','/disks/paradata/shimwell/LoTSS-DR2/archive_extract/')
