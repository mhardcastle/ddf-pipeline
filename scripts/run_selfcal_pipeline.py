#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

import time
from surveys_db import *
import sys
import os
import glob
from astropy.coordinates import SkyCoord
from astropy import units as u

def check_output_ada(cname):
	rclonepath = os.environ['RCLONE_CONFIG_DIR']
	os.system('ada --tokenfile  --config=%s/maca_sksp_disk_extract.conf --longlist /%s/* > extract_files.list'%(rclonepath,cname))
	tmpfile = open('extract_files.list','r')
	tmpobs = []
	for line in tmpfile:
		line = line[:-1]
		tmpobs.append(line)
	return tmpobs

def download_extract(cname,msfilename):
	rclonepath = os.environ['RCLONE_CONFIG_DIR']
	os.system('rclone  --multi-thread-streams 1 --config=%s/maca_sksp_disk_extract.conf copy maca_sksp_disk_extract:%s/%s .'%(rclonepath,cname,msfilename))
	print('rclone  --multi-thread-streams 1 --config=%s/maca_sksp_disk_extract.conf copy maca_sksp_disk_extract:%s/%s .'%(rclonepath,cname,msfilename))
	tarfiles = glob.glob('*tar')
	for tarfile in tarfiles:
		print('tar -xf %s'%tarfile)
		os.system('tar -xf %s'%tarfile)
		os.system('rm %s'%tarfile)

def upload_extract(cname,uploadfilename):
	rclonepath = os.environ['RCLONE_CONFIG_DIR']
	print('rclone  --multi-thread-streams 1 --config=%s/maca_sksp_disk_extract.conf copy %s maca_sksp_disk_extract:%s/selfcal/'%(rclonepath,uploadfilename,cname))
	os.system('rclone  --multi-thread-streams 1 --config=%s/maca_sksp_disk_extract.conf copy %s maca_sksp_disk_extract:%s/selfcal/'%(rclonepath,uploadfilename,cname))

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


def do_run_selfcal(name,basedir):
    startdir = os.getcwd()
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    sdb.close()
    fields = extractdict['fields'].split(',')
    selfcal_status = extractdict['selfcal_status']
    extract_status = extractdict['extract_status'].split(',')
    try:
        bad_pointings = extractdict['bad_pointings'].split(',')
    except AttributeError:
        bad_pointings = ['']

    print('Populating the selfcal pointings -- a copy of fields but excluding bad_pointings')
    selfcal_pointings = ''
    for field in fields:
        if field not in bad_pointings:
            selfcal_pointings+= '%s,'%field
    selfcal_pointings = selfcal_pointings[:-1]
    
    
    print('Working on ',name, 'in fields', fields,'bad pointings',bad_pointings,'selfcal_pointings',selfcal_pointings,'current selfcal status',selfcal_status)
    
  
    workdir=basedir+'/'+name
    try:
        os.mkdir(workdir)
    except OSError:
        print('Working directory already exists')
    print('In directory', os.getcwd())
    os.chdir(workdir)
    print('In directory', os.getcwd())
    # Update status to running here
    selfcal_status = 'STARTED'
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    extractdict['selfcal_status'] = selfcal_status
    extractdict['selfcal_pointings'] = selfcal_pointings
    sdb.db_set('reprocessing',extractdict)
    sdb.close()
    print('Updated status to STARTED for',name)
    
    print('Starting rsync')
    fieldstring = ''

    for fieldid, field in enumerate(fields):
        print(field, fields)
    
        sdb=SurveysDB()
        extractdict = sdb.get_reprocessing(name)
        sdb.close()
        extract_status = extractdict['extract_status'].split(',')

        print('FIELDS', fields)
        print('SELFCAL_POINTINGS', selfcal_pointings)
        print('BAD_POINTINGS',bad_pointings)
        print('EXTRACT STATUS', extract_status)
        
        print('Copying data from %s'%field)
        download_extract(name,field)
        fieldfiles = glob.glob('*%s*'%field)
        for fieldfile in fieldfiles:
	        fieldstring += fieldfile + ' '
    fieldstring = fieldstring[:-1]

    # Update status to copied here
    print('Updating %s status to copied'%name)
    selfcal_status = 'COPIED'
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    extractdict['selfcal_status'] = selfcal_status
    uvminstr =  str(extractdict['uvmin'])
    uvmin =  extractdict['uvmin']
    sdb.db_set('reprocessing',extractdict)
    sdb.close()
    
    
    # Create boxfile
    print('Create ds9 region file for extraction')
    create_ds9_region('%s.ds9.reg'%name,extractdict['ra'],extractdict['decl'],extractdict['size'])
    

    # Run subtract code
    print(os.getcwd(), 'working here')
    selfcalscriptpath = os.environ['SELFCAL_CONFIG_DIR']

    os.system('cp %s/*.py .'%selfcalscriptpath)

    if uvmin == None:
       uvmin = 0.0
    if uvmin > 0.0:
       print ('python facetselfcal.py --auto  --uvmin=%s --remove-flagged-from-startend --helperscriptspath %s -b  %s.ds9.reg -i %s %s'%(uvminstr,selfcalscriptpath,name,name+"_image",fieldstring))
       excom = 'python facetselfcal.p y--auto --uvmin=%s --remove-flagged-from-startend  --helperscriptspath %s -b  %s.ds9.reg -i %s %s'%(uvminstr,selfcalscriptpath,name,name+"_image",fieldstring)
    else:    
       print ('python  facetselfcal.py --auto  --remove-flagged-from-startend  --helperscriptspath %s -b  %s.ds9.reg -i %s %s'%(selfcalscriptpath,name,name+"_image",fieldstring))
       excom = 'python  facetselfcal.py --auto --remove-flagged-from-startend  --helperscriptspath %s -b  %s.ds9.reg -i %s %s'%(selfcalscriptpath,name,name+"_image",fieldstring)

    os.system(excom)

    print('Archiving the results to SURF')
    os.chdir(workdir)
    f = glob.glob('%s.ds9.tar.gz'%(name)) + glob.glob('%s_image_9.png'%(name))

    for uploadfilename in f:
         upload_extract(name,uploadfilename)
    #do_rsync_upload(name,outarchivedir,f)
    
    # update the database to give success
    selfcal_status = 'SDONE'
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    extractdict['selfcal_status'] = selfcal_status
    sdb.db_set('reprocessing',extractdict)
    sdb.close()
    print('Updated status to SDONE for',name)
    

if __name__=='__main__':
    target = get_next_selfcalibration()['id']
    print(target)
    # Takes the targetname, the current directory (the working directory), and the directory that contains the LoTSS-DR2 archive
    do_run_selfcal(target,os.getcwd())
