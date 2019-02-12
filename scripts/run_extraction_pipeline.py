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


def do_run_subtract(name,basedir,inarchivedir,outarchivedir):
    startdir = os.getcwd()
    sdb=SurveysDB()
    extractdict = sdb.get_reprocessing(name)
    sdb.close()
    fields = extractdict['fields'].split(',')
    extract_status = extractdict['extract_status'].split(',')

    print 'Working on ',name, 'in fields', fields,'which have status',extract_status
    
    for i in range(0,len(fields)):
        os.chdir(startdir)
        if extract_status[i] != 'EREADY':
            continue
        field = fields[i]

        workdir=basedir+'/'+name
        try:
            os.mkdir(workdir)
        except OSError:
            warn('Working directory already exists')
        print 'In directory', os.getcwd()
        os.chdir(workdir)
        # Update status to running here
        extract_status[i] = 'STARTED'
        sdb=SurveysDB()
        extractdict = sdb.get_reprocessing(name)
        extractdict['extract_status'] = ','.join(extract_status)
        sdb.db_set('reprocessing',extractdict)
        sdb.close()
        print 'Updated status to STARTED for',field,name
        time.sleep(2.0)
        report('Copying data from %s'%inarchivedir)
        
        # WANT TO MAKE THIS INTO A RSYNC SO THAT IT CAN BE DONE OUTSIDE LEIDEN
        os.system('cp -r %s/%s %s'%(inarchivedir,field,workdir))

        # Update status to copied here
        extract_status[i] = 'COPIED'
        sdb=SurveysDB()
        extractdict = sdb.get_reprocessing(name)
        extractdict['extract_status'] = ','.join(extract_status)
        sdb.db_set('reprocessing',extractdict)
        sdb.close()
        print 'Updated status to COPIED for',field,name


        # Create boxfile
        create_ds9_region('%s.ds9.reg'%name,extractdict['ra'],extractdict['decl'],extractdict['size'])


        # Run subtract code
        print os.getcwd(), 'working here'
        os.chdir(field)
        print ('sub-sources-outside-region.py -b %s/%s.ds9.reg -p %s'%(workdir,name,name))
        os.system('sub-sources-outside-region.py -b %s/%s.ds9.reg -p %s'%(workdir,name,name))

        # Archive the results need an rsync code this is just the *archive file that needs to be archived.
        os.system('mkdir %s/%s'%(outarchivedir,name))
        os.system('mkdir %s/%s/%s'%(outarchivedir,name,field))
        print  ('cp -r %s_%s.dysco.sub.shift.avg.weights.ms.archive %s/%s/%s'%(field,name,outarchivedir,name,field))
        os.system('cp -r %s_%s.dysco.sub.shift.avg.weights.ms.archive %s/%s/%s'%(field,name,outarchivedir,name,field))


        # update the database to give success
        extract_status[i] = 'EDONE'
        sdb=SurveysDB()
        extractdict = sdb.get_reprocessing(name)
        extractdict['extract_status'] = ','.join(extract_status)
        sdb.db_set('reprocessing',extractdict)
        sdb.close()
        print 'Updated status to EDONE for',field,name

if __name__=='__main__':
    target = get_next_extraction()['id']
    # Takes the targetname, the current directory (the working directory), and the directory that contains the LoTSS-DR2 archive
    do_run_subtract(target,os.getcwd(),'/disks/paradata/shimwell/LoTSS-DR2/archive/','/disks/paradata/shimwell/LoTSS-DR2/archive_extract/')
