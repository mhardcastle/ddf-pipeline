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
from reprocessing_utils import *


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
    executionstr = 'sub-sources-outside-region.py -b %s.ds9.reg -p %s'%(target,target)
    print(executionstr)
    result=os.system(executionstr)
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
    os.system('mkdir %s'%field)
    os.chdir(field)

    update_reprocessing_extract(target,field,'STARTING')

    prepare_field(field,startdir +'/'+target+'/'+field)

    update_reprocessing_extract(target,field,'STARTED')

    create_ds9_region('%s.ds9.reg'%(target),ra,dec,size)

    do_run_extract(field,target)

    resultfiles = glob.glob('*sub*archive*')
    resultfilestar = []
    for resultfile in resultfiles:
        d=os.system('tar -cvf %s.tar %s'%(resultfile,resultfile))
        if d!=0:
            raise RuntimeError('Tar of %s failed'%resultfile)	
        resultfilestar.append('%s.tar'%resultfile)

    do_rclone_extract_upload(field,os.getcwd(),resultfilestar,target)

    update_reprocessing_extract(target,field,'EDONE')

