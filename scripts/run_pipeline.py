#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from auxcodes import report,warn,die
from surveys_db import SurveysDB,update_status
from download import download_dataset
from download_field import download_field
from run_job import do_run_job
from unpack import unpack
from make_mslists import make_list,list_db_update
from average import average
from auxcodes import MSList
import numpy as np
import sys
import os
import glob

def make_custom_config(name,workdir,do_field,averaged=False):
    if do_field:
        with SurveysDB() as sdb:
            idd=sdb.get_field(name)
            no_wenss=((idd['decl']<32) | (idd['decl']>72))
    else:
        no_wenss=False

    if no_wenss:
        template=os.environ['DDF_DIR']+'/ddf-pipeline/examples/tier1-jul2018-NVSS.cfg'
    else:
        template=os.environ['DDF_DIR']+'/ddf-pipeline/examples/tier1-jul2018.cfg'

    lines=open(template).readlines()
    outfile=open(workdir+'/tier1-config.cfg','w')
    for l in lines:
        if 'colname' in l and averaged:
            outfile.write('colname=DATA\n')
        else:
            outfile.write(l)

def do_run_pipeline(name,basedir):

    if name[0]!='P' and name[0]!='L':
        die('This code should be used only with field or observation names',database=False)

    do_field=(name[0]=='P')

    try:
        qsubfile=sys.argv[2]
    except:
        qsubfile='/home/mjh/pipeline-master/ddf-pipeline/torque/pipeline.qsub'


    workdir=basedir+'/'+name
    try:
        os.mkdir(workdir)
    except OSError:
        warn('Working directory already exists')

    report('Downloading data')
    if do_field:
        success=download_field(name,basedir=basedir)
    else:
        success=download_dataset('https://lofar-webdav.grid.sara.nl','/SKSP/'+name+'/',basedir=basedir)

    if not success:
        die('Download failed, see earlier errors',database=False)

    report('Unpacking data')
    try:
        unpack(workdir=workdir)
    except RuntimeError:
        if do_field:
            update_status(name,'Unpack failed',workdir=workdir)
        raise
    if do_field:
        update_status(name,'Unpacked',workdir=workdir)

    report('Deleting tar files')
    os.system('rm '+workdir+'/*.tar.gz')
    os.system('rm '+workdir+'/*.tar')

    averaged=False
    report('Checking structure')
    g=glob.glob(workdir+'/*.ms')
    msl=MSList(None,mss=g)
    dysco=np.any(msl.dysco)
    uobsids=set(msl.obsids)
    for thisobs in uobsids:
        # check one MS with each ID
        for m,ch,o,hc in zip(msl.mss,msl.channels,msl.obsids,msl.hascorrected):
            if o==thisobs:
                if not(hc):
                    print 'MS',m,'has no corrected_data column, force use of DATA'
                    averaged=True
                channels=len(ch)
                print 'MS',m,'has',channels,'channels'
                if channels>20:
                    update_status(name,'Averaging',workdir=workdir)
                    print 'Averaging needed for',thisobs,'!'
                    averaged=True
                    average(wildcard=workdir+'/*'+thisobs+'*')
                    os.system('rm -r '+workdir+'/*'+thisobs+'*pre-cal.ms')
                break
    
    report('Making ms lists')
    success=make_list(workdir=workdir)
    if do_field:
        list_db_update(success,workdir=workdir)
    if not success:
        die('make_list could not construct the MS list',database=False)
        
    report('Creating custom config file from template')
    make_custom_config(name,workdir,do_field,averaged)
    
    # now run the job
    do_run_job(name,basedir=basedir,qsubfile=None,do_field=do_field,dysco=dysco)


if __name__=='__main__':
    do_run_pipeline(sys.argv[1],'/beegfs/car/mjh')
