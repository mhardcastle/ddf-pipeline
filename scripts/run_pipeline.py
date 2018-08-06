#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from auxcodes import report,warn,die
from surveys_db import use_database,update_status
from download import download_dataset
from download_field import download_field
from unpack import unpack
from make_mslists import make_list,list_db_update
import sys
import os

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
    unpack(workdir=workdir)
    if do_field:
        update_status(name,'Unpacked',workdir=workdir)

    report('Deleting tar files')
    os.system('rm '+workdir+'/*.tar.gz')

    report('Making ms lists')
    success=make_list(workdir=workdir)
    if do_field:
        list_db_update(success,workdir=workdir)

    if success:
        report('Submit job')
        os.system('qsub -N ddfp-'+name+' -v WD='+workdir+' '+qsubfile)
        if do_field:
            update_status(name,'Queued',workdir=workdir)

    else:
        die('make_list could not construct the MS list',database=False)

if __name__=='__main__':
    do_run_pipeline(sys.argv[1],'/beegfs/car/mjh')
