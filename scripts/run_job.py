#!/usr/bin/env python
#Run job

from __future__ import print_function
from surveys_db import update_status, SurveysDB
from auxcodes import report,warn,die
import sys
import os
import glob

def do_run_job(name,basedir,qsubfile=None,do_field=True,prefix='ddfp',dysco=False):
    config=''
    workdir=basedir+'/'+name
    g=glob.glob(workdir+'/tier1*.cfg')
    if len(g)>0:
        print('Local config file exists, using that')
        config=',CONFIG='+g[0]
    if qsubfile is None:
        qsubfile='/home/mjh/pipeline-master/ddf-pipeline/torque/pipeline.qsub'
    if dysco:
        qsubfile=qsubfile.replace('.qsub','-fdr14.qsub')
    report('Submit job')
    os.system('qsub -N '+prefix+'-'+name+' -v WD='+workdir+config+' '+qsubfile)
    if do_field:
        update_status(name,'Queued',workdir=workdir)

def rerun_select():
    with SurveysDB() as sdb:
        sdb.cur.execute('select * from fields where clustername="Herts" and status="Archived" and priority<10 and archive_version<3 order by priority desc')
        results=sdb.cur.fetchall()
    for r in results:
        name=r['id']
        print('Submitting job for',name)
        do_run_job(name,'/beegfs/car/mjh',qsubfile='/home/mjh/pipeline-master/ddf-pipeline/torque/rerun.qsub',prefix='ddfpr')
        
if __name__=='__main__':
    name=sys.argv[1]
    try:
        qsubfile=sys.argv[2]
    except:
        qsubfile=None

    do_run_job(name,'/beegfs/car/mjh',qsubfile)
