#!/usr/bin/env python
#Run job

from surveys_db import update_status
from auxcodes import report,warn,die
import sys
import os
import glob

def do_run_job(name,basedir,qsubfile=None,do_field=True):
    config=''
    workdir=basedir+'/'+name
    g=glob.glob(workdir+'/tier1*.cfg')
    if len(g)>0:
        print 'Local config file exists, using that'
        config=',CONFIG='+g[0]
    if qsubfile is None:
        qsubfile='/home/mjh/pipeline-master/ddf-pipeline/torque/pipeline.qsub'
    report('Submit job')
    os.system('qsub -N ddfp-'+name+' -v WD='+workdir+config+' '+qsubfile)
    if do_field:
        update_status(name,'Queued',workdir=workdir)
    
if __name__=='__main__':
    name=sys.argv[1]
    try:
        qsubfile=sys.argv[2]
    except:
        qsubfile=None

    do_run_job(name,'/beegfs/car/mjh',qsubfile)
