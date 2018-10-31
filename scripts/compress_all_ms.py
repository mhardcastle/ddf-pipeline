#!/usr/bin/python

# Torque-based script to generate compressed versions of all archived
# datasets that don't already have .archive files

import os
import glob
from surveys_db import SurveysDB

with SurveysDB() as sdb:
    sdb.cur.execute('select * from fields where status="Archived" and clustername="Herts"')
    results=sdb.cur.fetchall()

for r in results:
    id=r['id']
    location=r['location']
    if glob.glob(location+'/*.archive'):
        print id,'compressed files already exist'
    else:
        print id,'should be compressed at',location
        command='qsub -N cprs-%s -v WD=%s /home/mjh/pipeline-master/ddf-pipeline/torque/compress.qsub' % (id,location)
        print command
        os.system(command)
        
