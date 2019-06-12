#!/usr/bin/python

# Run the quality pipeline
# Herts-only code

from surveys_db import SurveysDB
import os
from subprocess import call,check_output
from time import sleep

queued=[]
queue=check_output('qstat -a',shell=True).split('\n')
for l in queue:
    if 'qual-' in l:
        bits=l.split()
        field=bits[3][5:]
        print 'Found',field,'already in queue'
        queued.append(field)

with SurveysDB() as sdb:
    sdb.cur.execute('select * from fields left join quality on quality.id=fields.id where status="Archived" and archive_version=4 and quality.dr is NULL')
    results=sdb.cur.fetchall()

for r in results:
    id=r['id']
    dir='/data/lofar/DR2/fields/'+id
    if os.path.isfile(dir+'/image_full_ampphase_di_m.NS.cat.reg'):
        print id,'has the quality catalogue'
    else:
        print id,'does not have the quality catalogue'
    if id in queued:
        print 'Not queueing it as it is already queued'
    else:
        os.system('qsub -l nodes=1:ppn=1 -l pmem=8gb -N qual-%s -v WD=%s ~/pipeline-master/ddf-pipeline/torque/quality.qsub' % (id,dir))
