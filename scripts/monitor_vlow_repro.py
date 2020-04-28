#!/usr/bin/python

# clone of monitor.py to do vlow monitoring

from __future__ import print_function
from time import sleep
import datetime
from surveys_db import SurveysDB,get_next
import os
import threading
from run_smooth_again import do_download,update_status
import glob

cluster=os.environ['DDF_PIPELINE_CLUSTER']
download_thread=None
download_name=None
tidy_thread=None
upload_name=None
basedir='/beegfs/car/mjh/reprocess'
totallimit=10

os.chdir(basedir)

while True:

    with SurveysDB(readonly=True) as sdb:
        sdb.cur.execute('select * from fields where vlow_reprocess!="Not needed" order by end_date desc')
        result=sdb.cur.fetchall()
        sdb.cur.execute('select * from fields where vlow_reprocess="Required" order by end_date desc')
        result2=sdb.cur.fetchall()
        if len(result2)>0:
            nextfield=result2[0]['id']
        else:
            nextfield=None

    d={}
    for r in result:
        status=r['vlow_reprocess']
        if status in d:
            d[status]=d[status]+1
        else:
            d[status]=1
    print('\n\n-----------------------------------------------\n\n')
    print('DDF-pipeline reprocess status on cluster',cluster)
    print(datetime.datetime.now())
    print()

    for k in sorted(d.keys()):
        print("%-20s : %i" % (k,d[k]))

    print()
    ksum=len(glob.glob(basedir+'/*'))
    print(ksum,'live files out of',totallimit)
    print('Next field to work on is',nextfield)
            
    if download_thread is not None:
        print('Download thread is running (%s)' % download_name)
        
    if download_thread is not None and not download_thread.isAlive():
        print('Download thread seems to have terminated')
        # run the qsub here!
        command="qsub -v FIELD=%s -N vlow-%s /home/mjh/pipeline-master/ddf-pipeline/torque/vlow_reprocess.qsub" % (download_name, download_name)
        if os.system(command)==0:
            update_status(download_name,"Queued")
        download_thread=None

    if ksum<totallimit and nextfield is not None and download_thread is None:
        download_name=nextfield
        print('We need to download a new file (%s)!' % download_name)
        download_thread=threading.Thread(target=do_download, args=(download_name,basedir))
        download_thread.start()
    
    if 'Uploaded' in d:
        for r in result:
            if r['vlow_reprocess']=='Uploaded':
                field=r['id']
                print('Tidying uploaded directory for',field)
                command='cp '+basedir+'/'+field+'/*vlow*fz /data/lofar/DR2/fields/'+field
                print('running',command)
                os.system(command)
                command='rm -r '+basedir+'/'+field
                os.system(command)
                update_status(field,'Archived')

    print('\n\n-----------------------------------------------\n\n')
        
    sleep(300)
