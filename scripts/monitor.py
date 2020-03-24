#!/usr/bin/env python

# Code to run continuously and keep an eye on the state of the queue
# Download new files if needed
# archive using the upload script

from __future__ import print_function
from __future__ import absolute_import
from time import sleep
import datetime
from surveys_db import SurveysDB,get_next
import os
import threading
from run_pipeline import do_run_pipeline
from upload import do_upload,do_upload_compressed,do_delete_keep
import glob
import MySQLdb

queuelimit=10
cluster=os.environ['DDF_PIPELINE_CLUSTER']
download_thread=None
download_name=None
upload_thread=None
upload_name=None
# check_dict contains times to check and repeat counts for potential
# problem fields
check_dict={}
check_interval=6000
maxcount=5
basedir='/beegfs/car/mjh'

while True:

    try:
        with SurveysDB(readonly=True) as sdb:
            sdb.cur.execute('select * from fields where status!="Not started" and clustername="%s" order by priority desc, end_date' % cluster)
            result=sdb.cur.fetchall()
    except MySQLdb.OperationalError as e:
        print('Database not available! -- sleeping',e)
        sleep(240)
        continue
    
    d={}
    for r in result:
        status=r['status']
        if status in d:
            d[status]=d[status]+1
        else:
            d[status]=1
    print('\n\n-----------------------------------------------\n\n')
    print('DDF-pipeline status on cluster',cluster)
    print(datetime.datetime.now())
    print()
    for k in sorted(d.keys()):
        print("%-20s : %i" % (k,d[k]))

    if len(check_dict)>0:
        print("%-20s : %i" % ('In check list',len(check_dict)))
    
    if download_thread is not None:
        print('Download thread is running (%s)' % download_name)
    if upload_thread is not None:
        print('Upload thread is running (%s)' % upload_name)
        
    if download_thread is not None and not download_thread.isAlive():
        print('Download thread seems to have terminated')
        if download_name in check_dict:
            _,count=check_dict[download_name]
            count+=1
        else:
            count=0
        check_dict[download_name]=(datetime.datetime.now(),count)
        download_thread=None

    if upload_thread is not None and not upload_thread.isAlive():
        print('Upload thread seems to have terminated')
        upload_thread=None

    # do the check for DL failed fields
    for r in result:
        if r['id'] in check_dict:
            if r['status']=='D/L failed':
                print('Field',r['id'],'in check_dict failed download!')
                # check the count and time
                ttime,count=check_dict[r['id']]
                dt=datetime.datetime.now()-ttime
                if dt.total_seconds()>check_interval and count<maxcount:
                    # reset status so it's eligible for retry
                    print('Recheck time -- resetting its status')
                    r['status']='Not started'
                    with SurveysDB() as sdb:
                        sdb.set_field(r)

            else:
                # fields with any other status should be removed from check_dict
                print('Removing',r['id'],'from check_dict')
                del check_dict[r['id']]
    
    if ('Queued' not in d or d['Queued']<queuelimit) and download_thread is None:
        download_name=get_next()
        if download_name is not None:
            print('We need to download a new file (%s)!' % download_name)
            download_thread=threading.Thread(target=do_run_pipeline, args=(download_name,basedir))
            download_thread.start()
    

    if 'Complete' in d and upload_thread is None:
        for r in result:
            if r['status']=='Complete' and r['archive_version']==4:
                upload_name=r['id']
                kw={}
                print('We need to upload a new file (%s)!' % upload_name)
                if os.path.isdir(basedir+'/'+upload_name+'/KEEP'):
                    # remade directory, remove previous files
                    do_delete_keep(upload_name,basedir)
                    kw['skipstokes']=True
                upload_thread=threading.Thread(target=do_upload, args=(upload_name,basedir),kwargs=kw)
                upload_thread.start()
                break
    '''
    if upload_thread is None:
        for r in result:
            if r['archive_version']<2 and len(glob.glob(basedir+'/'+r['id']+'/*.archive'))>0:
                upload_name=r['id']
                print 'We need to update the archive version for %s' % upload_name
                upload_thread=threading.Thread(target=do_upload_compressed, args=(upload_name,basedir))
                upload_thread.start()
                break
    '''        
    print('\n\n-----------------------------------------------\n\n')
        
    sleep(300)
