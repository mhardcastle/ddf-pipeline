#!/usr/bin/python

# download a field by finding all the observations from the database that go with that field

import sys
import os
from surveys_db import SurveysDB,tag_field
from download import download_dataset

def download_field(fname):

    # check database

    with SurveysDB() as sdb:
        result=sdb.get_field(fname)
        if result is None:
            print 'Field',fname,'does not exist in the database'
            sys.exit(1)
        if result['status']!='Not started':
            print 'Field',fname,'has status',result['status']
            sys.exit(2)
        # get the ids of the observations
        sdb.cur.execute('select * from observations where field=%s',(fname,))
        obs=sdb.cur.fetchall()
        if len(obs)>0:
            result['status']='Downloading'
            os.mkdir(fname)
            os.chdir(fname)
            tag_field(sdb,result)
            sdb.set_field(result)
        else:
            print 'No downloadable observations for this field'
            sys.exit(3)

    # now do the download for each field
    overall_success=True
    for o in obs:
        print 'Downloading observation ID L'+str(o['id'])
        success=download_dataset('https://lofar-webdav.grid.sara.nl','/SKSP/L'+str(o['id'])+'/')
        if success==False:
            print 'Download failed'
        overall_success=overall_success and success

    with SurveysDB() as sdb:
        if overall_success:
            result['status']='Downloaded'
            tag_field(sdb,result)
        else:
            result['status']='D/L failed'
        sdb.set_field(result)

if __name__=='__main__':

    download_field(sys.argv[1])
    
