#!/usr/bin/env python

# download a field by finding all the observations from the database that go with that field

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
import sys
import os
from surveys_db import SurveysDB,tag_field
from download import download_dataset

def download_field(fname,basedir=None,force=False):

    # check database
    if basedir is None:
        print('No basedir supplied, working in current directory...')
        basedir='.'
    workdir=basedir+'/'+fname
    with SurveysDB() as sdb:
        result=sdb.get_field(fname)
        if result is None:
            print('Field',fname,'does not exist in the database')
            sys.exit(1)
        if result['status']!='Not started' and result['status']!='D/L failed':
            print('Field',fname,'has status',result['status'])
            if not force:
                return False
        # get the ids of the observations
        sdb.cur.execute('select * from observations where field=%s and (status="DI_processed" or status="Archived")',(fname,))
        obs=sdb.cur.fetchall()
        if len(obs)>0:
            result['status']='Downloading'
            if not os.path.isdir(workdir):
                os.mkdir(workdir)
            tag_field(sdb,result)
            sdb.set_field(result)
        else:
            print('No downloadable observations for this field')
            sys.exit(3)

    # now do the download for each field
    overall_success=True
    for o in obs:
        print('Downloading observation ID L'+str(o['id']))
        for prefix in ['','prefactor_v1.0/']:
            success=download_dataset('https://lofar-webdav.grid.sara.nl','/SKSP/'+prefix+'L'+str(o['id'])+'/',workdir=workdir)
            if success:
                break
            else:
                print('URL failed, trying alternative')
            
        if not success:
            print('Download failed')
        overall_success=overall_success and success

    with SurveysDB() as sdb:
        if overall_success:
            result['status']='Downloaded'
            tag_field(sdb,result,workdir=workdir)
        else:
            result['status']='D/L failed'
            tag_field(sdb,result,workdir=workdir)
        sdb.set_field(result)

    return overall_success
        
if __name__=='__main__':
    download_field(sys.argv[1],force=(len(sys.argv)>2))
    
