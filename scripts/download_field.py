#!/usr/bin/env python

# download a field by finding all the observations from the database that go with that field

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
import sys
import os
from surveys_db import SurveysDB,tag_field
from download import download_dataset
from rclone import RClone

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
    rclone_works=True
    for o in obs:
        obsname='L'+str(o['id'])
        print('Downloading observation ID',obsname)

        # first try rclone

        try:
            rc=RClone('maca_sksp_tape_spiderpref.conf',debug=True)
        except RuntimeError as e:
            print('rclone setup failed, probably RCLONE_CONFIG_DIR not set:',e)
            rclone_works=False

        if rclone_works:
            try:
                remote_obs=rc.get_dirs()
            except OSError as e:
                print('rclone command failed, probably rclone not installed or RCLONE_COMMAND not set:',e)
                rclone_works=False
        
        if rclone_works and obsname in remote_obs:
            print('Data available in rclone repository, downloading!')
            d=rc.execute_live(['-P','copy',rc.remote+'/'+obsname,workdir])
            if d['err'] or d['code']!=0:
                print('rclone download failed')
                success=False
            else:
                print('rclone download succeeded')
                success=True

        else:
            # revert to download method
            for prefix in ['','prefactor_v1.0/','prefactor_v3.0/']:
                success=download_dataset('https://lofar-webdav.grid.sara.nl:2880','/SKSP/'+prefix+obsname+'/',workdir=workdir)
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
    
