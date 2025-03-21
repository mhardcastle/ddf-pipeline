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
import argparse

def download_field(fname,basedir=None,force=False,use_http=False,macaroons=['maca_sksp_tape_spiderlinc.conf','maca_sksp_tape_spiderpref3.conf','maca_sksp_distrib_Pref3.conf'],update_status=True):

    # check database
    if basedir is None:
        print('No basedir supplied, working in current directory...')
        basedir='.'
    workdir=basedir+'/'+fname
    with SurveysDB(readonly=not update_status) as sdb:
        result=sdb.get_field(fname)
        if result is None:
            print('Field',fname,'does not exist in the database')
            sys.exit(1)
        if result['status']!='Not started' and result['status']!='D/L failed':
            print('Field',fname,'has status',result['status'])
            if not force:
                return False
        # get the ids of the observations
        if force:
            sdb.cur.execute('select * from observations where field=%s',(fname,))
        else:
            sdb.cur.execute('select * from observations where field=%s and (status="DI_processed" or status="Archived")',(fname,))
        obs=sdb.cur.fetchall()
        if len(obs)>0:
            if not os.path.isdir(workdir):
                os.mkdir(workdir)
            if update_status:
                result['status']='Downloading'
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
        success=False # will be set to true if rclone works and we can
                      # find the dataset there
                      
        for macaroon in macaroons:
            try:
                rc=RClone(macaroon,debug=True)
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
                    print('rclone found data but download failed')
                else:
                    print('rclone download succeeded')
                    success=True
                    break # out of rclone loop
            
        if not success and use_http:
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

    if update_status:
        with SurveysDB() as sdb:
            if overall_success:
                print('Download overall successful!')
                result['status']='Downloaded'
                tag_field(sdb,result,workdir=workdir)
            else:
                print('Download overall failed! (some fields did not download or were not available.')
                result['status']='D/L failed'
                tag_field(sdb,result,workdir=workdir)
            sdb.set_field(result)

    return overall_success

def download_obsid(fname,basedir=None,force=False):

    # check database
    if basedir is None:
        print('No basedir supplied, working in current directory...')
        basedir='.'
    workdir=basedir+'/'+fname
    with SurveysDB() as sdb:
        result=sdb.get_observation(fname)
        if result is None:
            print('Observation',fname,'does not exist in the database')
            sys.exit(1)
        if result['status']!='Archived':
            print('Observation',fname,'has status',result['status'])
            if not force:
                return False
        # get the ids of the observations
        sdb.cur.execute('select * from observations where id=%s and (status="DI_processed" or status="Archived")',(fname,))
        obs=sdb.cur.fetchall()
        if len(obs)==0:
            print('No downloadable observations for this field')
            sys.exit(3)

    # now do the download for each field
    overall_success=True
    rclone_works=True
    print(obs)
    for o in obs:
        obsname='L'+str(o['id'])
        print('Downloading observation ID',obsname)
        success=False
        # first try rclone

        for macaroon in ['maca_sksp_tape_spiderlinc.conf','maca_sksp_tape_spiderpref.conf']:
            try:
                rc=RClone(macaroon,debug=True)
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
                    print('rclone found data but download failed')
                else:
                    print('rclone download succeeded')
                    success=True
                    break # out of rclone loop
            
        if not success:
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

    return overall_success
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Download LoTSS field')
    parser.add_argument('--force',action='store_true',help='Force download irrespective of field status')
    parser.add_argument('--no_update',action='store_true',help='No database update')
    parser.add_argument('--use_http',action='store_true',help='Use HTTP download method (obsolete)')
    parser.add_argument('--macaroons', metavar='filename.conf', nargs='+', help='Macaroons to use')
    parser.add_argument('field', help='field name')
    args = parser.parse_args()
    field=args.field
    update_status=not args.no_update
    del(args.field)
    del(args.no_update)
    if not args.macaroons:
        del(args.macaroons)
    download_field(field,update_status=update_status,**vars(args))
    #download_field(args.field,force=args.force,use_http=args.use_http,update_status=not args.no_update)
    
