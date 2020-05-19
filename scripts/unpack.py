#!/usr/bin/env python

from __future__ import print_function
import glob
import os
from surveys_db import update_status,use_database

def unpack_db_update():
    update_status(None,'Unpacked')
    
def unpack(workdir='.'):
    # unpack all the files in workdir
    destdir='prefactor/results/'
    files=glob.glob(workdir+'/*.tar.gz')+glob.glob(workdir+'/*.tar')
    if len(files)==0:
        raise RuntimeError('Cannot find files to unpack')
    for f in files:
        if 'tokens' in f:
            print('Skipping',f)
            continue
        fn=os.path.basename(f)
        print('Unpacking',fn)
        result=os.system('cd '+workdir+'; tar xf '+f)
        if result!=0:
            raise RuntimeError('Untar failed')
        if os.path.isdir(workdir+'/prefactor'):
            os.system('cd '+workdir+'; mv prefactor/results/*.ms .')
        elif os.path.isdir(workdir+'/scratch'):
            os.system('cd '+workdir+'; mv scratch/*/*/*/Output/*.ms .')
            os.system('cd '+workdir+'; mv scratch/*/*/Output/*.ms .')
            os.system('cd '+workdir+'; mv scratch/*/Output/*.ms .')
        elif os.path.isdir(workdir+'/results'):
            # prefactor3 can create things in 'results'
            os.system('cd '+workdir+'; mv results/*.ms .')            
        elif len(glob.glob(workdir+'/*.ms'))>0:
            # they can appear in the root directory!
            pass
        else:
            raise RuntimeError('Cannot find unpacked ms files')
            

if __name__=='__main__':
    unpack()
    if use_database():
        unpack_db_update()
