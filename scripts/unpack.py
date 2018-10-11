#!/usr/bin/python

import glob
import os
from surveys_db import update_status,use_database

def unpack_db_update():
    update_status(None,'Unpacked')
    
def unpack(workdir='.'):
    # unpack all the files in workdir
    destdir='prefactor/results/'
    files=glob.glob(workdir+'/*.tar.gz')
    if len(files)==0:
        # occasionally they're not compressed??
        files=glob.glob(workdir+'/*.tar')
    if len(files)==0:
        raise RuntimeError('Cannot find files to unpack')
    for f in files:
        fn=os.path.basename(f)
        print 'Unpacking',fn
        os.system('cd '+workdir+'; tar xf '+f)
        if os.path.isdir(workdir+'/prefactor'):
            os.system('cd '+workdir+'; mv prefactor/results/*.ms .')
        elif os.path.isdir(workdir+'/scratch'):
            os.system('cd '+workdir+'; mv scratch/*/*/*/Output/*.ms .')
        elif len(glob.glob('*.ms'))>0:
            # they can appear in the root directory!
            pass
        else:
            raise RuntimeError('Cannot find unpacked ms files')
            

if __name__=='__main__':
    unpack()
    if use_database():
        unpack_db_update()
