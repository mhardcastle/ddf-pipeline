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
    for f in files:
        fn=os.path.basename(f)
        print 'Unpacking',fn
        os.system('cd '+workdir+'; tar xf '+f)
        os.system('cd '+workdir+'; mv '+destdir+'* .')

if __name__=='__main__':
    unpack()
    if use_database():
        unpack_db_update()
