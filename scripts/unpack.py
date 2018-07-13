#!/usr/bin/python

import glob
import os
from surveys_db import update_status,use_database

def unpack_db_update():
    update_status(None,'Unpacked')
    
def unpack():
    # unpack all the files in current working directory
    destdir='prefactor/results/'
    files=glob.glob('*.tar.gz')
    for f in files:
        print 'Unpacking',f
        os.system('tar xf '+f)
        os.system('mv '+destdir+'* .')

if __name__=='__main__':
    unpack()
    if use_database():
        unpack_db_update()
