#!/usr/bin/python

import glob
import os
from surveys_db import SurveysDB,use_database,get_id,tag_idd

def unpack_db_update():
    id=get_id()
    sdb=SurveysDB()
    idd=sdb.get_id(id)
    idd['status']='Unpacked'
    tag_idd(sdb,idd)
    sdb.set_id(idd)
    
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
