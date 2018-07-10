#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from auxcodes import report,warn,die
from surveys_db import use_database,update_status
from download import download_dataset,download_db_create,download_db_update
from unpack import unpack,unpack_db_update
from make_mslists import make_list,list_db_update
import sys
import os

rootdir='/beegfs/car/mjh'
os.chdir(rootdir)

name=sys.argv[1]
try:
    qsubfile=sys.argv[2]
except:
    qsubfile='/home/mjh/git/ddf-pipeline/pipeline.qsub'

try:
    os.mkdir(name)
except OSError:
    warn('Working directory already exists')
    pass
os.chdir(name)
report('Downloading data')
if use_database:
    download_db_create(name)
success=download_dataset('https://lofar-webdav.grid.sara.nl','/SKSP/'+name+'/')
if use_database:
    download_db_update(name,success)
if not success:
    die('Download failed to get the right number of files')

    
report('Unpacking data')
unpack()
if use_database():
    unpack_db_update()
    
report('Deleting tar files')
os.system('rm *.tar.gz')

report('Making ms lists')
success=make_list()
if use_database():
        list_db_update(success)

if success:
    report('Submit job')
    os.system('qsub -N ddfp-'+name+' -v WD='+rootdir+'/'+name+' '+qsubfile)
    if use_database():
        update_status(name,'Queued')

else:
    die('make_list could not construct the MS list')
