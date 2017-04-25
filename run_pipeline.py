#!/usr/bin/env python
# Run pipeline download/unpack steps followed by the main job

from auxcodes import report,warn
from download import download_dataset
from unpack import unpack
from make_mslists import make_list
import sys
import os

rootdir='/data/lofar/mjh'
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
download_dataset('https://lofar-webdav.grid.sara.nl','/SKSP/'+name+'/')

report('Unpacking data')
unpack()

report('Deleting tar files')
os.system('rm *.tar.gz')

report('Making ms lists')
make_list()

report('Submit job')
os.system('qsub -N ddfp-'+name+' -v WD='+rootdir+'/'+name+' '+qsubfile)
