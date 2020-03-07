from __future__ import print_function
# Compress the QU cubes for all archived data

import os
import glob
from surveys_db import SurveysDB

def compress_fits(filename,q=4):
    print('Compressing and deleting',filename)
    if not os.path.isfile(filename+'.fz'):
        command='fpack -q %i %s' % (q,filename)
        os.system(command)
    else:
        print('Compressed file already exists, skipping')
    os.remove(filename)
    
from multiprocessing import Pool

with SurveysDB() as sdb:
    sdb.cur.execute('select * from fields where status="Archived" and clustername="Herts"')
    results=sdb.cur.fetchall()

files=[]
for r in results:
    id=r['id']
    for f in ['image_full_vlow_QU.cube.dirty.fits','image_full_vlow_QU.cube.dirty.corr.fits','image_full_low_QU.cube.dirty.fits','image_full_low_QU.cube.dirty.corr.fits']:
        fname=id+'/'+f
        if os.path.isfile(fname):
            files.append(fname)

p=Pool(16)
p.map(compress_fits,files)
p.terminate()


