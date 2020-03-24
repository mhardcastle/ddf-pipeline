from __future__ import print_function
from __future__ import division
# Check for a nearby LoTSS pointing (completed, archived)

from past.utils import old_div
from surveys_db import SurveysDB
from auxcodes import sepn
import numpy as np

deg2rad=old_div(np.pi,180)
rad2deg=180.0/np.pi


ra=[]
dec=[]
id=[]
with SurveysDB(readonly=True) as sdb:
    sdb.cur.execute('select * from fields where status="Archived" or status="Completed"')
    results=sdb.cur.fetchall()

for r in results:
    ra.append(r['ra'])
    dec.append(r['decl'])
    id.append(r['id'])

'''
for l in open('/home/mjh/pipeline-master/ddf-pipeline/misc/DR2-pointings.txt').readlines():
    bits=l.rstrip().split()
    id.append(bits[0])
    ra.append(float(bits[1]))
    dec.append(float(bits[2]))
'''

ra=np.array(ra)*deg2rad
dec=np.array(dec)*deg2rad

r_ra=270*deg2rad
r_dec=45*deg2rad

sep=sepn(ra,dec,r_ra,r_dec)
i=np.argmin(sep)
print(i,id[i],sep[i]*rad2deg)

# mosaic distance is 2.2 deg
limit=2.2*deg2rad

area=0
totarea=0
ras=np.linspace(0,360,700)*deg2rad
decs=np.linspace(0,90,700)*deg2rad
for d in decs:
    a=np.cos(d)
    print('.', end=' ')
    for r in ras:
        sep=sepn(ra,dec,r,d)
        i=np.argmin(sep)
        totarea+=a
        if sep[i]<limit:
            area+=a
print()
print(area,totarea,100.0*area/totarea)
