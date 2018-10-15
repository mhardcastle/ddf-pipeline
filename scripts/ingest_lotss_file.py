#!/usr/bin/python
# Ingest the LOTSS file provided by Tim

import sys
from surveys_db import SurveysDB
import requests

root='https://lofar-webdav.grid.sara.nl/SKSP/L'

lines=[l.rstrip() for l in open(sys.argv[1]).readlines()]

sdb=SurveysDB()
sdb.cur.execute('delete from fields_new')
sdb.cur.execute('delete from observations')
for l in lines[1:]:
    bits=l.split(',')
    if bits[1]=='WRONG':
        continue
    field=bits[0]
    ra=float(bits[1])
    dec=float(bits[2])
    time=float(bits[3])
    if bits[4]=='None':
        ids=[]
    else:
        ids=bits[4].split(' ')
        if ids[0]=='':
            ids=ids[1:]
    print time,field,ra,dec,ids
    sdb.cur.execute('insert into fields_new values ( %s,"Not started",%s,%s,NULL,NULL,NULL,NULL,0)', (field,ra,dec))
    for i in ids:
        page=requests.get(root+i)
        if page.status_code==200:
            status='Preprocessed'
        else:
            status='Observed'
        sdb.cur.execute('insert into observations values ( %s, %s, %s, %s, NULL, NULL )',(i,field,status,time/len(ids)))

sdb.close()
