#!/usr/bin/env python

#insert a field and priority into the list for processing

from __future__ import print_function
import sys
from surveys_db import SurveysDB


name=sys.argv[1]
if name[0]=='L':
    name=name[1:]
id=int(name)
priority=int(sys.argv[2])

sdb=SurveysDB()
idd=sdb.get_id(id)
if idd is not None:
    print('Field already exists in the database! (Status is %s)' % idd['status'])
else:
    idd=sdb.create_id(id)
    idd['status']='Preprocessed'
    idd['priority']=priority
    sdb.set_id(idd)

sdb.close()
