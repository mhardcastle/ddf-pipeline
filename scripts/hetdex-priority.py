#!/usr/bin/python
#Give HETDEX fields priority in the database

import sys
from surveys_db import SurveysDB

lines=[l.rstrip() for l in open(sys.argv[1]).readlines()]

sdb=SurveysDB()

for l in lines:
    id=int(l[1:])
    sdb.cur.execute('update fields,observations set fields.priority=10 where observations.field=fields.id and observations.id=%s',(id,))

sdb.close()
