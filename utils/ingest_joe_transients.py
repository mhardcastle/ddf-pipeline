from __future__ import absolute_import
from .surveys_db import SurveysDB
from astropy.table import Table

t=Table.read('/data/lofar/mjh/master_dynspec_target_list.fits')


with SurveysDB() as sdb:
    sdb.cur.execute('delete from transients')
    for r in t:
        name=r['name'].rstrip()
        count=0
        while True:
            tn=sdb.get_transient(name)
            if tn is None:
                break
            name=name+'_%i' % (count+2)
            count=count+1
        sdb.cur.execute('insert into transients values ("%s",%s,%s,"%s")' % (name,r['ra'],r['dec'],r['cat'].rstrip()))
