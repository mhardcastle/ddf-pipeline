from surveys_db import SurveysDB
from astropy.coordinates import SkyCoord
import astropy.units as u

with SurveysDB() as sdb:
    sdb.cur.execute('select * from fields order by ra')
    result=sdb.cur.fetchall()

    ra=[]
    dec=[]
    for r in result:
        ra.append(r['ra'])
        dec.append(r['decl'])

    sc=SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg),frame='icrs')
    ls=sc.galactic.l.value
    bs=sc.galactic.b.value
    
    for i,r in enumerate(result):
        r['gal_l']=ls[i]
        r['gal_b']=bs[i]
        sdb.set_field(r)
        

