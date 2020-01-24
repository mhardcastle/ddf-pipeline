# Check what needs doing next for RGZ.

from surveys_db import SurveysDB
from astropy.table import Table
import os

def create_rgz_input(id):
    t=Table.read('/data/lofar/DR2/temp/LoTSS_DR2_v0.9.srl.fits')
    filt=t['Total_flux']>8
    filt&=t['Mosaic_ID']==id
    filt&=t['Maj']>15
    filt&=t['Peak_flux']>2*t['Isl_rms']
    st=t[filt]
    print 'LGZ targets',len(st)
    tdir='/data/lofar/DR2/RGZ/'+id
    if not os.path.isdir(tdir):
        os.mkdir(tdir)
    st.write(tdir+'/'+id+'.fits',overwrite=True)
    return len(st)
    
with SurveysDB(readonly=True) as sdb:
    sdb.cur.execute('select id,gz_status,weave_priority from fields where dr2_final_mosaic=1 and dr2=1 and weave_priority is not NULL order by weave_priority')
    results=sdb.cur.fetchall()

for r in results:
    if r['gz_status'] is not None:
        print r['id'],r['gz_status']
    if r['gz_status'] is None:
        print 'First non-running field is',r['id']
        break
else:
    r=None

# Here we should:
# -- keep track of how many sources have gone into active fields
# -- if more are needed, make them, one mosaic at a time. update when the images are ready
#    -- steps are: download_image_files_legacy.py (prob about 10 mins)
#    -- make_images (make_overlays_legacy_scale.py)
#    NB IMAGEDIR and LOTSS_COMPONENT_CATALOGUE need to be set
# -- if images are ready to upload, upload them
# -- later we can also extract results from RGZ
# Do we want to update the component catalogue?

if r is not None:
    print 'Checking sources for LGZ input'
    n=create_rgz_input(r['id'])

