#!/usr/bin/python
# Check what needs doing next for RGZ.

from surveys_db import SurveysDB
from astropy.table import Table
from panoptes_client import Panoptes, Project, SubjectSet, Subject, Workflow
import os
import glob
from time import sleep

target='/data/lofar/DR2/RGZ/'

def update_status(id,**kwargs):
    with SurveysDB() as sdb:
        r=sdb.get_field(id)
        for k in kwargs:
            r[k]=kwargs[k]
        sdb.set_field(r)

def create_rgz_input(id):
    t=Table.read(os.environ['LOTSS_COMPONENT_CATALOGUE'])
    filt=t['Total_flux']>8
    filt&=t['Mosaic_ID']==id
    filt&=t['Maj']>15
    filt&=t['Peak_flux']>2*t['Isl_rms']
    st=t[filt]
    print 'LGZ targets',len(st)
    tdir=target+id
    if not os.path.isdir(tdir):
        os.mkdir(tdir)
    st.write(tdir+'/'+id+'.fits',overwrite=True)
    return len(st)

def download(id):
    print 'Checking sources in %s for LGZ input' % id
    n=create_rgz_input(id)
    print 'Downloading image files'
    update_status(id,rgz_sources=n,gz_status="Downloading")
    os.system('cd %s; python %s/utils/download_image_files_legacy.py %s' % (target+id, lgz, id+'.fits'))
    update_status(id,gz_status="Downloaded")

def make_images(id):
    print 'Making images for %s' % id
    update_status(id,gz_status="Creating images")
    os.system('cd %s; python %s/lgz_create/make_overlays_legacy_cscale.py %s 0 9999' % (target+id, lgz, id+'.fits'))
    update_status(id,gz_status='Created')

def upload_images(id):
    print 'Create subject set and upload images for',id
    update_status(id,gz_status='Uploading')
    wd=os.getcwd()
    Panoptes.connect(username='mjh22',password=os.environ['PANOPTES_PASSWORD'])
    os.chdir(target+id)
    project = Project.find(slug='chrismrp/radio-galaxy-zoo-lofar')
    subject_set = SubjectSet()

    subject_set.display_name=id
    subject_set.links.project=project
    subject_set.save()
    print 'Made subject set'
    new_subjects = []
    g=glob.glob('*-manifest.txt')
    for f in g:
        bits=open(f).readlines()[0].split(',')
        metadata={'subject_id':int(bits[0]),'ra':float(bits[5]),'dec':float(bits[6]),'#size':float(bits[7]),'source_name':bits[4]}
        print 'Doing',bits[4]
        subject = Subject()
        subject.links.project = project
        subject.metadata.update(metadata)
        for location in bits[1:4]:
            subject.add_location(location)
        subject.save()
        new_subjects.append(subject)

    subject_set.add(new_subjects)

    workflow=Workflow(11973)
    workflow.links.subject_sets.add(subject_set)
    update_status(id,gz_status='In progress')
    
    
lgz=os.environ['LGZPATH']
limit=1000

while True:

    with SurveysDB(readonly=True) as sdb:
        sdb.cur.execute('select id,gz_status,weave_priority,rgz_sources from fields where dr2_final_mosaic=1 and dr2=1 and weave_priority is not NULL order by weave_priority')
        results=sdb.cur.fetchall()

    total=0
    for r in results:
        if r['gz_status']=='In progress' and r['rgz_sources'] is not None:
            total+=r['rgz_sources']
    print 'Total sources in progress',total,'limit is',limit
    for r in results:
        if r['gz_status'] in ['Complete','In progress']:
            print r['id'],r['gz_status']
            continue
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

    if r is not None and total<limit:
        id=r['id']
        if r['gz_status'] is None:
            download(id)
        elif r['gz_status']=='Downloaded':
            make_images(id)
        elif r['gz_status']=='Created':
            upload_images(id)
        else:
            print 'Cannot do anything with status',r['gz_status']

    sleep(240)
