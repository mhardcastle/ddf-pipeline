#!/usr/bin/python
# Check what needs doing next for RGZ.

from surveys_db import SurveysDB
from astropy.table import Table
from panoptes_client import Panoptes, Project, SubjectSet, Subject, Workflow
import os
import glob
from time import sleep
import datetime
import threading

target='/data/lofar/DR2/RGZ/'
lgz=os.environ['LGZPATH']

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
    result=os.system('cd %s; python %s/utils/download_image_files_legacy.py %s' % (target+id, lgz, id+'.fits'))
    if result==0:
        update_status(id,gz_status="Downloaded")
    else:
        raise RuntimeError('Download failed!')
    
    

def make_images(id):
    print 'Making images for %s' % id
    update_status(id,gz_status="Creating images")
    result=os.system('cd %s; python %s/lgz_create/make_overlays_legacy_cscale.py %s 0 9999' % (target+id, lgz, id+'.fits'))
    if result==0:
        update_status(id,gz_status='Created')
    else:
        raise RuntimeError('Make images failed!')

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
    for i,f in enumerate(g):
        bits=open(f).readlines()[0].split(',')
        metadata={'subject_id':int(bits[0]),'ra':float(bits[5]),'dec':float(bits[6]),'#size':float(bits[7]),'source_name':bits[4]}
        print 'Upload doing',bits[4],'%i/%i' % (i,len(g))
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
    print 'Done!'
    
if __name__=='__main__':    
    download_thread=None
    download_name=None
    create_thread=None
    create_name=None
    upload_thread=None
    upload_name=None
    while True:

        with SurveysDB(readonly=True) as sdb:
            sdb.cur.execute('select id,gz_status,weave_priority,rgz_sources,rgz_complete from fields where dr2=1 order by weave_priority')
            results=sdb.cur.fetchall()

        d={}
        for r in results:
            status=r['gz_status']
            if status in d:
                d[status]=d[status]+1
            else:
                d[status]=1

        limit=int(open('/home/mjh/pipeline-master/ddf-pipeline/misc/lgz-limit.txt').readlines()[0].rstrip())
                
        print '\n\n-----------------------------------------------\n\n'
        print 'LGZ status'
        print(datetime.datetime.now())
        print
        for k in sorted(d.keys()):
            print "%-20s : %i" % (k,d[k])

        total=0
        tremain=0
        ctotal=0
        ftotal=0
        for r in results:
            if r['gz_status']=='In progress' and r['rgz_sources'] is not None:
                total+=r['rgz_sources']
                tremain+=r['rgz_sources']
                if r['rgz_complete'] is not None:
                    tremain-=r['rgz_complete']
            if r['gz_status']=='Created':
                ctotal+=r['rgz_sources']
            if r['rgz_complete'] is not None:
                ftotal+=r['rgz_complete']
                
                
        print 'Total sources in fields in progress',total,'of which',tremain,'are not retired'
        print 'Non-retired lower limit is',limit
        print 'Total sources created but not uploaded',ctotal
        print 'Total sources retired',ftotal
        
        if download_thread is not None:
            print 'Download thread is running (%s)' % download_name
        if create_thread is not None:
            print 'Create thread is running (%s)' % create_name
        if upload_thread is not None:
            print 'Upload thread is running (%s)' % upload_name

        if download_thread is not None and not download_thread.isAlive():
            print 'Download thread seems to have terminated'
            download_thread=None

        if create_thread is not None and not create_thread.isAlive():
            print 'Create thread seems to have terminated'
            create_thread=None

        if upload_thread is not None and not upload_thread.isAlive():
            print 'Upload thread seems to have terminated'
            upload_thread=None

            
        non_running=None
        for r in results:
            if r['gz_status'] is None:
                non_running=r['id']
                print 'First non-running field is',r['id']
                break

        if non_running is not None and ctotal<limit and download_thread is None:
            download_name=non_running
            print 'We need to download a new file (%s)!' % download_name
            download_thread=threading.Thread(target=download, args=(download_name,))
            download_thread.start()

        if 'Downloaded' in d and create_thread is None:
            for r in results:
                if r['gz_status']=='Downloaded':
                    create_name=r['id']
                    create_thread=threading.Thread(target=make_images, args=(create_name,))
                    create_thread.start()
                    break

        if 'Created' in d and tremain<limit and upload_thread is None:
            for r in results:
                if r['gz_status']=='Created':
                    upload_name=r['id']
                    upload_thread=threading.Thread(target=upload_images, args=(upload_name,))
                    upload_thread.start()
                    break
            
            
        # Here we should:
        # -- keep track of how many sources have gone into active fields
        # -- if more are needed, make them, one mosaic at a time. update when the images are ready
        #    -- steps are: download_image_files_legacy.py (prob about 10 mins)
        #    -- make_images (make_overlays_legacy_scale.py)
        #    NB IMAGEDIR and LOTSS_COMPONENT_CATALOGUE need to be set
        # -- if images are ready to upload, upload them
        # -- later we can also extract results from RGZ
        # Do we want to update the component catalogue?

        '''
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
        '''
        sleep(60)
