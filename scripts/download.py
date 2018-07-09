#!/usr/bin/env python

# Grab all the files from a SKSP URL

import requests
from lxml import html
import shutil
import os.path
from time import sleep
from surveys_db import SurveysDB,use_database,tag_idd
from auxcodes import die

def download_dataset(server,root):
    page=requests.get(server+root,verify=False)
    print page.status_code
    print page.headers['content-type']
    tree=html.fromstring(page.text)
    row = tree.xpath('//a')
    files=[]
    urls=[]
    for r in row:
        if 'title' in r.attrib and 'Download' in r.attrib['title']:
            files.append(r.attrib['download'])
            urls.append(r.attrib['href'].replace('../..',''))
    if len(files)!=25:
        print 'There should be 25 files but there are only %s! Check SARA manually.'%len(files)
        return False
    for f,u in zip(files,urls):
        if os.path.isfile(f):
            print 'File',f,'already exists, skipping'
        else:
            print 'Downloading',f
            url=server+u
            print url
            while True:
                try:
                    response = requests.get(url, stream=True,verify=False,timeout=300)
                except requests.exceptions.ConnectionError,requests.exceptions.Timeout:
                    print 'Connection error! sleeping 60 seconds before retry...'
                    sleep(60)
                else:
                    break
            with open(f, 'wb') as out_file:
                shutil.copyfileobj(response.raw, out_file)
            del response
    return True

def download_db_create(name):
    sdb=SurveysDB()
    id=int(name[1:]) # get the L out
    if sdb.get_id(id):
        die('ID to be downloaded already exists in database!')

    idd=sdb.create_id(id)
    idd['status']='Downloading'
    tag_idd(sdb,idd)
    sdb.set_id(idd)
    sdb.close()

def download_db_update(name,status):
    sdb=SurveysDB()
    id=int(name[1:])
    idd=sdb.get_id(id)
    idd['status']='Downloaded' if status else 'D/L failed'
    sdb.set_id(idd)
    sdb.close()
    
if __name__=='__main__':

    import sys
    name=sys.argv[1]
    try:
        os.mkdir(name)
    except OSError:
        pass
    os.chdir(name)
    if use_database:
        download_db_create(name)
    
    status=download_dataset('https://lofar-webdav.grid.sara.nl','/SKSP/'+name+'/')

    if use_database:
        download_db_update(name,status)

