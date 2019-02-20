#!/usr/bin/env python

# Grab all the files from a SKSP URL

import requests
from lxml import html
import shutil
import os.path
from time import sleep
import sys
from download_file import download_file

def download_dataset(server,root,workdir='.'):
    print server+root
    while True:
        try:
            print 'Downloading index page',server+root
            page=requests.get(server+root,verify=False,timeout=60)
        except requests.exceptions.ConnectionError,requests.exceptions.Timeout:
            print 'Connection error! sleeping 30 seconds before retry...'
            sleep(30)
        else:
            break
    
    print page.status_code
    if page.status_code!=200:
        print page.headers
        return False
    print page.headers['content-type']
    tree=html.fromstring(page.text)
    row = tree.xpath('//a')
    files=[]
    urls=[]
    for r in row:
        if 'title' in r.attrib and 'Download' in r.attrib['title']:
            # and 'GSM' in r.attrib['download']:
            files.append(r.attrib['download'])
            urls.append(r.attrib['href'].replace('../..',''))
    if len(files)<25:
        print 'There should be 25 files but there are only %s! Check SARA manually.'%len(files)
        return False
    else:
        print 'Downloading',len(files),'distinct files'
    for f,u in zip(files,urls):
        if os.path.isfile(workdir+'/'+f):
            print 'File',f,'already exists, skipping'
        else:
            print 'Downloading',f
            url=server+u
            print url
            filename=workdir+'/'+f
            download_file(url,filename)
            
    return True
    
if __name__=='__main__':

    import sys
    name=sys.argv[1]
    try:
        os.mkdir(name)
    except OSError:
        pass
    
    status=download_dataset('https://lofar-webdav.grid.sara.nl','/SKSP/'+name+'/',workdir='./'+name)

