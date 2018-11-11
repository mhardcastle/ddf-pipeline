#!/usr/bin/env python

# Grab all the files from a SKSP URL

import requests
from lxml import html
import shutil
import os.path
from time import sleep
import sys

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
        if 'title' in r.attrib and 'Download' in r.attrib['title'] and 'GSM' in r.attrib['download']:
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
            downloaded=False
            while not downloaded:
                connected=False
                while not connected:
                    try:
                        print 'Opening connection'
                        response = requests.get(url, stream=True,verify=False,timeout=60)
                        if response.status_code!=200:
                            print response.headers
                            raise RuntimeError('Code was %i' % response.status_code)
                        esize=long(response.headers['Content-Length'])
                    except requests.exceptions.ConnectionError,requests.exceptions.Timeout:
                        print 'Connection error! sleeping 30 seconds before retry...'
                        sleep(30)
                    else:
                        connected=True
                try:
                    print 'Downloading'
                    with open(filename, 'wb') as fd:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                fd.write(chunk)
                    fsize=os.path.getsize(filename)
                    if esize!=fsize:
                        print 'Download incomplete (expected %i, got %i)! Retrying' % (esize, fsize)
                    else:
                        print 'Download successful'
                        downloaded=True
                        
                except requests.exceptions.ConnectionError,requests.exceptions.Timeout:
                    print 'Connection error! sleeping 30 seconds before retry...'
                    sleep(30) # back to the connection


            del response
            
    return True
    
if __name__=='__main__':

    import sys
    name=sys.argv[1]
    try:
        os.mkdir(name)
    except OSError:
        pass
    
    status=download_dataset('https://lofar-webdav.grid.sara.nl','/SKSP/'+name+'/',workdir='./'+name)

