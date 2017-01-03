#!/usr/bin/env python

# Grab all the files from a SKSP URL

import requests
from lxml import html
import shutil
import os.path
from time import sleep

def download_dataset(server,root):
    page=requests.get(server+root,verify=False)
    print page.status_code
    print page.headers['content-type']
    tree=html.fromstring(page.text)
    row = tree.xpath('//a')
    for r in row:
        if r.text is not None and 'CAL' in r.text:
            if os.path.isfile(r.text):
                print 'File',r.text,'already exists, skipping'
            else:
                print 'Downloading',r.text
                url=server+r.attrib['href']
                while True:
                    try:
                        response = requests.get(url, stream=True,verify=False)
                    except requests.exceptions.ConnectionError:
                        print 'Connection error! sleeping 60 seconds before retry...'
                        sleep(60)
                    else:
                        break
                with open(r.text, 'wb') as out_file:
                    shutil.copyfileobj(response.raw, out_file)
                del response
            

if __name__=='__main__':

    import sys
    name=sys.argv[1]
    try:
        os.mkdir(name)
    except OSError:
        pass
    os.chdir(name)
    download_dataset('https://lofar-webdav.grid.sara.nl','/SKSP/'+name+'/')
