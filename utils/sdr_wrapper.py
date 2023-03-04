# Wrapper functions for the science data repository interaction
# Because obviously we need another way of downloading files...

from __future__ import print_function
import requests
import os
import sys
import json
from download_file import download_file
from time import sleep
import numpy as np

def request_with_retry(url,rfunction=requests.get):
    done=False
    while not done:
        try:
            r=rfunction(url,timeout=60)
        except requests.exceptions.ReadTimeout:
            print('Timeout fetching URL ',url,': retrying')
            # fuzzy sleep to avoid synchronized retries from multiple clients
            sleep(abs(np.random.normal(30,5)))
        else:
            done=True
    return r

class SDR(object):

    def __init__(self,url='https://repository.surfsara.nl/api/objects/lotss-dr2/',target='.'):
        if 'SDR_URL' in os.environ:
            url=os.environ['SDR_URL']
        self.url=url
        if 'SDR_TOKEN' in os.environ:
            self.token=os.environ['SDR_TOKEN']
        else:
            self.token=None

        if self.token is not None:
            self.tokenstr='?share-token='+self.token
        else:
            self.tokenstr=""
        self.target=target

    def fc(self,field):
        ''' Convert a field name from LoTSS format '''
        return field.replace('+','-')

    def get_files(self,field):
        # return mapping of file name to number
        field=self.fc(field)
        r=request_with_retry(self.url+field+'/files'+self.tokenstr)
        if r.status_code!=200:
            raise RuntimeError('Failed to get status! (Probably non-existent field')
        else:
            j=json.loads(r.text)
            files={}
            for i,d in enumerate(j):
                files[d['name']]=i
            return files
        
    def get_status(self,field):
        # return staging status of all files
        field=self.fc(field)
        r=request_with_retry(self.url+field+'/status'+self.tokenstr)
        if r.status_code!=200:
            raise RuntimeError('Failed to get status! (Probably non-existent field')
        else:
            j=json.loads(r.text)
            files={}
            for d in j:
                files[d['name']]=d['status']
            return files

    def stage(self,field,filename):
        field=self.fc(field)
        files=self.get_files(field)
        status=self.get_status(field)
        if filename not in files or filename not in status:
            raise RuntimeError('Staging nonexistent file')
        status=status[filename]
        number=files[filename]
        if status=='DUL':
            return 'Online'
        else:
            r=request_with_retry(self.url+field+'/stage/'+str(number)+self.tokenstr,rfunction=requests.post)
            return 'Staged'
        
    def download(self,field,filename,progress_bar=False):
        field=self.fc(field)
        files=self.get_status(field)
        if filename not in files:
            raise RuntimeError('File not found')
        status=files[filename]
        if status!='DUL':
            raise RuntimeError('File not online!')
        download_file(self.url+field+'/files/'+filename+self.tokenstr,self.target+'/'+filename,catch_codes=(500,),retry_partial=True,progress_bar=progress_bar)
        
    def download_and_stage(self,field,filenames,progress_bar=False):
        print('Progress bar is',progress_bar)
        field=self.fc(field)
        files=self.get_status(field)
        for f in filenames:
            if f not in files:
                raise RuntimeError('File not found!')
            else:
                if files[f]=='OFL':
                    self.stage(field,f)

        print('Waiting for files to be online:')

        while True:
            files=self.get_status(field)
            count=0
            for f in filenames:
                if files[f]=='DUL':
                    count+=1
            print('%i/%i... ' % (count,len(filenames)),end='')
            sys.stdout.flush()
            if count==len(filenames):
                print()
                break
            else:
                sleep(30)
            
        for f in filenames:
            print('Downloading',f)
            self.download(field,f,progress_bar=progress_bar)
            
