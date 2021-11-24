# Wrapper functions for the science data repository interaction
# Because obviously we need another way of downloading files...

from __future__ import print_function
import requests
import os
import sys
import json
from download_file import download_file
from time import sleep

class SDR(object):

    def __init__(self,url='https://tdr-test.surfsara.nl/api/objects/lotss-dr2/',target='.'):
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
            
    def get_status(self,field):
        field=self.fc(field)
        r=requests.get(self.url+field+'/status'+self.tokenstr)
        if r.status_code!=200:
            raise RuntimeError('Failed to get status! (Probably non-existent field')
        else:
            j=json.loads(r.text)
            files={}
            for i,d in enumerate(j):
                files[d['name']]=(d['status'],i)
            return files

    def stage(self,field,filename):
        field=self.fc(field)
        files=self.get_status(field)
        if filename not in files:
            raise RuntimeError('Staging nonexistent file')
        status,number=files[filename]
        if status=='DUL':
            return 'Online'
        else:
            r=requests.post(self.url+field+'/stage/'+str(number)+self.tokenstr)
            return 'Staged'
        
    def download(self,field,filename):
        field=self.fc(field)
        files=self.get_status(field)
        if filename not in files:
            raise RuntimeError('File not found')
        status,number=files[filename]
        if status!='DUL':
            raise RuntimeError('File not online!')
        download_file(self.url+field+'/files/'+filename+self.tokenstr,self.target+'/'+filename)
        
    def download_and_stage(self,field,filenames):
        field=self.fc(field)
        files=self.get_status(field)
        for f in filenames:
            if f not in files:
                raise RuntimeError('File not found!')
            else:
                if files[f][0]=='OFL':
                    self.stage(field,f)

        print('Waiting for files to be online:')

        while True:
            files=self.get_status(field)
            count=0
            for f in filenames:
                if files[f][0]=='DUL':
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
            self.download(field,f)
            
