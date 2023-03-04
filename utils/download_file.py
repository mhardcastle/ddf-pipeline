#!/usr/bin/python
# Generic file download with retry and check for length

from __future__ import print_function
import requests
from tqdm import tqdm
import os
from time import sleep,time

def download_file(url,filename,catch_codes=(),retry_interval=60,retry_partial=False,selected_range=None,force_retry_partial=False,chunk_size=8192,progress_bar=False):
    '''Download a file from URL url to file filename.  Optionally, specify
    a tuple of HTTP response codes in catch_codes where we will back
    off and retry rather than failing.

    retry_interval is the time in seconds to wait before retry if one
    of these codes is received.

    if retry_partial is set, then an incomplete file prompts a retry
    with HTTP Range keyword (should work on most servers). Otherwise
    incomplete files are retried from start.

    if force_retry_partial is set, then the target file must already
    exist, and the code will attempt to complete the download.

    if progress_bar is set then a TQDM progress bar is displayed (not
    the default since normally this is not used interactively.

    '''
    
    downloaded=False
    retrying_partial=force_retry_partial
    if force_retry_partial:
         fsize=os.path.getsize(filename)
         selected_range=(fsize,)
    while not downloaded:
        connected=False
        while not connected:
            try:
                print('Opening connection to',url)
                if selected_range is not None:
                    use_range=True
                    if len(selected_range)==2:
                        headers={'Range':'bytes=%i-%i' % selected_range}
                    else:
                        headers={'Range':'bytes=%i-' % selected_range}
                    expected_response=206
                else:
                    use_range=False
                    headers=None
                    expected_response=200
                response = requests.get(url, stream=True,verify=True,timeout=60,headers=headers)
                if response.status_code!=expected_response:
                    print('Unexpected response code received!')
                    print(response.headers)
                    if response.status_code in catch_codes:
                        print('Retrying in %i seconds' % retry_interval)
                        sleep(retry_interval)
                        continue
                    else:
                        raise RuntimeError('Download failed, code was %i' % response.status_code)
                if not use_range:
                    esize=int(response.headers['Content-Length'])
                    psize=esize
                else:
                    cr=response.headers['Content-Range'].split()[1]
                    bits=cr.split('/')
                    esize=int(bits[1])
                    pmin,pmax=[int(v) for v in bits[0].split('-')]
                    psize=pmax-pmin+1
            except requests.exceptions.ConnectionError:
                print('Connection error! sleeping 30 seconds before retry...')
                sleep(30)
            except (requests.exceptions.Timeout,requests.exceptions.ReadTimeout):
                print('Timeout! sleeping 30 seconds before retry...')
                sleep(30)
            else:
                connected=True
        try:
            print('Downloading %i bytes from a total of %i' % (psize,esize))
            starttime=time()
            if retrying_partial:
                mode='ab'
            else:
                mode='wb'
            with open(filename, mode) as fd:
                if progress_bar:
                    for chunk in tqdm(response.iter_content(chunk_size=chunk_size),total=int(psize/chunk_size),unit_scale=chunk_size*1.0/1024/1024,unit='MB'):
                        if chunk: fd.write(chunk)
                    print()
                else:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        if chunk: fd.write(chunk)
            fsize=os.path.getsize(filename)
            if use_range and not retrying_partial and psize!=fsize:
                # a range was specified as an argument: in this case we only retry if the file is not the same as the specified partial size
                print('Partial download incomplete (expected %i, got %i)! Retrying' % (psize, fsize))
            elif (not(use_range) or retrying_partial) and esize!=fsize:
                print('Download incomplete (expected %i, got %i)! Retrying' % (esize, fsize))
                if retry_partial:
                    retrying_partial=True
                    selected_range=(fsize,)
                    print('Retry partial range is:',selected_range)
            else:
                endtime=time()
                dt=endtime-starttime
                print('Download successful, %i of %i bytes received in %.2f seconds (%.2f MB/s)' % (fsize, esize, dt, fsize/(dt*1024*1024)))
                downloaded=True

        except (requests.exceptions.ConnectionError,requests.exceptions.Timeout,requests.exceptions.ChunkedEncodingError):
            print('Connection error! sleeping 30 seconds before retry...')
            if retry_partial:
                fsize=os.path.getsize(filename)
                retrying_partial=True
                selected_range=(fsize,)
                print('Retry partial range is:',selected_range)

            sleep(30) # back to the connection

    del response
    return downloaded

if __name__=='__main__':
    import sys
    download_file(sys.argv[1],sys.argv[2])
    
