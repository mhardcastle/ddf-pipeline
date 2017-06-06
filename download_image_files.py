# formerly get_panstarrs_file.py . Will be amended to allow downloads
# of images from WISE as well.

import requests
from astropy.io import fits
from astropy.table import Table
import shutil
import os.path
import sys

def download_file(url,outname):
    if os.path.isfile(outname):
        print 'File',outname,'already exists, skipping'
    else:
        print 'Downloading',outname
        while True:
            try:
                response = requests.get(url, stream=True,verify=False)
            except requests.exceptions.ConnectionError:
                print 'Connection error! sleeping 60 seconds before retry...'
                sleep(60)
            else:
                break
        with open(outname, 'wb') as out_file:
            shutil.copyfileobj(response.raw, out_file)
        del response

def get_panstarrs(ra,dec,psband):
    page=requests.get('http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra=%f&dec=%f' % (ra,dec),verify=False)
    print page.status_code
    print page.headers['content-type']
    lines=page.text.split('\n')
    downloads=[]
    for l in lines[1:]:
        bits=l.split()
        if bits[4] in psband or psband=='*':
            download_file('http://ps1images.stsci.edu/'+bits[7],bits[8])
            downloads.append(bits[8])
    return downloads

if __name__=='__main__':
    outfile=open('panstarrs-list.txt','w')
    t=Table.read(sys.argv[1])
    for r in t:
        psnames=get_panstarrs(r['RA'],r['DEC'],'i')
        print >>outfile,r['Source_id'],psnames[0]
    outfile.close()
