# formerly get_panstarrs_file.py . Will be amended to allow downloads
# of images from WISE as well.

import requests
from astropy.io import fits
from astropy.table import Table
import astropy.coordinates as coord
import astropy.units as u
from astroquery.ibe import Ibe
import shutil
import os.path
import sys
import re

def download_file(url,outname):
    if os.path.isfile(outname):
        print 'File',outname,'already exists, skipping'
    else:
        print 'Downloading',outname
        while True:
            try:
                response = requests.get(url, stream=True,verify=False)
                if response.status_code!=200:
                    print 'Warning, HTML status code',response.status_code
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

def get_wise(ra,dec,band):
    mission='wise'
    dataset='allwise'
    table='p3am_cdd'
    results=Ibe.query_region(coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'),mission=mission,dataset=dataset,table=table)
    url = 'http://irsa.ipac.caltech.edu/ibe/data/'+mission+'/'+dataset+'/'+table+'/'
    params = { 'coadd_id': results[results['band']==band]['coadd_id'][0],
           'band': band }
    params['coaddgrp'] = params['coadd_id'][:2]
    params['coadd_ra'] = params['coadd_id'][:4]
    path = str.format('{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id:s}-w{band:1d}-int-3.fits',**params)
    outname=path.split('/')[-1]
    download_file(url+path,outname)
    return outname

def get_first(ra,dec):
    url="http://archive.stsci.edu/"
    page=requests.get(url+"vlafirst/search.php?RA=%f&DEC=%f&Radius=30.0&action=Search" % (ra,dec),verify=False)
    print page.status_code
    lines=page.text.split('\n')
    path=None
    for l in lines:
        if 'href="/pub' in l:
            bits=l.split()
            for b in bits:
                m=re.search('href=\"(.*)\"',b)
                if m:
                  path=m.group(1)
                  break
    if path:
        outname=path.split('/')[-1]
        download_file(url+path,outname)
        return outname
    else:
        return None

if __name__=='__main__':
    outfile=open('panstarrs-list.txt','w')
    t=Table.read(sys.argv[1])
    for r in t:
        psnames=get_panstarrs(r['RA'],r['DEC'],'i')
        wisename=get_wise(r['RA'],r['DEC'],1)
        print >>outfile,r['Source_id'],psnames[0],wisename
    outfile.close()
