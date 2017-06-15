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
import numpy as np
from lxml import html
import glob

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
    page=requests.get(url+"vlafirst/search.php?RA=%.7f&DEC=%.6f&Radius=30.0&action=Search" % (ra,dec),verify=False)
    print page.status_code

    tree=html.fromstring(page.text)
    table=tree.xpath('//tbody')
    links=[]
    dists=[]
    for row in table[0].getchildren():
        td=row.getchildren()
        links.append(td[0].getchildren()[0].attrib['href'])
        dists.append(float(td[8].text))

    index=np.argmin(dists)
    path=links[index]
   
    outname=path.split('/')[-1]
    download_file(url+path,outname)
    return outname

if __name__=='__main__':
    t=Table.read(sys.argv[1])
    outfile=open(sys.argv[1].replace('.fits','-list.txt'),'w')

    # read the LOFAR map positions
    g=glob.glob('/data/lofar/mjh/hetdex_v3/mosaics/P*')

    files=[]
    ras=[]
    decs=[]
    for d in g:
        file=d+'/mosaic.fits'
        hdu=fits.open(file)
        ras.append(hdu[0].header['CRVAL1'])
        decs.append(hdu[0].header['CRVAL2'])
        files.append(file)
    ras=np.array(ras)
    decs=np.array(decs)

    for r in t:
        dist=np.cos(decs*np.pi/180.0)*(ras-r['RA'])**2.0 + (decs-r['DEC'])**2.0
        i=np.argmin(dist)
        lofarname=files[i]

        psnames=get_panstarrs(r['RA'],r['DEC'],'i')
        wisename=get_wise(r['RA'],r['DEC'],1)
        firstname=get_first(r['RA'],r['DEC'])
        print >>outfile,r['Source_id'],lofarname,psnames[0],wisename,firstname

    outfile.close()
