from __future__ import print_function
from __future__ import absolute_import
# Get an external catalogue by downloading

from builtins import str
from hextile import hextile
import requests
import astropy.coordinates as coord
import astropy.units as u
import os
from time import sleep
from download_file import download_file
from auxcodes import find_fullres_image

CSIZE=0.5
PSBASE='/data/lofar/panstarrs/healpix'

def tile(file):
    return hextile(file,CSIZE*0.9)

def download_required(method):
    if method=='pslocal':
        #pslocal makes the merged file directly
        return ~os.path.isfile(method+'/'+method+'.txt')
            
    ra_factor,pos=tile(find_fullres_image())

    for i,p in enumerate(pos):
        outfile=method+'/'+method+'-'+str(i)+'.vo'
        if not os.path.isfile(outfile):
            return True
    return False

def get_cat(method,retries=100):

    cwd=os.getcwd()
    try:
        os.mkdir(method)
    except OSError:
        pass

    if method=='pslocal':
        hplist=[]
    
    if method=='wise':
        from astroquery.irsa import Irsa
        Irsa.ROW_LIMIT=1000000

    ra_factor,pos=tile(find_fullres_image())
    print('Downloading catalogues for',len(pos),'sky positions')
    for i,p in enumerate(pos):
        outfile=method+'/'+method+'-'+str(i)+'.vo'
        if os.path.isfile(outfile):
            print('Catalogue at position',p,'already present')
            continue
        print('Downloading at position',p)
        if method=='panstarrs':
            count=0
            while True:
                try:
                    r = requests.post('http://archive.stsci.edu/panstarrs/search.php', data = {'ra':p[0],'dec':p[1],'SR':CSIZE,'max_records':100000,'nDetections':">+5",'action':'Search','selectedColumnsCsv':'objid,ramean,decmean'},timeout=300)
                except requests.exceptions.Timeout:
                    print('Timeout, retrying!')
                else:
                    if 'Warning' not in r.text and 'Please' not in r.text:
                        break
                    else:
                        # will go round the loop again
                        print('Bad response, retry download (%i)' % count)
                        sleep(5+count*15)
                count+=1
                if count>=retries:
                    raise RuntimeError('Number of retries exceeded for download')
                        
            f=open(outfile,'w')
            f.writelines(r.text)
            f.close()
        elif method=='wise':
            t=Irsa.query_region(coord.SkyCoord(p[0],p[1],unit=(u.deg,u.deg)), catalog='allwise_p3as_psd', radius='0d30m0s')
            t.write(outfile,format='votable')
        elif method=='pslocal':
            from astropy_healpix import HEALPix
            hp = HEALPix(nside=64)
            cs = hp.cone_search_lonlat(p[0]*u.deg, p[1]*u.deg, radius=CSIZE*u.deg)
            hplist += list(cs)
            if not os.path.isdir(PSBASE):
                # we don't have a local PS database, so download
                for pix in cs:
                    outfile=method+'/'+str(pix)
                    if not os.path.isfile(outfile):
                        print('Downloading healpix pixel',pix)
                        download_file('http://uhhpc.herts.ac.uk/panstarrs-healpix/'+str(pix),outfile)
        else:
            raise NotImplementedError('Method '+method)
    if method=='pslocal':
        hplist=list(set(hplist))
        print('Found',len(hplist),'unique healpix pixels')
        outname=method+'/'+method+'.txt'
        with open(outname,'w') as outfile:
            outfile.write('# RA DEC ObjID\n')
        for pixel in hplist:
            print('Appending pixel',pixel)
            if os.path.isdir(PSBASE):
                pixelfile=PSBASE+'/'+str(pixel)
            else:
                pixelfile=method+'/'+str(pixel)
            if not os.path.isfile(pixelfile):
                raise RuntimeError('Pixel file '+pixelfile+'does not exist')
            os.system('cat '+pixelfile+' >> '+outname)
        
if __name__=='__main__':
    import sys
    get_cat(sys.argv[1])
