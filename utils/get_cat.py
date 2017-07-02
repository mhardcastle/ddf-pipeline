# Get an external catalogue by downloading

from hextile import hextile
import requests
import astropy.coordinates as coord
import astropy.units as u
import os
CSIZE=0.5

def tile(file):
    return hextile(file,CSIZE*0.9)

def download_required(method):
    ra_factor,pos=tile('image_ampphase1.app.restored.fits')
    for i,p in enumerate(pos):
        outfile=method+'/'+method+'-'+str(i)+'.vo'
        if not os.path.isfile(outfile):
            return True
    return False

def get_cat(method):

    cwd=os.getcwd()
    try:
        os.mkdir(method)
    except OSError:
        pass

    if method=='wise':
        from astroquery.irsa import Irsa
        Irsa.ROW_LIMIT=1000000

    ra_factor,pos=tile(cwd+'/image_ampphase1.app.restored.fits')
    print 'Downloading catalogues for',len(pos),'sky positions'
    for i,p in enumerate(pos):
        outfile=method+'/'+method+'-'+str(i)+'.vo'
        if os.path.isfile(outfile):
            print 'Catalogue at position',p,'already present'
            continue
        print 'Downloading at position',p
        if method=='panstarrs':
            while True:
                try:
                    r = requests.post('http://archive.stsci.edu/panstarrs/search.php', data = {'ra':p[0],'dec':p[1],'SR':CSIZE,'max_records':100000,'nDetections':">+5",'action':'Search','selectedColumnsCsv':'objid,ramean,decmean'},timeout=300)
                except requests.exceptions.Timeout:
                    print 'Timeout, retrying!'
                else:
                    break
            f=open(outfile,'w')
            f.writelines(r.text)
            f.close()
        elif method=='wise':
            t=Irsa.query_region(coord.SkyCoord(p[0],p[1],unit=(u.deg,u.deg)), catalog='allwise_p3as_psd', radius='0d30m0s')
            t.write(outfile,format='votable')
        else:
            raise NotImplementedError('Method '+method)
        
if __name__=='__main__':
    import sys
    get_cat(sys.argv[1])
