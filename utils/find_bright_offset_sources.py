from __future__ import print_function
from __future__ import division
from builtins import str
from past.utils import old_div
from astropy.io import fits
from astropy.table import Table
from astropy import wcs
from astropy.coordinates import SkyCoord

def find_bright(root='image_full_ampphase_di_m.NS',cutoff=1):

    t=Table.read(root+'.offset_gaul.fits')
    filter=t['Peak_flux']>cutoff
    filter&=t['DC_Maj']<0.003
    t=t[filter]
    print('Found',len(t),'sources')

    i_int=fits.open(root+'.int.restored.fits')
    i_app=fits.open(root+'.app.restored.fits')
    scale=old_div(i_app[0].data,i_int[0].data)
    w=wcs.WCS(i_int[0])
    x,y,_,_=w.wcs_world2pix(t['RA'],t['DEC'],0,0,0)
    app=[]
    for i,r in enumerate(t):
        app.append(r['Peak_flux']*scale[0,0,int(y[i]),int(x[i])])
    t['App_flux']=app
    filter=t['App_flux']>cutoff
    t=t[filter]
    print('Now at',len(t),'sources')
    t.write('dynspec-targets.fits',overwrite=True)
    f=open('brightlist.csv','w')
    sc=SkyCoord(t['RA'],t['DEC'],frame='icrs')
    strings=sc.to_string(style='hmsdms',sep='',precision=2)
    ilt=[]
    for s in strings:
        ilt.append(str('ILTJ'+s).replace(' ','')[:-1])
    for i,r in enumerate(t):
        f.write('%s,%f,%f,Bright LOFAR\n' % (ilt[i],r['RA'],r['DEC']))
    f.close()
    return len(t)>0
        
if __name__=='__main__':
    find_bright()
    
    
