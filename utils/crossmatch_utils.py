from __future__ import print_function
from __future__ import division
# Tools for cross-matching catalogues

from builtins import str
from builtins import range
from past.utils import old_div
import numpy as np

def bootstrap(data,function,iters):
    result=np.zeros(iters)
    for i in range(iters):
        x=data[np.random.randint(low=0,high=len(data),
                                 size=len(data))]
        result[i]=function(x)
    return result

def separation(c_ra,c_dec,ra,dec):
    # all values in degrees
    return np.sqrt((np.cos(c_dec*np.pi/180.0)*(ra-c_ra))**2.0+(dec-c_dec)**2.0)

def filter_catalogue(t,c_ra,c_dec,radius):
#    r=np.sqrt((np.cos(c_dec*np.pi/180.0)*(t['RA']-c_ra))**2.0+(t['DEC']-c_dec)**2.0)
    r=separation(c_ra,c_dec,t['RA'],t['DEC'])
    return t[r<radius]

def select_isolated_sources(t,radius):
    t['NN_dist']=np.nan
    for r in t:
        dist=3600.0*separation(r['RA'],r['DEC'],t['RA'],t['DEC'])
        # dist=np.sqrt((np.cos(c_dec*np.pi/180.0)*(t['RA']-r['RA']))**2.0+(t['DEC']-r['DEC'])**2.0)*3600.0
        dist.sort()
        r['NN_dist']=dist[1]

    t=t[t['NN_dist']>radius]
    return t

def match_catalogues(t,tab,radius,label,group=None):
    # a replacement for STILTS
    # t is the original catalogue to which we append results from tab labelled by label
    # group is a label which is used in bootstrap: if set then t['g_count_'+str(label)]
    # must exist

    oldv=['Total_flux','E_Total_flux','Peak_flux','E_Peak_flux','RA','DEC']
    newv=['_'+s for s in oldv]
    blankv=['_separation','_dRA','_dDEC']
    for s in newv+blankv:
        t[label+s]=np.nan
    rdeg=radius/3600.0
    minra=np.min(t['RA']-rdeg)
    maxra=np.max(t['RA']+rdeg)
    mindec=np.min(t['DEC']-rdeg)
    maxdec=np.max(t['DEC']+rdeg)
    # pre-filter tab, which may be all-sky
    tab=tab[(tab['RA']>minra) & (tab['RA']<maxra) & (tab['DEC']>mindec) & (tab['DEC']<maxdec)]
    matches=0
    for r in t:
        dist=3600.0*separation(r['RA'],r['DEC'],tab['RA'],tab['DEC'])
        stab=tab[dist<radius]
        df=dist[dist<radius]
        if len(stab)==1:
            # got a unique match
            matches+=1
            for i in range(len(oldv)):
                if oldv[i] in stab.colnames:
                    r[label+newv[i]]=stab[0][oldv[i]]

            r[label+'_separation']=df[0]
            r[label+'_dRA']=3600.0*np.cos(np.pi*r['DEC']/180.0)*(r['RA']-stab[0]['RA'])
            r[label+'_dDEC']=3600.0*(r['DEC']-stab[0]['DEC'])

            if group is not None:
                r['g_count_'+str(group)]+=1

    return matches
                
if __name__=='__main__':
    from astropy.table import Table

    t=Table.read('image_full_ampphase1m.smooth.int.restored.pybdsm.srl.fits')
    ref_ra=np.mean(t['RA'])
    ref_dec=np.mean(t['DEC'])
    print('original length:',len(t))
    t=filter_catalogue(t,ref_ra,ref_dec,2.5)
    print('filter to 2.5 deg:',len(t))
    t=select_isolated_sources(t,30)
    print('isolated sources',len(t))
    t=t[t['Total_flux']>8e-3]
    nvss=Table.read('/stri-data/mjh/bootstrap/NVSS.fits')
    nvss=filter_catalogue(nvss,ref_ra,ref_dec,2.5)
    print('There are',len(t),'sources and',len(nvss),'NVSS sources')
    t.write('lofar_in.fits',overwrite=True)
    nvss.write('nvss_in.fits',overwrite=True)
    match_catalogues(t,nvss,10,'NVSS')
    t=t[~np.isnan(t['NVSS_separation'])]
    ratios=old_div(t['Total_flux'],t['NVSS_Total_flux'])
    print('Median flux ratio:',np.median(ratios))
    print('Bootstrap range:',np.percentile(bootstrap(ratios,np.median,1000),(16,84)))
    t.write('matched.fits',overwrite=True)
