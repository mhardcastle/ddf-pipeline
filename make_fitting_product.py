#!/usr/bin/python

# Make the catalogue needed for the MCMC scaling factor fitting

from astropy.table import Table
import numpy as np

def filter_catalogue(t,c_ra,c_dec,radius):
    r=np.sqrt((np.cos(c_dec*np.pi/180.0)*(t['RA']-c_ra))**2.0+(t['DEC']-c_dec)**2.0)
    return t[r<radius]

def make_catalogue(name,c_ra,c_dec,radius,cats):

    # cats needs to be a list of catalogues with a filename, short
    # name, group ID and matching radius in arcsec.
    # group IDs are a set of objects -- can be anything -- such that we require at least one flux from each group. 
    # Each catalogue needs RA, DEC, Total_flux and E_Total_flux.

    t=Table.read(name,format='ascii.commented_header',header_start=-1)
    print 'Total table length is',len(t)
    t=filter_catalogue(t,c_ra,c_dec,radius)
    print 'Filtered within',radius,'degrees:',len(t)
    t=t[t['Total_flux']>0.15]
    print 'Bright sources:',len(t)
    t['NN_dist']=np.nan

    # Filter for isolated sources
    for r in t:
        dist=np.sqrt((np.cos(c_dec*np.pi/180.0)*(t['RA']-r['RA']))**2.0+(t['DEC']-r['DEC'])**2.0)*3600.0
        dist.sort()
        r['NN_dist']=dist[1]

    t=t[t['NN_dist']>100]
    print 'Remove close neighbours:',len(t)

    ctab=[]
    groups=[]
    for n,sh,group,cmrad in cats:
        tab=Table.read(n)
        ctab.append(filter_catalogue(tab,c_ra,c_dec,radius))
        groups.append(group)
        print 'Table',sh,'has',len(ctab[-1]),'entries'

    groups=set(groups)
    # now do cross-matching
    for g in groups:
        t['g_count_'+str(g)]=0
    for i,(n,sh,group,cmrad) in enumerate(cats):
        tab=ctab[i]
        t[sh+'_flux']=np.nan
        t[sh+'_e_flux']=np.nan
        for r in t:
            dist=np.sqrt((np.cos(c_dec*np.pi/180.0)*(tab['RA']-r['RA']))**2.0+(tab['DEC']-r['DEC'])**2.0)*3600.0
            stab=tab[dist<cmrad]
            if len(stab)==1:
                # got a unique match
                r[sh+'_flux']=stab[0]['Total_flux']
                r[sh+'_e_flux']=stab[0]['E_Total_flux']
                r['g_count_'+str(group)]+=1
    # Now reject sources that have no match in a given group
    for g in groups:
        t=t[t['g_count_'+str(g)]>0]

    t.write('crossmatch-1.fits',overwrite=True)
                          

if __name__=='__main__':
    cats=[['/stri-data/mjh/bootstrap/VLSS.fits','VLSS',1,40.0],
          ['/stri-data/mjh/bootstrap/wenss.fits','WENSS',2,10.0]]
    make_catalogue('cube.pybdsm.srl',209.56,54.92,3.0,cats)
