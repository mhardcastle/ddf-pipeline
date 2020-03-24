#!/usr/bin/python

# Make the catalogue needed for the MCMC scaling factor fitting

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from astropy.table import Table
import numpy as np
from crossmatch_utils import filter_catalogue,select_isolated_sources,match_catalogues

def make_catalogue(name,c_ra,c_dec,radius,cats,outnameprefix=''):

    # cats needs to be a list of catalogues with a filename, short
    # name, group ID and matching radius in arcsec.
    # group IDs are a set of objects -- can be anything -- such that we require at least one flux from each group. 
    # Each catalogue needs RA, DEC, Total_flux and E_Total_flux.

    t=Table.read(name,format='ascii.commented_header',header_start=-1)
    print('Total table length is',len(t))
    if len(t)==0:
        raise RuntimeError('No sources in table from pybdsm')
    t=filter_catalogue(t,c_ra,c_dec,radius)
    print('Filtered within',radius,'degrees:',len(t))
    if len(t)==0:
        raise RuntimeError('No sources in central part of image')
    t=t[t['Total_flux']>0.15]
    print('Bright sources:',len(t))
    if len(t)==0:
        raise RuntimeError('No bright sources for crossmatching')

    # Filter for isolated sources
    t=select_isolated_sources(t,100)
    print('Remove close neighbours:',len(t))
    if len(t)==0:
        raise RuntimeError('No sources in table before crossmatching')

    ctab=[]
    groups=[]
    for n,sh,group,cmrad in cats:
        tab=Table.read(n)
        ctab.append(filter_catalogue(tab,c_ra,c_dec,radius))
        groups.append(group)
        print('Table',sh,'has',len(ctab[-1]),'entries')

    groups=set(groups)
    # now do cross-matching
    for g in groups:
        t['g_count_'+str(g)]=0
    for i,(n,sh,group,cmrad) in enumerate(cats):
        tab=ctab[i]
        match_catalogues(t,tab,cmrad,sh,group=group)
#        
#        t[sh+'_flux']=np.nan
#        t[sh+'_e_flux']=np.nan
#        for r in t:
#            dist=np.sqrt((np.cos(c_dec*np.pi/180.0)*(tab['RA']-r['RA']))**2.0+(tab['DEC']-r['DEC'])**2.0)*3600.0
#            stab=tab[dist<cmrad]
#            if len(stab)==1:
#                # got a unique match
#                r[sh+'_flux']=stab[0]['Total_flux']
#                r[sh+'_e_flux']=stab[0]['E_Total_flux']
#                r['g_count_'+str(group)]+=1
    # Now reject sources that have no match in a given group
    for g in groups:
        t=t[t['g_count_'+str(g)]>0]

    if len(t)==0:
        raise RuntimeError('No crossmatches exist after group matching')
    t.write(outnameprefix+'crossmatch-1.fits',overwrite=True)
                          

if __name__=='__main__':
    cats=[['/stri-data/mjh/bootstrap/VLSS.fits','VLSS',1,40.0],
          ['/stri-data/mjh/bootstrap/wenss.fits','WENSS',2,10.0]]
    make_catalogue('cube.pybdsm.srl',209.56,54.92,3.0,cats)
