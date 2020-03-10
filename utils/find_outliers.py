#!/usr/bin/python

# Determine scale factors by fitting polynomials to data

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from scipy.optimize import curve_fit
#from astropy.table import Table
import numpy as np
import sys
from fitting_factors import pl,chi2,read_frequencies_fluxes

def fitall(scale,frequencies,fluxes,errors,smask):

    alphas=[]
    norms=[]
    chiv=[]
    d,f=fluxes.shape
    for i in range(d):
        sf=np.copy(fluxes[i])
        ef=np.copy(errors[i])
        sf[smask]*=scale
        ef[smask]*=scale
        popt, pcov = curve_fit(pl, frequencies, sf, sigma=ef, p0=[1.0,-1],maxfev=20000)
        chiv.append(chi2(frequencies,sf,ef,*popt))
        alphas.append(popt[1])
        norms.append(popt[0])
    return np.array((norms,alphas,chiv))

def run_all(run, name=''):

    frequencies,fluxes,errors,smask,data=read_frequencies_fluxes(name+'crossmatch-'+str(run)+'.fits',name=name)

    try:
        scale=np.load(name+'crossmatch-results-'+str(run)+'.npy')[:,0]
    except IOError:
        print('Can\'t load results file')
        return False
    print(scale)

    a=fitall(scale,frequencies,fluxes,errors,smask)
    print('Number of sources is',len(data))
    print('Mean flux density is',np.mean(a[0]))
    print('Mean spectral index  is',np.mean(a[1]))
    print('Median chi^2 is',np.median(a[2]))
    print('Mean chi^2 is',np.mean(a[2]))
    print('Max chi^2 is',np.max(a[2]))

    threshold=100
    print('Number of sources rejected',np.sum(a[2]>threshold))
    filtered=data[a[2]<threshold]

    filtered.write(name+'crossmatch-'+str(run+1)+'.fits',overwrite=True)
    return np.sum(a[2]>threshold)

if __name__=='__main__':
    if len(sys.argv) == 2:
        run_all(int(sys.argv[1]))
    else:
        run_all(int(sys.argv[1]),name=sys.argv[2])
