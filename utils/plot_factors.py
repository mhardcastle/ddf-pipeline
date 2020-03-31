#!/usr/bin/env python

# Determine scale factors by fitting polynomials to data

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
from scipy.optimize import curve_fit
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib
from fitting_factors import pl,chi2,read_frequencies_fluxes

def fitall(scale,frequencies,fluxes,errors):
    alpha=[]
    cnames=[]
    for name in matplotlib.colors.cnames:
        if 'white' not in name:
            cnames.append(name)
    d,f=fluxes.shape
    for i in range(d):
        sf=np.copy(fluxes[i])
        ef=np.copy(errors[i])
        sf[smask]*=scale
        ef[smask]*=scale
        popt, pcov = curve_fit(pl, frequencies, sf, sigma=ef, p0=[frequencies[0],-1],maxfev=20000)
        chi=chi2(frequencies,sf,ef,*popt)
        sf[ef==1e6]=np.nan
        plt.errorbar(frequencies,sf,yerr=ef,linestyle='none',color=cnames[i])
        plt.plot(frequencies,pl(frequencies,*popt),color=cnames[i])
#        print i,chi,popt[1]
        alpha.append(popt[1])
    return alpha

def run_all(run, name=''):

    global smask
    frequencies,fluxes,errors,smask,data=read_frequencies_fluxes(name+'crossmatch-'+str(run)+'.fits',name=name)

    print('About to plot',len(data),'data points')

    scale=np.load(name+'crossmatch-results-'+str(run)+'.npy')[:,0]
    print('Scaling factors applied are',scale)
    plt.xscale('log')
    plt.yscale('log')
    a=fitall(scale,frequencies,fluxes,errors)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Flux (Jy)')
    plt.show()
    print('mean alpha is',np.mean(a),'error on mean',old_div(np.std(a),np.sqrt(len(a))))
    plt.hist(a,bins=10)
    plt.show()

if __name__=='__main__':
    if len(sys.argv) == 2:
        run_all(int(sys.argv[1]))
    else:
        run_all(int(sys.argv[1]),name=sys.argv[2])
