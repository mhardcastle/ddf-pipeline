#!/usr/bin/env python

# Determine scale factors by fitting power laws to data
# This has to be run as a standalone code since it may use MPI

from __future__ import print_function
from __future__ import division
from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div
from scipy.optimize import curve_fit
from astropy.table import Table
import numpy as np
#import corner
#import matplotlib.pyplot as plt
import sys

def pl(freq,norm,alpha):
     normfreq=140e6
     return norm*(old_div(freq,normfreq))**alpha

def chi2(freq,data,errors,norm,alpha):
    return np.sum(old_div(((data-pl(freq,norm,alpha))**2.0),(2*errors)**2.0))

def lnlike(scale,frequencies,fluxes,errors):
    cs=0
    d,f=fluxes.shape
    for i in range(d):
        sf=np.copy(fluxes[i])
        ef=np.copy(errors[i])
        sf[smask]*=scale
        ef[smask]*=scale
        try:
            popt, pcov = curve_fit(pl, frequencies, sf, sigma=ef, p0=[sf[4],-0.8],maxfev=20000)
            chi=chi2(frequencies,sf,ef,*popt)
            cs+=chi
        except RuntimeError:
            print('Caught maxfev error:',scale)
            cs+=1e6
#        print i,pcov[1],chi
    retval=-0.5*cs-d*np.sum(np.log(scale))
    if np.isnan(retval):
         return -np.inf
    else:
         return retval

def lnprior(X):
    for s in X:
        if s < 0.5 or s > 2.0:
            return -np.inf
    return np.sum(-np.log(X))

def lnpost(scale,x,y,yerr):
    return lnprior(scale)+lnlike(scale,frequencies,fluxes,errors)

def read_frequencies_fluxes(intable,name=''):
    lines=open(name+'frequencies.txt').readlines()
    frequencies=[]
    keywords=[]
    e_keywords=[]
    smask=[]
    for l in lines:
        bits=l.split()
        frequencies.append(float(bits[0]))
        keywords.append(bits[1])
        e_keywords.append(bits[2])
        smask.append(bits[3][:4]=='True')

    frequencies=np.array(frequencies)
    smask=np.array(smask)

    data=Table.read(intable)
    nf=len(frequencies)
    fluxes=np.zeros((len(data),nf))
    errors=np.zeros((len(data),nf))
    for i,r in enumerate(data):
        for j,(k,k_e) in enumerate(zip(keywords,e_keywords)):
             fluxes[i,j]=r[k]
             errors[i,j]=r[k_e]
             if np.isnan(fluxes[i,j]):
                  fluxes[i,j]=1.0
                  errors[i,j]=1e6

    return frequencies,fluxes,errors,smask,data

def run_all(run, name=''):

    global fluxes
    global errors
    global frequencies
    global smask
    import emcee

    frequencies,fluxes,errors,smask,data=read_frequencies_fluxes(name+'crossmatch-'+str(run)+'.fits',name=name)

    print('About to fit to',len(data),'data points')

    ndim=np.sum(smask)
    nwalkers=2*(ndim+1)
    if nwalkers<24:
         nwalkers=24
    print('Fitting',ndim,'scale factors')
    if run>1:
        oscale=np.load(name+'crossmatch-results-'+str(run-1)+'.npy')[:,0]
        scale=[oscale+np.random.normal(loc=0.0,scale=0.02,size=ndim)
               for i in range(nwalkers)]
    else:
        scale=[np.abs(np.random.normal(loc=1.0,scale=0.1,size=ndim))
               for i in range(nwalkers)]

    # run MCMC
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost,
                                    args=(frequencies,fluxes,errors))
    sampler.run_mcmc(scale, 1000)

    samples=sampler.chain[:, 400:, :].reshape((-1, ndim))
    samplest=samples.transpose()

    means=np.mean(samplest,axis=1)
    errors=np.percentile(samplest,(16,84),axis=1)-means

    for i in range(ndim):
        print(i,means[i],errors[0,i],errors[1,i])

    output=np.vstack((means,errors)).T
    np.save(name+'crossmatch-results-'+str(run)+'.npy',output)

if __name__=='__main__':
    if len(sys.argv) == 2:
        run_all(int(sys.argv[1]))
    else:
        run_all(int(sys.argv[1]),name=sys.argv[2])
