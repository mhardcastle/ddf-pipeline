#!/usr/bin/python

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div
from astropy.table import Table
from astropy.io import fits
import numpy as np
from scipy.special import gammaln
import emcee
from auxcodes import get_rms

def model(cbins,norm,alpha):
    return (10**norm)*ds*(old_div(cbins,fluxnorm))**(-alpha)

def loglf(cbins,hist,norm,alpha):
    mu=model(cbins,norm,alpha)
    likelihood=hist*np.log(mu)-mu-gammaln(hist)
    return np.sum(likelihood)
    
def lnprior(X):
    if X[0]<0 or X[0]>5 or X[1]<0 or X[1]>4.0:
        return -np.inf
    return old_div((-(X[1]-1.753)**2.0),(0.01**2.0))

def lnpost(X,x,y):
    return lnprior(X)+loglf(x,y,X[0],X[1])

def do_fit_sourcecounts(t=None, rms=None,do_plots=False,sfindarea=17.09):

    global fluxnorm,ds
    if t is None:
        t=Table.read('image_full_ampphase_di_m.NS.cat.fits')
    print('Number of sources in full table is',len(t))
    if rms is None:
        rms=get_rms(fits.open('image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'))

    cutoff=rms*20
    print('Cutoff will be',cutoff,'Jy')
    #cutoff=1e-3
    maxflux=np.max(t['Total_flux'])
    t=t[t['Total_flux']>cutoff]
    print('Number of sources after completeness cut is',len(t))

    bins=np.logspace(np.log10(cutoff),np.log10(maxflux)*1.01,20)
    cbins=0.5*(bins[:-1]+bins[1:])
    ds=bins[1:]-bins[:-1]

    hist,_=np.histogram(t['Total_flux'],bins=bins)

    zval=None
    for i in range(len(hist)):
        if hist[i]==0:
            zval=i
            break

    if zval is not None:
        cbins=cbins[:zval]
        hist=hist[:zval]
        ds=ds[:zval]

    nwalkers=10
    ndim=2
    fluxnorm=0.1 # find normalization at this flux in Jy

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost,
                                    args=(cbins,hist))

    pos=[[1.6,0.86]+0.01*np.random.normal(size=ndim)
         for i in range(nwalkers)]

    sampler.run_mcmc(pos, 1000)
    samples=sampler.chain[:, 200:, :].reshape((-1, ndim))

    samplest=samples.transpose()

    means=np.mean(samplest,axis=1)
    errors=np.percentile(samplest,(16,84),axis=1)-means

    for i in range(ndim):
        print(i,means[i],errors[0,i],errors[1,i])

    fnorm=means[0]
    falpha=means[1]

    C=[3.5142,0.3738,-0.3138,-0.0717,0.0213,0.0097] # Intema+

    fn=fluxnorm
    totfactor=1.0
    for i in range(10):
        lf=np.log10(fn)
        ncn=0
        for i in range(6):
            ncn+=C[i]*lf**i
        print('for %.3f Jy number count norm should be %f' % (fn,10**ncn))
        measured_ncn=old_div(10**fnorm*fluxnorm*(fn**1.5)*3282.8,sfindarea) # check precise area
        print('measured number count norm is',measured_ncn) 
        scale=(old_div(measured_ncn,10**ncn))#**(1.0/1.5)
        print('scaling factor should be',scale)
        totfactor*=scale
        print('total factor is',totfactor)
        fn=old_div(fluxnorm,totfactor)
        print('New flux norm value is',fn)
        if abs(scale-1.0)<1e-4:
            print('Converged, stopping')
            break

    if do_plots:
        import matplotlib.pyplot as plt
        import corner
        
        plt.scatter(cbins,hist)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(cutoff,maxflux)
        plt.ylim(0.5,np.max(hist)*1.3)

        yv=model(cbins,fnorm,falpha)
        plt.plot(cbins,yv)

        fig = corner.corner(samples)
        plt.show()

    return fnorm,falpha,totfactor

def pl(a, b, g, size=1):
    """Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b"""
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)

def test_sourcecounts(niter=100,size=13780,rms=1e-4):
    factors=[]
    scale=0.7
    for i in range(niter):
        print('Iteration',i)
        fluxes=pl(3*rms,1,-0.753,size=size)
        print(fluxes)
        fluxes+=rms*np.random.normal(size=size)
        fluxes*=scale
        t=Table([fluxes],names=('Total_flux',))
        _,_,totfactor=do_fit_sourcecounts(t=t,rms=rms*scale,do_plots=False)
        factors.append(totfactor)
    print(np.mean(factors), np.std(factors))
        
if __name__=='__main__':
    test_sourcecounts()
    #print do_fit_sourcecounts(do_plots=True)
    
