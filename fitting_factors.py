#!/usr/bin/python

# Determine scale factors by fitting power laws to data
# This has to be run as a standalone code since it may use MPI

from scipy.optimize import curve_fit
from astropy.table import Table
import numpy as np
import emcee
#import corner
#import matplotlib.pyplot as plt
import sys

def check_mpi():
    have_mpi=False
    rank=0
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank=comm.Get_rank()
        size=comm.Get_size()
        if size>1:
            if rank==0:
                print 'Using MPI with',size,'threads'
            have_mpi=True
        else:
            print 'MPI is present but inactive'
            have_mpi=False
    except:
        print 'MPI not present or not working'
    return have_mpi,rank
    
def pl(freq,norm,alpha):
     normfreq=140e6
     return norm*(freq/normfreq)**alpha

def chi2(freq,data,errors,norm,alpha):
    return np.sum(((data-pl(freq,norm,alpha))**2.0)/(2*errors)**2.0)

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
            print 'Caught maxfev error:',scale
            cs+=1e6
#        print i,pcov[1],chi
    return -0.5*cs-d*np.sum(np.log(scale))

def lnprior(X):
    for s in X:
        if s < 0.5 or s > 2.0:
            return -np.inf
    return np.sum(-np.log(X))

def lnpost(scale,x,y,yerr):
    return lnprior(scale)+lnlike(scale,frequencies,fluxes,errors)

def read_frequencies_fluxes(intable):
    lines=open('frequencies.txt').readlines()
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
            if np.isnan(r[k]):
                fluxes[i,j]=1.0
                errors[i,j]=1e6
            else:
                fluxes[i,j]=r[k]
                errors[i,j]=r[k_e]

    return frequencies,fluxes,errors,smask,data

def run_all(run):

    global fluxes
    global errors
    global frequencies
    global smask

    frequencies,fluxes,errors,smask,data=read_frequencies_fluxes('crossmatch-'+str(run)+'.fits')

    have_mpi,rank=check_mpi()

    if rank==0:
        print 'About to fit to',len(data),'data points'

    nwalkers=24
    ndim=np.sum(smask)
    if rank==0:
        print 'Fitting',ndim,'scale factors'
    if run>1:
        oscale=np.load('crossmatch-results-'+str(run-1)+'.npy')[:,0]
        scale=[oscale+np.random.normal(loc=0.0,scale=0.02,size=ndim)
               for i in range(nwalkers)]
    else:
        scale=[np.abs(np.random.normal(loc=1.0,scale=0.1,size=ndim))
               for i in range(nwalkers)]

    if have_mpi:
        threads=None
        from emcee.utils import MPIPool
        pool = MPIPool()
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    else:
        pool=None

    # run MCMC
    print 'run the sampler, pool is',pool
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost,
                                    args=(frequencies,fluxes,errors), pool=pool)
    sampler.run_mcmc(scale, 1000)

    if pool:
        pool.close()

    samples=sampler.chain[:, 400:, :].reshape((-1, ndim))
    samplest=samples.transpose()

    means=np.mean(samplest,axis=1)
    errors=np.percentile(samplest,(16,84),axis=1)-means

    for i in range(ndim):
        print i,means[i],errors[0,i],errors[1,i]

    output=np.vstack((means,errors)).T
    np.save('crossmatch-results-'+str(run)+'.npy',output)

    # plot the sampler chain
#    for i in range(ndim):
#        plt.subplot(ndim,1,i+1)
#        plt.plot(sampler.chain[:,:,i].transpose())
#        plt.ylabel(str(i))
#    plt.xlabel('Samples')
#    plt.savefig('walkers-'+str(run)+'.pdf')

    # make triangle plot
#    fig = corner.corner(samples)
#    plt.savefig('corner-'+str(run)+'.pdf')

if __name__=='__main__':
    run_all(int(sys.argv[1]))
