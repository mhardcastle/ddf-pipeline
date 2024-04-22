import glob
from astropy.io import fits
from astropy.table import Table
import numpy as np
import math
from itertools import chain
import random
import matplotlib.pyplot as plt 
from matplotlib import rc
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import scipy
import os,sys
from auxcodes import sepn
# Get all scipy distirubtions as _distn_names
from scipy.stats._continuous_distns import _distn_names
from scipy import stats

#Define various angle conversion factors (multiply to undertake operation)
arcsec2deg=1.0/3600
arcmin2deg=1.0/60
deg2rad=np.pi/180
deg2arcsec = 1.0/arcsec2deg
rad2deg=180.0/np.pi
arcmin2rad=arcmin2deg*deg2rad
arcsec2rad=arcsec2deg*deg2rad
rad2arcmin=1.0/arcmin2rad
rad2arcsec=1.0/arcsec2rad

def envelope_model(t, coeffs):
    #return abs(coeffs[3] + (coeffs[0] + coeffs[3])/(1+(t/coeffs[2])**(coeffs[1])))
    return coeffs[3] + (coeffs[0] + coeffs[3])/(1+(t/coeffs[2])**(coeffs[1]))

def envelope_residuals(coeffs, y, t):
    return y - envelope_model(t, coeffs)


def histedges_equalN(x, nbin):
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

def find_bestfit_scipyfunction(data,binvals):
    # Got from https://stackoverflow.com/questions/6620471/fitting-empirical-distribution-to-theoretical-ones-with-scipy-python
    pltyvals,pltxvals  = np.histogram(data,density=True,bins=binvals)
    pltxvals = (pltxvals + np.roll(pltxvals, -1))[:-1] / 2.0
    skipping = ['levy_stable','erlang','tukeylambda','ncf','nct']
    sses = []
    used_distnames = []
    notskipping = ['norminvgauss'] # Found to often be best (not always but nice to use the same function)
    for distribution in _distn_names:
        if distribution in skipping or distribution not in notskipping:
            continue
        used_distnames.append(distribution)
        print('Trying to fit function',distribution, 'Total number of distributions',len(_distn_names))
        distribution = getattr(stats,distribution)
        fitdistrib = distribution.fit(data)
        
        sse = np.sum(np.power(pltyvals - distribution.pdf(pltxvals,*fitdistrib), 2.0))
        print('SSE',sse)
        if not np.isnan(sse):
            sses.append(sse)
        #plt.plot(binvals,distribution.pdf(binvals,*fitdistrib),alpha=0.2)
    bestfitting = np.where(sses==np.min(sses))[0][0]
    distribution= used_distnames[bestfitting]
    print('Best fitting distribution is a %s function - %s'%(distribution,np.min(sses)))
    return distribution

def filter_sources(opencat,majcut):

    print('Original length %s'%(np.shape(opencat)))
    delete_entries = np.array([],dtype=np.int8)

    # Find sources that will be filtered out for the fitting.

    # Filter based on nearest neighbours 
    delete_entries = np.append(delete_entries,np.where(opencat['Neighbour_15'] > 1.0)[0])
    delete_entries = np.unique(delete_entries)
    print(len(opencat)-len(delete_entries),'after neighbour')
    neighbournum = len(opencat)-len(delete_entries)

    # Filter based on source size type
    delete_entries = np.append(delete_entries,np.where(opencat['S_Code'] != 'S'))
    delete_entries = np.unique(delete_entries)
    print(len(opencat)-len(delete_entries),'after S_Code cut')
    sourcetypenum = len(opencat)-len(delete_entries)

    # Crude filter to remove large sources
    delete_entries = np.append(delete_entries,np.where(opencat['Maj']*deg2arcsec > majcut))
    delete_entries = np.unique(delete_entries)
    print(len(opencat)-len(delete_entries),'after removing large sources cut')
    sourcesizenum = len(opencat)-len(delete_entries)

    # Remove low SNR ones
    delete_entries = np.append(delete_entries,np.where(opencat['Total_flux']/opencat['E_Total_flux'] < 5.0))
    delete_entries = np.unique(delete_entries)
    print(len(opencat)-len(delete_entries),'after removing <5 snr sources sources')

    opencat.remove_rows(delete_entries)
    print('Size of filterred cat',len(opencat))

    return opencat

def find_nearest_neighbours(opencatfile,ra,dec,radius):

    separations = sepn(ra*deg2rad,dec*deg2rad,opencatfile[1].data['RA']*deg2rad,opencatfile[1].data['DEC']*deg2rad)*rad2arcsec

    separations = separations[np.where(separations < radius)]
    
    return len(separations)

def radial_model(t, a,b):
    return a*np.exp(b*t)


def radial_correction(incatname,inimagename,debug=True):
    # Find some things from image
    f = fits.open(inimagename)
    racent = f[0].header['CRVAL1']
    deccent = f[0].header['CRVAL2']
    beammaj = f[0].header['BMAJ']

    # Add nearest neighbours to cat.
    cat = fits.open(incatname)
    nearestneighbour = np.zeros(len(cat[1].data))
    cols = []
    cols.append(
        fits.Column(name='Neighbour_15', format='F8', array=nearestneighbour)
    )
    orig_cols = cat[1].columns
    new_cols = fits.ColDefs(cols)
    hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
    nearcat = incatname.replace('.fits','_wnear.fits')
    hdu.writeto(nearcat,overwrite=True)
    cat = fits.open(nearcat)
    for i in range(0,len(cat[1].data)):
        ra,dec = cat[1].data['RA'][i],cat[1].data['DEC'][i]
        cat[1].data['Neighbour_15'][i] = find_nearest_neighbours(cat,ra,dec,15.0)
    cat.writeto(nearcat,clobber=True)
    print('Created nearest neighbour cat')
    cat.close()

    # Plot radial correction
    opencat = Table.read(nearcat)
    opencat = filter_sources(opencat,2*beammaj*deg2arcsec)
    separation = sepn(opencat['RA']*deg2rad,opencat['DEC']*deg2rad,racent*deg2rad,deccent*deg2rad)*rad2deg
    radbins = np.arange(0.0,5.0,0.2)
    xvals = []
    yvals = []
    for i in range(0,len(radbins)-1):
        binmin = radbins[i]
        binmax = radbins[i+1]
        sourcemin = np.where(separation > binmin)
        sourcemax = np.where(separation < binmax)
        meetscriterea = np.intersect1d(sourcemin,sourcemax)  
        print('Bin: %s deg - %s deg contains %s sources'%(binmin,binmax,len(meetscriterea)))
        medval = np.median(opencat['Total_flux'][meetscriterea]/opencat['Peak_flux'][meetscriterea])
        plt.plot((binmin+binmax)/2.0,medval,'g+',markersize=10)
        plt.plot(separation[meetscriterea],opencat['Total_flux'][meetscriterea]/opencat['Peak_flux'][meetscriterea],'b.',alpha=0.1)
        xvals.append((binmin+binmax)/2.0)
        yvals.append(medval)

    # Fit and apply radial correction
    x0 = [0.1,1.0]
    xfit, flag = scipy.optimize.curve_fit(radial_model,xvals,yvals)
    ans = xfit[0]*np.array(xvals) + xfit[1]
    print('Fitted model is',xfit)
    plt.plot(xvals,radial_model(np.array(xvals),xfit[0],xfit[1]),'gray',linestyle='-',linewidth=3)
    plt.ylim(ymin=0.8,ymax=2.0)
    plt.xlabel('Distance from pointing centre (deg)')
    plt.ylabel('Total/Peak flux')
    plt.savefig('radial-correction.png')
    plt.close()
    plt.cla()

    opencat = Table.read(nearcat)
    separation = sepn(opencat['RA']*deg2rad,opencat['DEC']*deg2rad,racent*deg2rad,deccent*deg2rad)*rad2deg
    for i in range(0,len(opencat)):
        #corval = model(separation[i],xfit[0],xfit[1])/np.min(model(np.array(xvals),xfit[0],xfit[1])) # To normalise to the value at min X
        corval = radial_model(separation[i],xfit[0],xfit[1]) # To normalise at 1.0
        newSI = opencat[i]['Peak_flux']*corval
        opencat[i]['Peak_flux'] = newSI
    radcorcat =  incatname.replace('.fits','_radcat.fits')
    opencat.write(radcorcat,overwrite=True)

    # Check corrections applied ok
    if debug == True:
        opencat = Table.read(radcorcat)
        opencat = filter_sources(opencat,2*beammaj*deg2arcsec)
        separation = sepn(opencat['RA']*deg2rad,opencat['DEC']*deg2rad,racent*deg2rad,deccent*deg2rad)*rad2deg
        xvals = []
        yvals = []
        for i in range(0,len(radbins)-1):
            binmin = radbins[i]
            binmax = radbins[i+1]
            print('Bin: %s - %s'%(binmin,binmax))
            sourcemin = np.where(separation > binmin)
            sourcemax = np.where(separation < binmax)
            meetscriterea = np.intersect1d(sourcemin,sourcemax)  
            medval = np.median(opencat['Total_flux'][meetscriterea]/opencat['Peak_flux'][meetscriterea])
            plt.plot((binmin+binmax)/2.0,medval,'g+',markersize=10)
            plt.plot(separation[meetscriterea],opencat['Total_flux'][meetscriterea]/opencat['Peak_flux'][meetscriterea],'b.',alpha=0.1)
            xvals.append((binmin+binmax)/2.0)
            yvals.append(medval)
        x0 = [0.1,1.0]
        xfit, flag = scipy.optimize.curve_fit(radial_model,xvals,yvals)
        ans = xfit[0]*np.array(xvals) + xfit[1]
        print('Fitted model is',xfit)
        plt.plot(xvals,radial_model(np.array(xvals),xfit[0],xfit[1]),'gray',linestyle='-',linewidth=3)
        plt.ylim(ymin=0.8,ymax=2.0)
        plt.xlabel('Distance from pointing centre (deg)')
        plt.ylabel('Total/Peak flux')
        plt.savefig('radial-correction-corrected.png')
        plt.close()
        plt.cla()
    return radcorcat

def find_only_compact(incatname,inimagename):
    
    envelopevals_point = []

    # Find some things from image
    f = fits.open(inimagename)
    racent = f[0].header['CRVAL1']
    deccent = f[0].header['CRVAL2']
    beammaj = f[0].header['BMAJ']
    
    
    opencat = Table.read(incatname)
    full_snr = opencat['Total_flux']/opencat['E_Total_flux']
    full_fratio = opencat['Total_flux']/opencat['Peak_flux']
    opencat = filter_sources(opencat,2*beammaj*deg2arcsec)
    filt_snr = opencat['Total_flux']/opencat['E_Total_flux']
    filt_fratio = opencat['Total_flux']/opencat['Peak_flux']


    steps = 6
    bins = histedges_equalN(filt_snr,steps)
    print('Bins',bins)

    indexes_real = []
    xvals = []
    yvals = []

    plthists = True

    for i in range(0,len(bins)-1):
        binmin = bins[i]
        binmax = bins[i+1]
        binmean = (binmax+binmin)/2.0
        print('Bin: %s - %s'%(binmin,binmax))

        sourcemin = np.where(filt_snr > binmin)
        sourcemax = np.where(filt_snr < binmax)
        meetscriterea_real = np.intersect1d(sourcemin,sourcemax)  

        # Change binvals dependning on the number of sources
        # Generally most things fall between -0.5 and 1.5. Ideally want like an average of 30 sources a bit or something.
        # Find binsize that gives 30 as a peak
        possiblebinsize = np.arange(0.05,0.3,0.01)
        for binsize in possiblebinsize:
            binvals = np.arange(-10,10.0,binsize)
            print(binvals)
            numberbin = np.max(np.histogram(np.array(filt_fratio)[meetscriterea_real],bins=binvals)[0])
            print(binvals,numberbin)
            if numberbin > 30:
                break
        print('Number of sources in compact cat in flux bin',len(meetscriterea_real))
        print('Max number at different Total/Int',numberbin)

        plt.hist(np.log(np.array(filt_fratio)[meetscriterea_real]),bins=binvals,alpha=0.5,density=True,histtype='step',color='g')

        distribution = find_bestfit_scipyfunction(np.log(np.array(filt_fratio)[meetscriterea_real]),binvals)
        distribution = getattr(stats,distribution)
        fitdistrib = distribution.fit(np.log(np.array(filt_fratio)[meetscriterea_real]))
        plt.plot(binvals,distribution.pdf(binvals,*fitdistrib),alpha=0.2,color='g')

        pointSFs = []
        for j in range(0,len(binvals)-1):
            pointSFs.append(distribution.sf(binvals[j],*fitdistrib))

        # Also plot the full LoTSS cat which contains extended sources.
        sourcemin = np.where(full_snr > binmin)
        sourcemax = np.where(full_snr < binmax)
        meetscriterea_real = np.intersect1d(sourcemin,sourcemax)  

        print('Number of sources in full cat in bin',len(meetscriterea_real))
        plt.hist(np.log(np.array(full_fratio)[meetscriterea_real]),bins=binvals,alpha=0.5,density=True,histtype='step',color='r')
        distribution = find_bestfit_scipyfunction(np.log(np.array(full_fratio)[meetscriterea_real]),binvals)
        distribution = getattr(stats,distribution)
        fitdistrib = distribution.fit(np.log(np.array(full_fratio)[meetscriterea_real]))
        plt.plot(binvals,distribution.pdf(binvals,*fitdistrib),alpha=0.2,color='r')

        # Calculate survival functions
        realSFs = []
        for j in range(0,len(binvals)-1):
            realSFs.append(distribution.sf(binvals[j],*fitdistrib))

        for j in range(0,len(binvals)-1):
            if realSFs[j] > 5*pointSFs[j]:
                plt.plot([binvals[j],binvals[j]],[0,10],'g-') # plotting position at which the real ones are 5 times more likely
                envelopevals_point.append(binvals[j])
                break
        print('All-point 5 times val found at %s corresponds to SF of %s and CDF of %s'%(binvals[j],distribution.sf(binvals[j],*fitdistrib),distribution.cdf(binvals[j],*fitdistrib)))

        fontsize=20
        plt.yticks(fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.xlim(xmin=-0.5,xmax=1.25)
        plt.ylim(ymin=0.0,ymax=6.0)
        plt.ylabel('Probability density',fontsize=fontsize)
        plt.xlabel(r'$\ln \left( \frac{S_I}{S_P} \right)$ ',fontsize=fontsize)
        plt.tight_layout()
        plt.savefig('realandsimsources-bin-%s-%s.png'%(np.round(binmin,2),np.round(binmax,2)))
        plt.cla()
        plt.close()

    plt.close()
    plt.cla()

    binmeans = (bins + np.roll(bins, -1))[:-1] / 2.0
    x0 = [0.7,2.0,105,0.4]
    xfit, flag = scipy.optimize.leastsq(envelope_residuals, x0, args=(envelopevals_point,binmeans))
    resid = envelope_residuals(xfit,binmeans,envelopevals_point)
    chisq = sum(resid*resid)
    xvals = np.arange(5,5000,1)
    print('Fitted model is',xfit)
    plt.plot(xvals,envelope_model(xvals,xfit),'gray',linestyle='-',linewidth=3)
    
    for i in range(0,len(bins)-1):
        binmin = bins[i]
        binmax = bins[i+1]
        binmean = (binmax+binmin)/2.0
        plt.plot(binmean,envelopevals_point[i],'r+',markersize=15)

    compact = 0
    extended = 0
    opencat = Table.read(incatname)
    extended_sources = []
    for i in range(0,len(full_snr)):
        Rval = envelope_model(full_snr[i],xfit)
        if np.log(full_fratio[i]) > Rval:
            extended +=1
            extended_sources.append(i)
            plt.plot(np.array(full_snr)[i],np.log(np.array(full_fratio))[i],'r.',markeredgecolor='r',markerfacecolor='r',markersize=2,alpha=0.1)
        else:
            compact +=1
            plt.plot(np.array(full_snr)[i],np.log(np.array(full_fratio))[i],'g.',markeredgecolor='g',markerfacecolor='g',markersize=2,alpha=0.1)

    print('%s extended and %s compact'%(extended,compact))

    opencat.remove_rows(extended_sources)
    opencat.write(incatname.replace('.fits','_compact.fits'),overwrite=True)
    fontsize=20
    plt.yticks(fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.xlim(xmin=5,xmax=300.0)
    plt.ylim(ymin=-1.0,ymax=3.0)
    plt.semilogx()
    plt.xlabel(r'$\frac{S_I}{\sigma_I}$',fontsize=fontsize)
    plt.ylabel(r'$\ln \left( \frac{S_I}{S_P} \right)$ ',fontsize=fontsize)
    plt.tight_layout()
    plt.savefig('fitted-skew-envelope.png')
    plt.close()
    plt.cla()

    print('The fitted envelope is ',xfit)


#inimagename = '/disks/paradata/shimwell/Beyond-DR2/archive_images/testing-astrometry/P014+08/image_full_ampphase_di_m.NS.app.restored.fits'
#incatname = 'image_full_ampphase_di_m.NS.offset_cat.fits'

#radcorcat = radial_correction(incatname,inimagename)
#compactcat = find_only_compact(radcorcat,inimagename)
