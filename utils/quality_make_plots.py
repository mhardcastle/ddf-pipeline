from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from crossmatch_utils import separation
from matplotlib.ticker import NullFormatter
import numpy.random as npr
from astropy.table import Table
from astropy.io import fits
from auxcodes import sepn,warn
import os

#Define various angle conversion factors
arcsec2deg=1.0/3600
arcmin2deg=1.0/60
deg2rad=old_div(np.pi,180)
deg2arcsec = 1.0/arcsec2deg
rad2deg=180.0/np.pi
arcmin2rad=arcmin2deg*deg2rad
arcsec2rad=arcsec2deg*deg2rad
rad2arcmin=1.0/arcmin2rad
rad2arcsec=1.0/arcsec2rad
steradians2degsquared = (180.0/np.pi)**2.0
degsquared2steradians = 1.0/steradians2degsquared

def bootstrap(data, num_samples, statistic, alpha):
    """Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic."""
    n = len(data)
    idx = npr.randint(0, n, (num_samples, n))
    samples = data[idx]
    stat = np.sort(statistic(samples, 1))
    return (stat[int((alpha/2.0)*num_samples)],
            stat[int((1-alpha/2.0)*num_samples)])

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))


def plot_flux_errors(catalog,fitsimage,outname,auxcatname,options=None):
    if options is None:
        options = o
    if os.path.isfile(outname.replace('.png','_integrated.png')):
        warn('Plot file %s exists, not making it' % outname)
    else:

        scat = Table.read(catalog)

        fitsimage = fits.open(fitsimage)
        fieldra = fitsimage[0].header['CRVAL1']
        fielddec = fitsimage[0].header['CRVAL2']
        fitsimage.close()
        radialseps = sepn(scat['RA']*deg2rad,scat['DEC']*deg2rad,fieldra*deg2rad,fielddec*deg2rad)*rad2deg

        plt.plot(np.sort(scat['Total_flux']),np.sort(scat['Total_flux'])*np.median(old_div(np.array(old_div(scat[auxcatname+'_Total_flux'],options[auxcatname+'_fluxfactor'])),np.array(scat['Total_flux']))),'r--')
        plt.plot(scat['Total_flux'],old_div(scat[auxcatname+'_Total_flux'],options[auxcatname+'_fluxfactor']),'ko')
        plt.xlabel('Integrated LOFAR flux (Jy)')
        plt.ylabel('Integrated '+auxcatname+' flux (Jy)')
        plt.xlim(xmin=0.01,xmax=0.5)
        plt.ylim(ymin=0.01,ymax=0.5)
        plt.semilogx()
        plt.semilogy()
        equality = np.arange(0,10,0.01)
        plt.plot(equality,equality,'k-')
        plt.savefig(outname.replace('.png','_integrated.png'))
        plt.close('all')
        plt.cla()
        plt.clf()

        plt.plot(np.sort(scat['Peak_flux']),np.sort(scat['Peak_flux'])*np.median(old_div(np.array(old_div(scat[auxcatname+'_Peak_flux'],options['%s_fluxfactor'%auxcatname])),np.array(scat['Peak_flux']))),'r--')
        plt.plot(scat['Peak_flux'],old_div(scat[auxcatname+'_Peak_flux'],options[auxcatname+'_fluxfactor']),'ko')
        plt.xlabel('Peak LOFAR flux (Jy)')
        plt.ylabel('Peak '+auxcatname+' flux (Jy)')
        plt.xlim(xmin=0.01,xmax=0.5)
        plt.ylim(ymin=0.01,ymax=0.5)
        plt.semilogx()
        plt.semilogy()
        equality = np.arange(0,10,0.01)
        plt.plot(equality,equality,'k-')
        plt.savefig(outname.replace('.png','_peak.png'))
        plt.close('all')
        plt.cla()
        plt.clf()

        plt.plot(radialseps,old_div(scat['Total_flux'],(old_div(scat[auxcatname+'_Total_flux'],options[auxcatname+'_fluxfactor']))),'bo',alpha=0.4,markersize=3)
        equality = np.arange(-0.01,10,0.01)
        fractionrange = np.arange(0.99,1.01,0.002)
        for i in fractionrange:
            plt.plot(equality,1.0*i+0*equality,'k-')
        plt.plot(equality,1.0+0*equality,'k-')
        plt.xlabel('Distance from pointing centre (deg)')
        plt.ylabel('Integrated LOFAR flux / Integrated '+auxcatname+' flux')

        plt.xlim(xmin=0.0,xmax=2)
        plt.ylim(ymin=0.5,ymax=2.0)    

        distancerange = np.arange(0,np.max(radialseps)+0.15,0.15)

        for i in range(0,len(distancerange)-1):
            distancemin = distancerange[i]
            distancemax = distancerange[i+1]
            binvals = np.array([])
            for j in range(0,len(radialseps)):
                if distancemin < radialseps[j] < distancemax:
                    binvals = np.append(binvals,old_div(scat['Total_flux'][j],(old_div(scat[auxcatname+'_Total_flux'][j],options[auxcatname+'_fluxfactor']))))
            midpoint = (distancemin+distancemax)/2.0
            if len(binvals) > 0:
                booterrl,booterrh = bootstrap(binvals, 100000, np.mean, 0.05)
                booterrl,booterrh = bootstrap(binvals, 100000, np.median, 0.05)
    #            print booterrl,booterrh, binvals
                plt.errorbar(midpoint,np.median(binvals),xerr=(distancemin-distancemax)/2.0,yerr=[[np.median(binvals)-booterrl],[booterrh-np.median(binvals)]],fmt='--o',ecolor='b',color='b',zorder=999999)
        plt.savefig(outname.replace('.png','_total_radial.png'))
        plt.close('all')
        plt.cla()
        plt.clf()


def plot_flux_ratios(catalog,fitsimage,outname,options=None):
    if options is None:
        options = o

    if os.path.isfile(outname):
        warn('Plot file %s exists, not making it' % outname)
    else:
        scat = Table.read(catalog)
        fluxratios = old_div(scat['Total_flux'],scat['Peak_flux'])

        fitsimage = fits.open(fitsimage)
        fieldra = fitsimage[0].header['CRVAL1']
        fielddec = fitsimage[0].header['CRVAL2']
        fitsimage.close()
        radialseps = sepn(scat['RA']*deg2rad,scat['DEC']*deg2rad,fieldra*deg2rad,fielddec*deg2rad)*rad2deg

        nullfmt = NullFormatter()         # no labels

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.8
        bottom_h = left_h = left + width + 0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.2, height]

        # start with a rectangular Figure
        plt.figure(1, figsize=(8, 8))
        axScatter = plt.axes(rect_scatter)
        plt.xlabel('Distance from pointing centre (deg)')
        plt.ylabel('Integrated flux / Peak flux')
        plt.xticks(rotation=270)
        axScatter.plot([0,3.0],[1.0,1.0],'k-')

        smids = []
        svals = []
        distancerange = np.arange(0,max(radialseps),0.3)
        for i in range(0,len(distancerange)-1):
            distancemin = distancerange[i]
            distancemax = distancerange[i+1]
            binvals = np.array([])
            for j in range(0,len(radialseps)):
                if distancemin < radialseps[j] < distancemax:
                    binvals = np.append(binvals,fluxratios[j])
            midpoint = (distancemin+distancemax)/2.0
            #print i,binvals
            #bsmear = bandwidth_smearing2(6*arcsec2deg,150E6,midpoint,4*48.8E3)
            #tsmear = time_smearing2(16.0,midpoint,6*arcsec2deg)
            #smids.append(midpoint)
            #svals.append(1.0/(bsmear*tsmear))
            #booterrl,booterrh = bootstrap(binvals, 100000, np.mean, 0.05)
            if len(binvals)>0:
                booterrl,booterrh = bootstrap(binvals, 100000, np.median, 0.05)
                axScatter.errorbar(midpoint,np.median(binvals),xerr=(distancemin-distancemax)/2.0,yerr=[[np.median(binvals)-booterrl],[booterrh-np.median(binvals)]],fmt='--o',ecolor='b',color='b',zorder=999999)

        #axScatter.plot(smids,svals,'k--')

        axHisty = plt.axes(rect_histy)
        plt.xticks(rotation=270)

        # no labels
        axHisty.yaxis.set_major_formatter(nullfmt)

        # the scatter plot:

        axScatter.scatter(radialseps,fluxratios,color='blue',marker='+',s=10,alpha=0.5)

        binwidth = 0.05

        axScatter.set_xlim((0, 3.0))
        axScatter.set_ylim((0.8, 2.5))

        bins = np.arange(0.8, 3, binwidth)
        axHisty.hist(fluxratios, bins=bins, orientation='horizontal',color='blue',alpha=0.5)
        axHisty.set_ylim(axScatter.get_ylim())
        plt.savefig(outname)
        plt.close('all')
        plt.cla()
        plt.clf()

def plot_position_offset(catalog,fitsimage,outname,auxcatname,options=None):
    if options is None:
        options = o
    if os.path.isfile(outname):
        warn('Plot file %s exists, not making it' % outname)
    else:
        scat = Table.read(catalog)

        fitsimage = fits.open(fitsimage)
        fieldra = fitsimage[0].header['CRVAL1']
        fielddec = fitsimage[0].header['CRVAL2']
        fitsimage.close()

        x = scat[auxcatname+'_dRA']
        y = scat[auxcatname+'_dDEC']

        nullfmt = NullFormatter()         # no labels

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left + width + 0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.2, height]

        # start with a rectangular Figure
        plt.figure(1, figsize=(8, 8))

        axScatter = plt.axes(rect_scatter)
        plt.xlabel('$RA_{\\rm LOFAR} - RA_{\\rm %s}$ (arcsec)' % auxcatname)
        plt.ylabel('$DEC_{\\rm LOFAR} - DEC_{\\rm %s}$ (arcsec)' % auxcatname)
        plt.xticks(rotation=270)
        axScatter.plot(np.zeros(40),np.arange(-20,20),'k--')
        axScatter.plot(np.arange(-20,20),np.zeros(40),'k--')
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        plt.xticks(rotation=270)

        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)

        # the scatter plot:
        axScatter.scatter(x, y,marker='+',s=2,alpha=0.3)

        binwidth = 0.25
        xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
        lim = (int(old_div(xymax,binwidth)) + 1) * binwidth

        axScatter.set_xlim((-lim, lim))
        axScatter.set_ylim((-lim, lim))

        bins = np.arange(-lim, lim + binwidth, binwidth)
        axHistx.hist(x, bins=bins)
        axHisty.hist(y, bins=bins, orientation='horizontal')

        axHistx.set_xlim(axScatter.get_xlim())
        axHisty.set_ylim(axScatter.get_ylim())

        axScatter.add_artist(Ellipse((np.median(x),np.median(y)),np.std(x),np.std(y),angle=0,linewidth=2,fill=False,color='k',zorder=10000))

        plt.savefig(outname)
        plt.close('all')
        plt.cla()
        plt.clf()
