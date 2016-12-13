import pylab
import numpy as np
import pyfits
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy.random as npr

#Define various angle conversion factors
arcsec2deg=1.0/3600
arcmin2deg=1.0/60
deg2rad=np.pi/180
deg2arcsec = 1.0/arcsec2deg
rad2deg=180.0/np.pi
arcmin2rad=arcmin2deg*deg2rad
arcsec2rad=arcsec2deg*deg2rad
rad2arcmin=1.0/arcmin2rad
rad2arcsec=1.0/arcsec2rad
steradians2degsquared = (180.0/np.pi)**2.0
degsquared2steradians = 1.0/steradians2degsquared

def sepn(r1,d1,r2,d2):
    """
    Calculate the separation between 2 sources, RA and Dec must be
    given in radians. Returns the separation in radians
    """
    # NB slalib sla_dsep does this
    # www.starlink.rl.ac.uk/star/docs/sun67.htx/node72.html
    cos_sepn=np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
    sepn = np.arccos(cos_sepn)
    return sepn

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
    scat = pyfits.open(catalog)

    fitsimage = pyfits.open(fitsimage)
    fieldra = fitsimage[0].header['CRVAL1']
    fielddec = fitsimage[0].header['CRVAL2']
    fitsimage.close()
    radialseps = sepn(scat[1].data['RA_1']*deg2rad,scat[1].data['DEC_1']*deg2rad,fieldra*deg2rad,fielddec*deg2rad)*rad2deg

    pylab.plot(np.sort(scat[1].data['Total_flux_1']),np.sort(scat[1].data['Total_flux_1'])*np.median(np.array(scat[1].data['Total_flux_2']/options['%s_fluxfactor'%auxcatname])/np.array(scat[1].data['Total_flux_1'])),'r--')
    pylab.plot(scat[1].data['Total_flux_1'],scat[1].data['Total_flux_2']/options['%s_fluxfactor'%auxcatname],'ko')
    pylab.xlabel('Integrated LOFAR flux (Jy)')
    pylab.ylabel('Integrated TGSS flux (Jy)')
    pylab.xlim(xmin=0.01,xmax=0.5)
    pylab.ylim(ymin=0.01,ymax=0.5)
    pylab.semilogx()
    pylab.semilogy()
    equality = np.arange(0,10,0.01)
    pylab.plot(equality,equality,'k-')
    pylab.savefig(outname.replace('.png','_integrated.png'))
    pylab.close('all')
    pylab.cla()
    plt.clf()
    
    pylab.plot(np.sort(scat[1].data['Peak_flux_1']),np.sort(scat[1].data['Peak_flux_1'])*np.median(np.array(scat[1].data['Peak_flux_2']/options['%s_fluxfactor'%auxcatname])/np.array(scat[1].data['Peak_flux_1'])),'r--')
    pylab.plot(scat[1].data['Peak_flux_1'],scat[1].data['Peak_flux_2']/options['%s_fluxfactor'%auxcatname],'ko')
    pylab.xlabel('Peak LOFAR flux (Jy)')
    pylab.ylabel('Peak TGSS flux (Jy)')
    pylab.xlim(xmin=0.01,xmax=0.5)
    pylab.ylim(ymin=0.01,ymax=0.5)
    pylab.semilogx()
    pylab.semilogy()
    equality = np.arange(0,10,0.01)
    pylab.plot(equality,equality,'k-')
    pylab.savefig(outname.replace('.png','_peak.png'))
    pylab.close('all')
    pylab.cla()
    plt.clf()
    
    pylab.plot(radialseps,scat[1].data['Total_flux_1']/(scat[1].data['Total_flux_2']/options['%s_fluxfactor'%auxcatname]),'bo',alpha=0.4,markersize=3)
    equality = np.arange(-0.01,10,0.01)
    fractionrange = np.arange(0.99,1.01,0.002)
    for i in fractionrange:
        pylab.plot(equality,(equality/equality)*i,'k-')
    pylab.plot(equality,equality/equality,'k-')
    pylab.xlabel('Distance from pointing centre (deg)')
    pylab.ylabel('Integrated LOFAR flux / Integrated TGSS flux')
    
    pylab.xlim(xmin=0.0,xmax=2)
    pylab.ylim(ymin=0.5,ymax=2.0)    

    distancerange = np.arange(0,np.max(radialseps)+0.15,0.15)

    for i in range(0,len(distancerange)-1):
        distancemin = distancerange[i]
        distancemax = distancerange[i+1]
        binvals = np.array([])
        for j in range(0,len(radialseps)):
            if distancemin < radialseps[j] < distancemax:
                binvals = np.append(binvals,scat[1].data['Total_flux_1'][j]/(scat[1].data['Total_flux_2'][j]/options['%s_fluxfactor'%auxcatname]))
        midpoint = (distancemin+distancemax)/2.0
        if len(binvals) > 0.0:
            booterrl,booterrh = bootstrap(binvals, 100000, np.mean, 0.05)
            booterrl,booterrh = bootstrap(binvals, 100000, np.median, 0.05)
            print booterrl,booterrh, binvals
            pylab.errorbar(midpoint,np.median(binvals),xerr=(distancemin-distancemax)/2.0,yerr=[[np.median(binvals)-booterrl],[booterrh-np.median(binvals)]],fmt='--o',ecolor='b',color='b',zorder=999999)
    pylab.savefig(outname.replace('.png','_total_radial.png'))
    pylab.close('all')
    pylab.cla()
    plt.clf()


def plot_flux_ratios(catalog,fitsimage,outname,options=None):
    if options is None:
        options = o
    scat = pyfits.open(catalog)
    fluxratios = scat[1].data['Total_flux']/scat[1].data['Peak_flux']

    fitsimage = pyfits.open(fitsimage)
    fieldra = fitsimage[0].header['CRVAL1']
    fielddec = fitsimage[0].header['CRVAL2']
    fitsimage.close()
    radialseps = sepn(scat[1].data['RA_1']*deg2rad,scat[1].data['DEC_1']*deg2rad,fieldra*deg2rad,fielddec*deg2rad)*rad2deg

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
        #bsmear = bandwidth_smearing2(6*arcsec2deg,150E6,midpoint,4*48.8E3)
        #tsmear = time_smearing2(16.0,midpoint,6*arcsec2deg)
        #smids.append(midpoint)
        #svals.append(1.0/(bsmear*tsmear))
        booterrl,booterrh = bootstrap(binvals, 100000, np.mean, 0.05)
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
    scat = pyfits.open(catalog)

    fitsimage = pyfits.open(fitsimage)
    fieldra = fitsimage[0].header['CRVAL1']
    fielddec = fitsimage[0].header['CRVAL2']
    fitsimage.close()

    x = (scat[1].data['RA_1']-scat[1].data['RA_2'])*np.cos(fielddec*deg2rad)*deg2arcsec
    y = (scat[1].data['DEC_1']-scat[1].data['DEC_2'])*deg2arcsec

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
    plt.xlabel('$RA_{LOFAR}$ - $RA_{FIRST}$ (arcsec)')
    plt.ylabel('$DEC_{LOFAR}$ - $DEC_{FIRST}$ (arcsec)')
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
    lim = (int(xymax/binwidth) + 1) * binwidth

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
