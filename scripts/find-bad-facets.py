from auxcodes import *
import scipy.ndimage as nd
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.stats._continuous_distns import_distn_names
from astropy.table import Table

## ANOTHER SCRIPT TO FIND ENTIRE FILEDS THAT ARE BAD.

# Define a gaussian function with offset
def gaussian_func(x, a, x0, sigma,c):
    # Fix x0 to be 0.0
    #x0 = 0.0
    #c = 0.0
    return a * np.exp(-(x-x0)**2/(2*sigma**2)) + c

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.

    med = np.median(arr)
    return np.median(np.abs(arr - med))


def filterarray_outliers(dataarray,niter=20,eps=1e-6):
    rms=1
    oldrms=1
    for i in range(niter):
        if i == 1.0:
            rms=1E-3
        rms=np.std(dataarray)
        if old_div(np.abs(oldrms-rms),rms) < eps:
            return dataarray
        if len(dataarray) < 10000:
            return dataarray # Just return something if everything getting removed
        print('rms vals',len(dataarray),rms,oldrms)
        dataarray=dataarray[np.abs(dataarray)<5*rms]
        oldrms=rms
    print(dataarray,old_div(np.abs(oldrms-rms),rms)*1E6,old_div(np.abs(oldrms-rms),rms) < eps,eps,rms,oldrms,i)
    return dataarray
    raise RuntimeError('Failed to converge') # seems unreachable


# Make facet images
ds9region = 'image_full_ampphase_di_m.NS.tessel.reg'
infilename = 'image_full_low_m.app.restored.fits' # Note that in some tests the low res image was both faster to work with and easier to identify bad facets (might not always be true)
#infilename = 'image_full_ampphase_di_m.NS.app.restored.fits'
finaloutfilename = 'image_full_low_m_blanked.app.restored.fits' # Also could do high res image here
maskfile = 'image_full_low_m.mask01.fits'
cleanmask = fits.open(maskfile)
cleanmask[0].data[np.where(cleanmask[0].data==1)] = 2
cleanmask[0].data[np.where(cleanmask[0].data==0)] = 1
cleanmask[0].data[np.where(cleanmask[0].data==2)] = 0

polylist = convert_regionfile_to_poly(ds9region)
allastrooffsets = Table.read('pslocal-facet_offsets.fits')
debug = False

facetdata = {}
fullfacetdata = {}
fullfacetextendeddata = {}

facetsizes = {}
height_diffs = []
position_diffs = []
std_diffs = []
offset_diffs = []
missingfacets = []
facetids = []

for direction,ds9region in enumerate(polylist):

    r = pyregion.parse(ds9region)

    # Mask so that hduflat.data only contains a particular facet.
    hdu=fits.open(infilename)
    hduflat = flatten(hdu) # depending on MakeMask version may be 2D or 4D
    template=fits.open(infilename) # will be 4D
    outfilename = 'facet_%s.fits'%direction
    manualmask = r.get_mask(hdu=hduflat)
    facetsizes[direction] = np.sum(manualmask.astype(int))*abs(hdu[0].header['CDELT1']**2.0)
    manualmask_clean = cleanmask[0].data[0,0,:,:].astype(int) + manualmask.astype(int) - 1
    manualmask_clean[np.where(manualmask_clean == -1)] = 0
    manualmask_clean = manualmask_clean.astype(bool)
    fullfacetdata[direction] = hduflat.data[manualmask_clean]

    # Extend the manual mask a bit and remove facet (i.e. get other facet info)
    largemanualmask = nd.gaussian_filter(manualmask.astype(float),sigma=6) # Using a float makes the island of non zero values larger
    largemanualmask[largemanualmask>0.0] = 1.0
    largemanualmask = largemanualmask.astype(bool)
    fullfacetextendeddata[direction] = hduflat.data[largemanualmask]
    largemanualmask = largemanualmask.astype(int)-manualmask.astype(int)
    largemanualmask = largemanualmask.astype(int)+cleanmask[0].data[0,0,:,:].astype(int)-1
    largemanualmask[np.where(largemanualmask == -1)] = 0
    largemanualmask = largemanualmask.astype(bool)
    
    # make a small manual mask
    smallmanualmask = nd.gaussian_filter(manualmask,sigma=6) # This makes it smaller 
    smallmanualmask = manualmask.astype(int)-smallmanualmask.astype(int)
    smallmanualmask = smallmanualmask.astype(int)+cleanmask[0].data[0,0,:,:].astype(int)-1
    smallmanualmask[np.where(smallmanualmask == -1)] = 0
    smallmanualmask = smallmanualmask.astype(bool)

    
    # Write out shifted data
    if debug == True:
        outfile = fits.open(infilename)
        outfile[0].data[0,0][~smallmanualmask.astype(bool)] = np.nan
        outfile.writeto('facet_%s_small.fits'%direction,overwrite=True)
        outfile = fits.open(infilename)
        outfile[0].data[0,0][~largemanualmask.astype(bool)] = np.nan
        outfile.writeto('facet_%s_large.fits'%direction,overwrite=True)
        outfile = fits.open(infilename)
        outfile[0].data[0,0][~manualmask_clean.astype(bool)] = np.nan
        outfile.writeto('facet_%s_normal.fits'%direction,overwrite=True)
        
    # Filter data inside and outside facet
    print('Length before filter',np.shape(template[0].data[0,0][smallmanualmask.astype(bool)]))
    dataarrayin = filterarray_outliers(template[0].data[0,0][smallmanualmask.astype(bool)],niter=20,eps=1e-6)
    dataarrayout = filterarray_outliers(template[0].data[0,0][largemanualmask.astype(bool)],niter=20,eps=1e-6)
    
    #ksresult = ks_2samp(template[0].data[0,0][smallmanualmask.astype(bool)],template[0].data[0,0][largemanualmask.astype(bool)]) # This just gives pretty much 0 all the time...

    # Plot the data for inside and outside the facet
    print(len(dataarrayin),'length')
    if len(dataarrayin) == 0:
        missingfacets.append(direction)
        continue
    else:
        facetids.append(direction)
          
          
    bins = np.arange(-5*np.std(dataarrayin),5*np.std(dataarrayin),np.std(dataarrayin)/10.0)
    plt.hist(template[0].data[0,0][smallmanualmask.astype(bool)],bins=bins,histtype='step',color='g',density=True)
    plt.hist(template[0].data[0,0][largemanualmask.astype(bool)],bins=bins,histtype='step',color='r',density=True)

    # Fit pixels inside facet with gaussian
    yvals,xvals = np.histogram(dataarrayin,bins=bins,density=True)
    xcenters = (xvals[:-1]+xvals[1:])/2
    initialguess = [np.max(yvals), np.mean(xvals), np.std(xvals),0.0]
    xpopt1,pcov = curve_fit(gaussian_func,xcenters,yvals,p0=initialguess)
    plottingvals = np.arange(np.min(bins),np.max(bins),(np.max(bins)-np.min(bins))/1000.0)
    plt.plot(plottingvals,gaussian_func(np.array(plottingvals),xpopt1[0],xpopt1[1],xpopt1[2],xpopt1[3]),'g--')
    height_inner,position_inner,std_inner,offset_inner = xpopt1

    
    # Fit pixels outside facet with gaussian
    yvals,xvals = np.histogram(dataarrayout,bins=bins,density=True)
    xcenters = (xvals[:-1]+xvals[1:])/2
    initialguess = [np.max(yvals), np.mean(xvals), np.std(xvals),0.0]
    xpopt1,pcov = curve_fit(gaussian_func,xcenters,yvals,p0=initialguess)
    plottingvals = np.arange(np.min(bins),np.max(bins),(np.max(bins)-np.min(bins))/1000.0)
    plt.plot(plottingvals,gaussian_func(np.array(plottingvals),xpopt1[0],xpopt1[1],xpopt1[2],xpopt1[3]),'r--')
    height_outer,position_outer,std_outer,offset_outer = xpopt1

    # Save fig and data for inside facet boarder
    plt.savefig('facet_%s_hist.png'%(direction))
    plt.close()
    plt.cla()
    facetdata[direction] = dataarrayin
    print('Facet',direction,'done')

    # Detect if there was a jump
    height_diff = height_outer - height_inner # High positive bad (i.e. lower gaussian)
    position_diff = position_outer - position_inner # High negative bad (i.e. wider gaussian)
    std_diff = std_outer - std_inner
    offset_diff = offset_outer - offset_inner

    height_diffs.append(height_diff)
    position_diffs.append(position_diff)
    std_diffs.append(std_diff)
    offset_diffs.append(offset_diff)


    print('Height diff %s, position diff %s, std diff %s, offset diff %s'%(height_diff,position_diff,std_diff,offset_diff))
# Search for jumps between the facets that are significantly larger than normal jumps between facets

jump_height = np.where(height_diffs > (np.median(height_diffs) + 3*mad(height_diffs)))[0]
jump_std = np.where(std_diffs < (np.median(std_diffs) - 3*mad(std_diffs)))[0]

if len(jump_height) == 0 and len(jump_std) == 0:
    print('No jumps between facets found')
    
# Search for facets that are outliers.
stds = []
positions = []
heights = []
offsets = []
astrooffsets = []
astrooffset_errs = []
dynamicranges = []

#for direction in facetdata:
for direction,ds9region in enumerate(polylist):

    if direction in missingfacets:
        continue
    # Either choose just the outside bit (facetdata[direction]) or just use the full data from the facet (fullfacetdata[direction])
    dataarrayin = facetdata[direction]
    dataarrayin = filterarray_outliers(fullfacetdata[direction],niter=20,eps=1e-6)

    # Fit with gaussian
    yvals,xvals = np.histogram(dataarrayin,bins=bins,density=True)
    xcenters = (xvals[:-1]+xvals[1:])/2
    initialguess = [np.max(yvals), np.mean(xvals), np.std(xvals),0.0]

    xpopt1,pcov = curve_fit(gaussian_func,xcenters,yvals,p0=initialguess)
    plottingvals = np.arange(np.min(bins),np.max(bins),(np.max(bins)-np.min(bins))/1000.0)
    plt.plot(plottingvals,gaussian_func(np.array(plottingvals),xpopt1[0],xpopt1[1],xpopt1[2],xpopt1[3]),'g--')
    plt.hist(dataarrayin,bins=bins,density=True,histtype='step',color='g')

    height,position,std,offset = xpopt1

    # Change the dynamic range so it looks in an extended region around the facet.
    dynamicrange= abs(np.max(fullfacetextendeddata[direction])/np.min(fullfacetextendeddata[direction]))


    stds.append(std)
    positions.append(position)
    heights.append(height)
    offsets.append(offset)
    dynamicranges.append(dynamicrange)
    astrooffset = np.sqrt(allastrooffsets['RA_offset'][direction]**2.0 + allastrooffsets['DEC_offset'][direction]**2.0)
    astrooffset_err = np.sqrt(allastrooffsets['RA_error'][direction]**2.0 + allastrooffsets['DEC_error'][direction]**2.0)
    astrooffsets.append(astrooffset)
    astrooffset_errs.append(astrooffset_err)
    print('Facet %s: std %s, postion %s, hieght %s, noise offset %s, astro offset %s dynamci range %s'%(direction,std,position,height,offset,astrooffset,dynamicrange))
    #print('astro offset',np.sqrt(allastrooffsets[direction][2]**2.0 + allastrooffsets[direction][3]**2.0))

print('MAD std %s, position %s, height %s, offset %s, astro %s dynamc range %s'%(mad(stds),mad(positions),mad(heights),mad(offsets),mad(astrooffsets),mad(dynamicranges)))
badstd = np.where((abs(stds-np.median(stds)))>3*mad(stds))[0]
badposition= np.where((abs(positions-np.median(positions)))>3*mad(positions))[0]
badheight = np.where((abs(heights-np.median(heights)))>3*mad(heights))[0]
badoffset = np.where((abs(offsets-np.median(offsets)))>3*mad(offsets))[0]
badastro =  np.where((abs(astrooffsets-np.median(astrooffsets)))>3*mad(astrooffsets))[0]
badastro_err =  np.where((abs(astrooffset_errs-np.median(astrooffset_errs)))>3*mad(astrooffset_errs))[0]
baddyncut = np.max([np.median(dynamicranges)- 3*mad(dynamicranges),30])
baddynamic =  np.where(dynamicranges < baddyncut)[0] # np.where(np.array(dynamicranges) < 30)[0]

gooddynamic = np.where(dynamicranges > (np.median(dynamicranges)+mad(dynamicranges)))[0] #np.where(np.array(dynamicranges) > 200)[0]
                        
allbad = []
smallestfacets = []
recordsizes = []
for i in facetsizes:
    recordsizes.append(facetsizes[i])
smallfacets = np.sort(recordsizes)[:3]
for i in range(0,len(recordsizes)):
    if facetsizes[i] in smallfacets:
        smallestfacets.append(i)
print('bad noise std',badstd)
print('bad noise position',badposition)
print('bad noise height',badheight)
print('bad noise offset',badoffset)
print('bad astrometry',badastro)
print('bad astrometry err',badastro_err)
print('bad dynamic range',baddynamic)
print('Good dynamic range',gooddynamic)
print('Smallest facets',smallestfacets)
print('Facets with std jumps on boundaries',jump_std)
print('Facets with height jumps on boundaries',jump_height)

for i in range(0,len(polylist)):
    counterbad = 0
    if i in badstd:
        counterbad +=1
    if i in badoffset:
        counterbad +=1
    if i in badheight:
        counterbad +=1
    if i in badposition:
        counterbad +=1
    if i in badastro:
        counterbad +=1
    if i in badastro_err:
        counterbad +=1
    if i in baddynamic:
        counterbad +=1
    if i in jump_std and i in badstd:
        counterbad +=2
    if i in jump_height and i in badheight:
        counterbad+=2
    if i in gooddynamic: # Always keep ones with a good dynamic range
        counterbad = 0
    if i in smallestfacets and i not in baddynamic:
        counterbad = 0 # Smallest ones often around bright sources.
    print('Facet %s: %s points'%(i,counterbad))
    if counterbad > 5:
        print('Facet %s: %s points ***BAD***'%(i,counterbad))
        allbad.append(i)
print('Final Bad',allbad)

plt.savefig('all_facet_hists.png')
plt.close()
plt.cla()

hdu = fits.open(infilename)
hdu2 = fits.open('image_full_low_m.int.restored.fits')
fpb = hdu[0].data[0,0]/hdu2[0].data[0,0]
fpb[np.where(fpb < 0.3)] = 0.0
within03 = []
hdu[0].data[0,0] = fpb
#f.writeto('pbcutfile.fits',overwrite='True')
#hdu=fits.open(infilename)
#hdu=fits.open('pbcutfile.fits')
hduflat = flatten(hdu) # depending on MakeMask version may be 2D or 4D
for direction,ds9region in enumerate(polylist):
    #print(direction,ds9region)
    r = pyregion.parse(ds9region)
    if direction in missingfacets:
        continue
    # Mask so that hduflat.data only contains a particular facet.
    manualmask = r.get_mask(hdu=hduflat)
    try:
        fraczero = 1.0-len(np.where(hduflat.data[manualmask]==0.0)[0])/len(hduflat.data[manualmask])
    except ZeroDivisionError:
        #continue
        fraczero = 0.0 # possible for facet to be outside low image. If so s
    print(fraczero,direction)
    within03.append(fraczero)

hdu=fits.open(infilename)
hduflat = flatten(hdu) # depending on MakeMask version may be 2D or 4D
for direction,ds9region in enumerate(polylist):
    print(direction,ds9region)
    r = pyregion.parse(ds9region)
    # Mask so that hduflat.data only contains a particular facet.
    manualmask = r.get_mask(hdu=hduflat)
    if direction in allbad:
        hduflat.data[manualmask] = np.nan
hdu[0].data[0,0]= hduflat.data
hdu.writeto(finaloutfilename,overwrite=True)


saveT = Table()
saveT = Table(names=('Facet_id','Facet_GausFit_std','Facet_GausFit_xoffset','Facet_GausFit_yoffset','Facet_GausFit_amp','Facet_astro_offset','Facet_astro_err','Facet_dynrange','Facet_size','Boundary_jump_std','Boundary_height_std','fracwithin03'),dtype=('i4','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
for direction in range(0,len(facetids)):#,ds9region in enumerate(polylist):
    facetid = facetids[direction]
    facetsize = facetsizes[facetid]
    saveT.add_row((facetid,stds[direction],positions[direction],offsets[direction],heights[direction],astrooffsets[direction],astrooffset_errs[direction],dynamicranges[direction],facetsize,std_diffs[direction],height_diffs[direction],within03[direction]))
saveT.write('Facet_statistics.fits',overwrite=True)
