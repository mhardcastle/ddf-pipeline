import os,sys
from auxcodes import *
import scipy.ndimage as nd
import numpy as np


# Use scipy.ndimage.shift
# First get pixels for each facet with a bit of padding around the edge. Then apply the shift for that particular facet.
# Can probably do the facets in parallel but not done
# Then add together to make final image. If 2 pixels contribute then average.
# Make facet images

offsets = np.load('pslocal-facet_offsets.npy')
ds9region = 'image_full_ampphase_di_m.NS.tessel.reg'
infilename = 'image_full_ampphase_di_m.NS.app.restored.fits'
finaloutname = 'image_full_ampphase_di_m.NS_newshift.app.restored.fits'
polylist = convert_regionfile_to_poly(ds9region)
finalout = fits.open(infilename)
finalout[0].data[0,0] = finalout[0].data[0,0]*0.0
averagefilename = 'image_full_ampphase_di_m.NS_average.app.restored.fits'
averageout = fits.open(infilename)
averageout[0].data[0,0] = finalout[0].data[0,0]*0.0
debug = False

template=fits.open(infilename) # will be 4D

for direction,ds9region in enumerate(polylist):
    print(direction,ds9region)
    r = pyregion.parse(ds9region)

    # Mask so that hduflat.data only contains a particular facet.
    hdu=fits.open(infilename)
    hduflat = flatten(hdu) # depending on MakeMask version may be 2D or 4D
    template=fits.open(infilename) # will be 4D
    outfilename = 'facet_%s.fits'%direction
    manualmask = r.get_mask(hdu=hduflat)

    # Extende the manual mask a bit
    manualmask = nd.gaussian_filter(manualmask.astype(float),sigma=3)
    manualmask[manualmask>0.0] = 1.0
    manualmask = manualmask.astype(bool)
    np.where(manualmask==True)

    # Set data not in facet to 0
    hduflat.data[~manualmask] = 0.0 #np.random.normal(0.0,np.std(hduflat.data[~manualmask]),len(hduflat.data[~manualmask])) # Could have it at np.nan but that messes up interpolation. Or could have it to random noise. But 0 seems fine when also having the padding.

    if debug == True:
        # Write out the iamge before shifting.
        template[0].data[0,0]=hduflat.data
        template.writeto(outfilename,overwrite=True)
        tmp = fits.open(outfilename)
        tmp[0].data[0,0][~manualmask] = np.nan
        tmp.writeto(outfilename,overwrite=True) #

    # Do shifting
    outfilename = 'facet_%s_shifted.fits'%direction
    raoffset = offsets[direction][0]/1.5 # Still need to test what happens for fields around RA=0...
    decoffset = -offsets[direction][1]/1.5 # Tested and was indeed minus.

    raoffset_e = offsets[direction][2]/1.5
    decoffset_e = offsets[direction][3]/1.5

    move_amount = np.sqrt(raoffset**2.0 + decoffset**2.0)
    move_error = np.sqrt(raoffset_e**2.0 + decoffset_e**2.0)
    if move_amount < move_error:
        print('Not applying RA offset of %s pixels and DEC offset %s pixels to facet %s as errors are %s (RA) %s (DEC)'%(raoffset,decoffset,direction,raoffset_e,decoffset_e))
	continue
    print('Applying RA offset %s pixels and DEC offset %s pixels to facet %s'%(raoffset,decoffset,direction))
    shifteddata = nd.shift(hduflat.data,[decoffset,raoffset],order=3,cval=np.nan) # Shift in pixels

    # Shrink mask again to get rid of edge effets
    manualmask = nd.gaussian_filter(manualmask,sigma=1)

    # Write out image after shifting.
    if debug == True:
        template[0].data[0,0] = shifteddata
        template.writeto(outfilename,overwrite=True)
        tmp = fits.open(outfilename)
        tmp[0].data[0,0][~manualmask] = np.nan
        tmp.writeto(outfilename,overwrite=True)
    
 
    # Add shifted data to final image and create image so that we know where more than 1 facetimage goes into the final image.
    finalout[0].data[0,0][manualmask] += shifteddata[manualmask]
    averageout[0].data[0,0][manualmask] += 1.0


finalout[0].data[0,0] = finalout[0].data[0,0]/averageout[0].data[0,0]
finalout.writeto(finaloutname,overwrite=True)
averageout.writeto(averagefilename,overwrite=True)
