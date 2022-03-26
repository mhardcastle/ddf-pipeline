#!/usr/bin/env python

# Script to align an extracted field with the LoTSS-DR2 catalogue or an a catalogue produced from running the lotss quality-pipeline.py
# Inputs are:
# fitsimage -- this is the image from the extract pipeline (should be about 6" resolution)
# catalogue -- this is the catalogue of LoTSS-DR2 (if field in lotss-DR2 or this is the catalogue from the lotss quality-pipeline if not in LoTSS-DR2)
# regionfile -- this is the region file that comes with the extract tar.gz file (so outside this region sources are subtracted)
# fieldname -- ONLY USE IF NOT USING LOTSS DR2 - this is the name of the LoTSS field if using the LoTSS database (i.e. PXXX+XX) or any fieldname if not using that.
# nodatabase -- ONLY USE IF NOT USING LOTSS DR2 - this means that  you have to provide a scaling factor rather than take from the LoTSS database automatically
# fieldfactor -- ONLY USE IF NOT USING LOTSS DR2 - this is the scaling factor that you multiply you catalogue by (so from the quality-pipeline the nvss_factor/5.9124 - if a lotss field this can be found on https://lofar-surveys.org/widefields.html and is the "scale" number in the table)


# intended as a one-stop shop for mosaicing
# contains some of the same arguments as mosaic.py
from astropy.coordinates import SkyCoord
from astropy.table import Table
import argparse
import os
import glob
from astropy.io import fits
import numpy as np
import pickle
import os,sys
try:
    import bdsf as bdsm
except ImportError:
    import lofar.bdsm as bdsm
from crossmatch_utils import *
import pyregion
from astropy import wcs
from astropy.wcs import WCS
from random import random
from random import randint
from random import seed
from numpy import arange
from numpy import mean
from numpy import std
from numpy import absolute
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import HuberRegressor
from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import TheilSenRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from matplotlib import pyplot

def flatten(f):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    naxis=f[0].header['NAXIS']
    if naxis<2:
        print('Can\'t make map from this')
        sys.exit(0)
    if naxis==2:
        return fits.PrimaryHDU(header=f[0].header,data=f[0].data)

    w = WCS(f[0].header)
    wn=WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2
    copy=('EQUINOX','EPOCH','BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
    for k in copy:
        r=f[0].header.get(k)
        if r is not None:
            header[k]=r

    slice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            slice.append(np.s_[:],)
        else:
            slice.append(0)
        
    hdu = fits.PrimaryHDU(header=header,data=f[0].data[slice])
    return hdu


def filter_outside_extract(ds9region,infilename,catalogue):

    hdu=fits.open(infilename)
    hduflat = flatten(hdu)
    map=hdu[0].data
    w = WCS(flatten(hdu).header)

    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)
    inregion = []
    for element in catalogue:
        i,j = w.wcs_world2pix(element['RA'],element['DEC'],0)
        #print element['RA'], element['DEC'],'RA,DEC',i,j
        if i < 0.0 or j < 0.0:
            # outside map
            inregion.append(False)
            continue
        try:
            
            inregion.append(manualmask[int(j),int(i)])
        except IndexError:
            #print 'Going into exception'
            inregion.append(False)
    return catalogue[inregion]

def filter_inside_extract(ds9region,infilename,catalogue):

    hdu=fits.open(infilename)
    hduflat = flatten(hdu)
    map=hdu[0].data
    w = WCS(flatten(hdu).header)

    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)
    inregion = []
    for element in catalogue:
        i,j = w.wcs_world2pix(element['RA'],element['DEC'],0)
        if i < 0.0 or j < 0.0:
            # outside map
            inregion.append(False)
            continue
        try:
            inregion.append(~manualmask[int(j),int(i)])
        except IndexError:
            #print 'Going into exception'
            inregion.append(False)
    return catalogue[inregion]

def evaluate_model(X, y, model):
	# define model evaluation method
	cv = RepeatedKFold(n_splits=5, n_repeats=3, random_state=1)
	# evaluate model
	scores = cross_val_score(model, X, y, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1)
	# force scores to be positive
	return absolute(scores)

# dictionary of model names and model objects
def get_models():
	models = list()
	models.append(LinearRegression(fit_intercept=False))
	models.append(HuberRegressor(fit_intercept=False))
	#models.append(RANSACRegressor())#fit_intercept=False)) # Doesnt have option to not fit the intercept
	models.append(TheilSenRegressor(fit_intercept=False)) # Strunggling a bit with this one as the output varies a lot given n_samples (if n_samples=1 then it returns the median of the ratio, if it equals the number of data points then it returns essentially the output of least square fitting)
	return models
 
# plot the dataset and the model's line of best fit
def plot_best_fit(X, y, xaxis, model):
    # fit the model on all data
    model.fit(X, y)
    # calculate outputs for grid across the domain
    yaxis = model.predict(xaxis.reshape((len(xaxis), 1)))
    # plot the line of best fit
    pyplot.plot(xaxis, yaxis, label=type(model).__name__)
    results = evaluate_model(X, y, model)
    gradient = model.coef_
    intercept = model.intercept_
    modelname =  type(model).__name__
    print(modelname,'Gradient',gradient,'intercept',intercept,'Mean MAE: %.3f (%.3f)' % (mean(results), std(results)))
        
    return modelname,gradient,intercept,mean(results),std(results)

# Run bdsm on the image
parser = argparse.ArgumentParser(description='fitsimage')
parser.add_argument('fitsimage', type=str, help='fitsimage')
parser.add_argument('catalogue', type=str, help='The LoTSS-DR2 catalogue or catalogue from another field')
parser.add_argument('regionfile', type=str, help='extractionregion')
parser.add_argument('--fieldname', type=str,default='LoTSS-DR2',help='field name if not LoTSS-DR2 but using surveys database')
parser.add_argument('--nodatabase', help='Do not use LOFAR surveys fields database', action='store_true')
parser.add_argument('--fieldfactor', help='Scaling factor of the field if not using surveys database',type=float)
args = parser.parse_args()

infile = args.fitsimage
catalogue = args.catalogue
regionfile = args.regionfile
fieldname = args.fieldname
nodatabase = args.nodatabase
fieldfactor = args.fieldfactor

restfrq=143.65e6 # should work this out from the FITS headers eventually

f = fits.open(infile)
ref_ra = f[0].header['CRVAL1']
if ref_ra < 0.0:
	ref_ra  = 360+ref_ra
ref_dec = f[0].header['CRVAL2']

if  not os.path.exists(infile.replace('.fits','cat.srl.fits')):
    img = bdsm.process_image(infile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=True, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None, frequency=restfrq)
    img.write_catalog(outfile=infile.replace('.fits','cat.srl.fits'),catalog_type='srl',format='fits',correct_proj='True')
    img.write_catalog(outfile=infile.replace('.fits','cat.srl.reg'),catalog_type='srl',format='ds9',correct_proj='True')
# LoTSS-DR2 catalogue
# Filter DR2 cat

lotssdr2=Table.read(catalogue)
print('original length:',len(lotssdr2))
lotssdr2=filter_catalogue(lotssdr2,ref_ra,ref_dec,1.0)
print('filter around',ref_ra,ref_dec)
print('filter to 1.0 deg:',len(lotssdr2))
# Cut to match extraction region
lotssdr2 = filter_outside_extract(regionfile,infile,lotssdr2)
print(len(lotssdr2),'region filtered')
if fieldname == 'LoTSS-DR2':
    print(lotssdr2['Mosaic_ID'][0])
else:
    if nodatabase:
        print('Scaling factor being used for field catalogue',fieldfactor)
        lotssdr2['Peak_flux'] = lotssdr2['Peak_flux']*fieldfactor*1000.0 #1000 scaling is to match mJy which is the units of LoTSS-DR2 cat
        lotssdr2['Total_flux'] = lotssdr2['Total_flux']*fieldfactor*1000.0 #1000 scaling is to match mJy which is the units of LoTSS-DR2 cat
        lotssdr2['E_Peak_flux'] = lotssdr2['E_Peak_flux']*fieldfactor*1000.0 #1000 scaling is to match mJy which is the units of LoTSS-DR2 cat
        lotssdr2['E_Total_flux'] = lotssdr2['E_Total_flux']*fieldfactor*1000.0 #1000 scaling is to match mJy which is the units of LoTSS-DR2 cat
    else:
        from surveys_db import SurveysDB
        sdb = SurveysDB()
        qualitydict = sdb.get_quality(fieldname)
        sdb.close()
        nvssfactor= 1.0/(qualitydict['nvss_scale']/5.9124)
        tgssscale = qualitydict['tgss_scale']
        print('Scaling comparison catalogue by nvssfactor', nvssfactor)
        lotssdr2['Peak_flux'] = lotssdr2['Peak_flux']*nvssfactor*1000.0 #1000 scaling is to match mJy which is the units of LoTSS-DR2 cat
        lotssdr2['Total_flux'] = lotssdr2['Total_flux']*nvssfactor*1000.0 #1000 scaling is to match mJy which is the units of LoTSS-DR2 cat
        lotssdr2['E_Peak_flux'] = lotssdr2['E_Peak_flux']*nvssfactor*1000.0 #1000 scaling is to match mJy which is the units of LoTSS-DR2 cat
        lotssdr2['E_Total_flux'] = lotssdr2['E_Total_flux']*nvssfactor*1000.0 #1000 scaling is to match mJy which is the units of LoTSS-DR2 cat



#lotssdr2=select_isolated_sources(lotssdr2,45)
print('isolated sources',len(lotssdr2))
lotsssnr = lotssdr2['Peak_flux']/(2.0*lotssdr2['E_Peak_flux']) + lotssdr2['Total_flux']/(2.0*lotssdr2['E_Total_flux'])
lotssdr2=lotssdr2[lotsssnr > 7.0]
print('snr  more than 7 sources',len(lotssdr2))
lotsssnr = lotssdr2['Peak_flux']/(2.0*lotssdr2['E_Peak_flux']) + lotssdr2['Total_flux']/(2.0*lotssdr2['E_Total_flux'])
lotsscompact = 0.41 + (1.10/(1.0+(lotsssnr/104.32)**2.03))
lotssR = np.log(lotssdr2['Total_flux']/lotssdr2['Peak_flux'])
lotssdr2=lotssdr2[lotssR < lotsscompact]
print('Compact sources len',len(lotssdr2))
#print 'Compact sources',len(lotssdr2)

# Filter cutout cat
cutout=Table.read(infile.replace('.fits','cat.srl.fits'))
if cutout['RA'][0] < 0.0:
	cutout['RA'] = 360.0+cutout['RA']
print('original length:',len(cutout))
cutout=filter_catalogue(cutout,ref_ra,ref_dec,1.0)
print('filter to 1.0 deg:',len(cutout))
cutout=select_isolated_sources(cutout,30)
print('isolated sources',len(cutout))
cutoutsnr = cutout['Peak_flux']/(2.0*cutout['E_Peak_flux']) + cutout['Total_flux']/(2.0*cutout['E_Total_flux'])
cutout=cutout[cutoutsnr > 7.0]
print('snr  more than 7 sources',len(cutout))
#cutout=cutout[cutout['S_Code'] == 'S']
#print 'S_Code = S sources',len(cutout)
#cutout = cutout[ cutout['Total_flux']/cutout['Peak_flux'] < 1.25 + 3.1*(cutout['Peak_flux']/cutout['Isl_rms'])**-0.53]
#print 'Compact sources',len(cutout)


# Simply nearest neighbour match
matched = match_catalogues(lotssdr2,cutout,5,'cutout')
lotssdr2=lotssdr2[~np.isnan(lotssdr2['cutout_separation'])]
print('After cross match',len(lotssdr2))

# Do a few different fitting techniques (see https://machinelearningmastery.com/robust-regression-for-machine-learning-in-python/)


X = lotssdr2['cutout_Total_flux']*1000.0
X=X.reshape(len(X),1)
y = lotssdr2['Total_flux']
# define a uniform grid across the input domain
xaxis = arange(X.min(), X.max(), 0.01)

bestmae = 10000.0
for model in get_models():
    # plot the line of best fit
    modelname,gradient,intercept,mae,maestd = plot_best_fit(X, y, xaxis, model)
    if mae < bestmae:
        bestmodel,bestgradient,bestintercept,bestmae,bestmaestd = modelname,gradient,intercept,mae,maestd 

# plot the dataset
pyplot.scatter(X, y)
# show the plot
if np.min(X) < np.min(y):
    minval = np.min(X)*0.9
else:
    minval = np.min(y)*0.9
pyplot.loglog()
#pyplot.xlim(xmin=minval)
#pyplot.ylim(ymin=minval)
print(minval,'minval')
pyplot.title('Robust (and non-Robust) Regression')
pyplot.legend()
pyplot.savefig(infile.split('_')[0]+'_fitted.png')


ratios=lotssdr2['Total_flux']/lotssdr2['cutout_Total_flux']/1000.0
lotssdr2.write(infile.replace('.fits','cat.srl.matched.fits'),overwrite=True)

print('MODELNAME','FITTEDGRAD','FITTEDINTERCEPT','FIT_MEAN_ABS_ERROR','FIT_MEAN_ABS_ERROR_STD','NUM SOURCES FIT','MED SOURCE RATIOS','MEAN SOURCE RATIOS','STD SOURCE RATIOS')
print('BEST',bestmodel,bestgradient,bestintercept,bestmae,bestmaestd,len(y),np.median(ratios),np.mean(ratios),np.std(ratios))
print('MULTIPLY IMAGE BY',bestgradient)
#print 'Median,mean,std',np.median(ratios),np.mean(ratios),np.std(ratios)
#print 'Multiply image by ',np.median(ratios)
ds9file = open(infile.replace('.fits','cat.srl.matched.reg'),'w')
ds9file.write('# Region file format: DS9 version 4.1 \n')
ds9file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
ds9file.write('fk5\n')

for source in lotssdr2:
	ds9file.write('circle(%s,%s,20")\n'%(source['RA'],source['DEC']))
ds9file.close()

