#!/usr/bin/env python

import os,sys
import glob
import pyfits
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from auxcodes import sepn, getposim
import argparse
import time
import random

# Code to concatenate catalogs made from mosaics, remove repeated sources and to manipulate catalog entries to the final catalog format.

## To do
# -- Need to finalise the column choices
# -- We probably want to check the masks for each image use those mask islands to determine which mosaic to select the source from as perhaps there  are very extneded sources approximately midway between pointings (especially for gaussian components list)
# -- We probably want an entry in the catalog to say if the source was in mask used for the final deconvolution
# -- Needs speeding up -- each pointing separately
# -- update the errors that are put on the RA and DEC (could use astrometry error maps?)
# -- update the source extension definition

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

def find_median_astrometry(astromaps,pointingra,pointingdec):

    foundpointing = False
    for astromap in astromaps:
        ra,dec=getposim('%s/astromap.fits'%astromap)
        if sepn(ra*deg2rad,dec*deg2rad,pointingra,pointingdec)*rad2deg < 0.6:        
            foundpointing = True
            f = pyfits.open('%s/astromap.fits'%astromap)
            _,_,ys,xs=f[0].data.shape
            # Use central 20% of the image
            subim=f[0].data[0][0][ys/2-ys/5:ys/2+ys/5,xs/2-xs/5:xs/2+xs/5].flatten()
    if foundpointing:
        subim = np.nan_to_num(subim)
        return np.median(subim)
    else:
        print 'Cant find astrometry near %s, %s'%(pointingra*rad2deg,pointingdec*rad2deg)
        sys.exit(0)
   

def concat_catalogs(cats,outconcatcat):
    # Use the first catalogue as a dummy and then just update the entries

    f = pyfits.open(cats[0])

    nrows = f[1].data.shape[0]
    for cat in cats[1:]:
        f2 = pyfits.open(cat)
        nrows += f2[1].data.shape[0]
        f2.close()
    print nrows
    hdu = pyfits.BinTableHDU.from_columns(f[1].columns, nrows=nrows)

    nrows1 = f[1].data.shape[0]
    for cat in cats[1:]:
        f2 = pyfits.open(cat)
        nrows2 = nrows1 + f2[1].data.shape[0]
        for colname in f[1].columns.names:
            hdu.data[colname][nrows1:nrows2] = f2[1].data[colname]
        f2.close()
        nrows1 += f2[1].data.shape[0]

    hdu.writeto(outconcatcat)

def find_pointing_coords(directories):

    mosdirectories = args.mosdirectories
    mosaiccats = []
    for directory in mosdirectories:
        dircats = glob.glob('%s/*cat.fits'%directory)
        for mosaiccat in dircats:
            mosaiccats.append(mosaiccat)

    # Determine all pointing coordinates
    pointingras = np.array([])
    pointingdecs = np.array([])
    for mosaiccat in mosaiccats:
        pointing = mosaiccat.replace('.cat.fits','-blanked.fits')
        f = pyfits.open(pointing)
        pointingras = np.append(pointingras,f[0].header['CRVAL1']*deg2rad)
        pointingdecs = np.append(pointingdecs,f[0].header['CRVAL2']*deg2rad)
        f.close()
        
    return pointingras,pointingdecs,mosaiccats

def filter_catalogs(args,pointingras,pointingdecs,mosaiccat,outname,dessourcenums,cattype):

    pointdirectories = args.pointdirectories

    sourceids = np.array([])
    sourceresolved = np.array([])
    sourcera = np.array([])
    e_sourcera = np.array([])
    e_sourcera_tot = np.array([])
    sourcedec = np.array([])
    e_sourcedec = np.array([])
    e_sourcedec_tot = np.array([])
    sint = np.array([])
    e_sint = np.array([])
    e_sint_tot = np.array([])
    speak = np.array([])
    e_speak = np.array([])
    e_speak_tot = np.array([])
    maj = np.array([])
    e_maj = np.array([])
    smin = np.array([])
    e_smin = np.array([])
    pa = np.array([])
    e_pa = np.array([])
    rms_noise = np.array([])
    stype = np.array([])
    mosaic_identifier  = np.array([])
    gausid = np.array([])
    islid = np.array([])
    sourcenum = np.array([])

    
    pointing = mosaiccat.replace('.cat.fits','-blanked.fits')
    if cattype == 'gaus':
        mosaiccat = glob.glob(mosaiccat.replace('mosaic.cat.fits','mosaic-blanked_pybdsm/*/catalogues/mosaic-blanked.pybdsm.gaul.FITS'))[0]
    cat = pyfits.open(mosaiccat)
    
    f = pyfits.open(pointing)
    rapointing = f[0].header['CRVAL1']*deg2rad
    decpointing = f[0].header['CRVAL2']*deg2rad
    f.close()

    numsources = len(cat[1].data['RA'])

    closepointingindex = np.where(sepn(pointingras,pointingdecs,rapointing,decpointing)*rad2deg < 5.0)

    astromed = find_median_astrometry(pointdirectories,rapointing,decpointing)

    keepindices = []
    time1 = time.time()
    for i in range(0,numsources):
        
        if dessourcenums == None:
            allsep = sepn(pointingras[closepointingindex],pointingdecs[closepointingindex],cat[1].data['RA'][i]*deg2rad,cat[1].data['DEC'][i]*deg2rad)
            centsep =  sepn(rapointing,decpointing,cat[1].data['RA'][i]*deg2rad,cat[1].data['DEC'][i]*deg2rad)
            if min(allsep) != centsep:
                continue
            else:
                keepindices.append(i)
        else:
            if cat[1].data['Source_id'][i] in dessourcenums:
                keepindices.append(i)
            else:
                continue
            
        sc=SkyCoord(cat[1].data['RA'][i]*deg2rad*u.rad,cat[1].data['DEC'][i]*deg2rad*u.rad,frame='icrs')
        s=sc.to_string(style='hmsdms',sep='',precision=2)
        identity = str('ILTJ'+s).replace(' ','')[:-1]

        sourceids = np.append(sourceids,identity)
        mosaic_identifier = np.append(mosaic_identifier,mosaiccat.split('/')[-2])
        fluxratio = cat[1].data['Total_flux'][i]/cat[1].data['Peak_flux'][i]
        snr  = cat[1].data['Peak_flux'][i]/cat[1].data['Isl_rms'][i]

        # Some equation to figure out if the source is resolved -- leave these dummy values for now.
        if fluxratio > (1.483 + 1000.4/(snr**3.94)):
            #if fluxratio > (1.50341355 + 1.78467767/(snr**0.78385826)):
            sourceresolved = np.append(sourceresolved,'R')
        else:
            sourceresolved = np.append(sourceresolved,'U')

    print 'Keeping %s sources for %s -- took %s'%(len(keepindices),pointing,time.time()-time1)

    sourcera = np.append(sourcera,cat[1].data[keepindices]['RA'])
    e_sourcera = np.append(e_sourcera,cat[1].data[keepindices]['E_RA'])
    e_sourcera_tot = np.append(e_sourcera_tot,(cat[1].data[keepindices]['E_RA']**2.0 + (astromed*arcsec2deg)**2.0)**0.5) #$ ADD SOME ERROR TO THE SOURCE POSITIONS
    
    sourcedec = np.append(sourcedec,cat[1].data[keepindices]['DEC'])
    e_sourcedec = np.append(e_sourcedec,cat[1].data[keepindices]['E_DEC'])
    e_sourcedec_tot = np.append(e_sourcedec_tot,(cat[1].data[keepindices]['E_DEC']**2.0 + (astromed*arcsec2deg)**2.0)**0.5) # ADD SOME ERROR TO THE SOURCE POSITIONS
    
    sint = np.append(sint,cat[1].data[keepindices]['Total_flux'])
    e_sint =np.append(e_sint,cat[1].data[keepindices]['E_Total_flux'])
    e_sint_tot =np.append(e_sint_tot,(cat[1].data[keepindices]['E_Total_flux']**2.0 + (cat[1].data[keepindices]['Total_flux']*0.2)**2.0)**0.5)
    
    speak = np.append(speak,cat[1].data[keepindices]['Peak_flux'])
    e_speak =np.append(e_speak,cat[1].data[keepindices]['E_Peak_flux'])
    e_speak_tot =np.append(e_speak_tot,(cat[1].data[keepindices]['E_Peak_flux']**2.0 + (cat[1].data[keepindices]['Peak_flux']*0.2)**2.0)**0.5)
    
    maj = np.append(maj,cat[1].data[keepindices]['Maj'])
    e_maj =np.append(e_maj,cat[1].data[keepindices]['E_Maj'])
    smin = np.append(smin,cat[1].data[keepindices]['Min'])
    e_smin = np.append(e_smin,cat[1].data[keepindices]['E_Min'])
    pa = np.append(pa,cat[1].data[keepindices]['PA'])
    e_pa = np.append(e_pa,cat[1].data[keepindices]['E_PA'])
    rms_noise = np.append(rms_noise,cat[1].data[keepindices]['Isl_rms'])
    stype = np.append(stype,cat[1].data[keepindices]['S_Code'])
    islid = np.append(islid,cat[1].data[keepindices]['Isl_id'])
    sourcenum = np.append(sourcenum,cat[1].data[keepindices]['Source_id'])
    
    col1 = pyfits.Column(name='Source_Name',format='24A',unit='',array=sourceids)
    col2 = pyfits.Column(name='RA',format='f8',unit='deg',array=sourcera)
    col3 = pyfits.Column(name='E_RA',format='f8',unit='arcsec',array=e_sourcera*deg2arcsec)
    col4 = pyfits.Column(name='E_RA_tot',format='f8',unit='arcsec',array=e_sourcera_tot*deg2arcsec)
    
    col5 = pyfits.Column(name='DEC',format='f8',unit='deg',array=sourcedec)
    col6 = pyfits.Column(name='E_DEC',format='f8',unit='arcsec',array=e_sourcedec*deg2arcsec)
    col7 = pyfits.Column(name='E_DEC_tot',format='f8',unit='arcsec',array=e_sourcedec_tot*deg2arcsec)
    
    col8 = pyfits.Column(name='Peak_flux',format='f8',unit='beam-1 mJy',array=speak*1000.0)
    col9 = pyfits.Column(name='E_Peak_flux',format='f8',unit='beam-1 mJy',array=e_speak*1000.0)
    col10 = pyfits.Column(name='E_Peak_flux_tot',format='f8',unit='beam-1 mJy',array=e_speak_tot*1000.0)
    
    col11 = pyfits.Column(name='Total_flux',format='f8',unit='mJy',array=sint*1000.0)
    col12 = pyfits.Column(name='E_Total_flux',format='f8',unit='mJy',array=e_sint*1000.0)
    col13 = pyfits.Column(name='E_Total_flux_tot',format='f8',unit='mJy',array=e_sint_tot*1000.0)
    
    #maj[np.where(sourceresolved=='U')] = np.nan
    #e_maj[np.where(sourceresolved=='U')] = np.nan
    #smin[np.where(sourceresolved=='U')] = np.nan
    #e_smin[np.where(sourceresolved=='U')] = np.nan
    #pa[np.where(sourceresolved=='U')] = np.nan
    #e_pa[np.where(sourceresolved=='U')] = np.nan
    
    col14 =  pyfits.Column(name='Maj',format='f8',unit='arcsec',array=maj*deg2arcsec)
    col15 =  pyfits.Column(name='E_Maj',format='f8',unit='arcsec',array=e_maj*deg2arcsec)
    
    col16 =  pyfits.Column(name='Min',format='f8',unit='arcsec',array=smin*deg2arcsec)
    col17 =  pyfits.Column(name='E_Min',format='f8',unit='arcsec',array=e_smin*deg2arcsec)
    
    col18 =  pyfits.Column(name='PA',format='f8',unit='deg',array=pa)
    col19 =  pyfits.Column(name='E_PA',format='f8',unit='deg',array=e_pa)

    #col20 = pyfits.Column(name='Resolved',format='1A',unit='',array=sourceresolved)
    
    col20 = pyfits.Column(name='Isl_rms',format='f8',unit='beam-1 mJy',array=rms_noise*1000.0)
    col21 = pyfits.Column(name='S_Code',format='1A',unit='',array=stype)
    
    col22 = pyfits.Column(name='Mosaic_ID',format='8A',unit='',array=mosaic_identifier)
    
    col23 = pyfits.Column(name='Isl_id',format='I8',unit='',array=islid)
    col24 = pyfits.Column(name='Source_id',format='I8',unit='',array=sourcenum)
    
    if cattype == 'gaus':
        col25 = pyfits.Column(name='Gaus_id',format='I8',unit='',array=gausid)

    if cattype == 'srl':    
        cols = pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24])
    if cattype == 'gaus':
        cols = pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25])
        
    tbhdu = pyfits.BinTableHDU.from_columns(cols)

    if cattype == 'gaus':
        regionfile = open('%s.gaus.reg'%outname,'w')
    if cattype == 'srl':
        regionfile = open('%s.srl.reg'%outname,'w')
    regionfile.write('# Region file format: DS9 version 4.0\n')
    regionfile.write('global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    regionfile.write('fk5\n')
    for i in range(0,len(sourceids)):
	if not np.isnan(maj[i]):
            regionfile.write('ellipse(%s,%s,%s,%s,%s)\n'%(sourcera[i],sourcedec[i],maj[i],smin[i],pa[i]+90))
	else:
            regionfile.write('box(%s,%s,5.0",5.0",0.0)\n'%(sourcera[i],sourcedec[i]))
    regionfile.close()

    prihdr = pyfits.Header()
    prihdr['NAME'] = outname
    prihdu = pyfits.PrimaryHDU(header=prihdr)
    tbhdulist = pyfits.HDUList([prihdu, tbhdu])
    if cattype == 'srl':
        outcat = outname +'.srl.fits'
    if cattype == 'gaus':
        outcat = outname +'.gaus.fits'
    tbhdulist.writeto(outcat)

    return sourcenum,outcat

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Concatenate ddf-pipeline mosaic directories')
    parser.add_argument('--mosdirectories', metavar='D', nargs='+',help='mosaic directory name')
    parser.add_argument('--pointdirectories', metavar='D', nargs='+',help='pointing directory name')
    args = parser.parse_args()

    pointingras,pointingdecs,mosaiccats = find_pointing_coords(args)

    srlcatnames = []
    gauscatnames = []
    random.shuffle(mosaiccats)
    for mosaiccat in mosaiccats:
        print 'Working on %s'%mosaiccat
        outname = mosaiccat.split('/')[-2] + 'cat'
        if not os.path.exists(outname +'.srl.fits'):
            pointingsourcenums,srlcat = filter_catalogs(args,pointingras,pointingdecs,mosaiccat,outname,None,'srl')
            pointingsourcenums,gauscat = filter_catalogs(args,pointingras,pointingdecs,mosaiccat,outname,pointingsourcenums,'gaus')
        else:
            srlcat = outname +'.srl.fits'
            gauscat = outname +'.gaus.fits'
            
        srlcatnames.append(srlcat)
        gauscatnames.append(gauscat)
    print 'Concatenating %s files'%len(srlcatnames)
    concat_catalogs(srlcatnames,'LOFAR_HBA_T1_DR1_catalog_v0.1.srl.fits')
    concat_catalogs(gauscatnames,'LOFAR_HBA_T1_DR1_catalog_v0.1.gaus.fits')
