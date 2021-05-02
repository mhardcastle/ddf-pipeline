from __future__ import print_function
from __future__ import absolute_import
from builtins import range
import numpy as np
from pipeline import ddf_image,ddf_shift
from pyrap.tables import table
from astropy.io import fits
import glob
from auxcodes import warn

def update_frequencies(rootname,freq):
    g=glob.glob(rootname+'*.fits')
    for f in g:
        with fits.open(f) as hdu:
            if 'CRVAL4' in hdu[0].header and hdu[0].header['CRVAL4']!=freq:
                warn('Updating FITS header for %s to freq of %f' % (f,freq))
                hdu[0].header['CRVAL4']=freq
                hdu[0].header['RESTFRQ']=freq
                hdu.writeto(f,overwrite=True)
            del hdu[0].data

def do_spectral_restored(colname,
                         CurrentMaskName,
                         CurrentBaseDicoModelName,
                         CurrentDDkMSSolName,
                         uvrange,
                         ddf_kw,
                         facet_offset_file,
                         options=None,catcher=None):
    o=options
    if o['centralfreqs']:
        CentralFreqs=np.array(o['centralfreqs'])*1e6
    else:
        # old hard-wired behaviour
        CentralFreqs=np.array([128.02581787109375, 143.65081787109375, 160.25238037109375])*1e6

    mslist=o['full_mslist']
    filenames=[l.strip() for l in open(mslist,'r').readlines()]

    
    fMS=[table("%s::SPECTRAL_WINDOW"%MSName,ack=False).getcol("CHAN_FREQ").mean() for MSName in filenames]
    
    # f0=np.mean(fMS[0:9])/1e6
    # f1=np.mean(fMS[9:17])/1e6
    # f2=np.mean(fMS[17::])/1e6
    # CentralFreqs=np.array([f0,f1,f2])*1e6

    df=np.array(fMS).reshape((-1,1))-CentralFreqs.reshape((1,-1))
    iBandMapping=np.argmin(np.abs(df),axis=1)
    print(df.shape)
    print(iBandMapping)

    print()
    mslistnames=[]
    for iBand in range(CentralFreqs.size):
        ind=np.where(iBandMapping==iBand)[0]
        print("iBand %i -> %i ms"%(iBand,np.count_nonzero(iBandMapping==iBand)))
        MSListFileName="mslist_iBand_%i.txt"%iBand
        f=open(MSListFileName,"w")
        for iMSName in ind:
            f.write("%s\n"%filenames[iMSName])
        f.close()
        mslistnames.append(MSListFileName)
        freq=CentralFreqs[iBand]
        
        ThisImageName='image_full_ampphase_di_m.NS_Band%i'%iBand
        
        ddf_image(ThisImageName,
                  MSListFileName,
                  cleanmask=CurrentMaskName,
                  reuse_psf=False,
                  cleanmode='SSD',
                  ddsols=CurrentDDkMSSolName,
                  applysols=o['apply_sols'][6],
                  majorcycles=0,robust=o['final_robust'],
                  colname=colname,use_dicomodel=True,
                  dicomodel_base=CurrentBaseDicoModelName,
                  AllowNegativeInitHMP=True,
                  peakfactor=0.001,automask=True,automask_threshold=o['thresholds'][2],
                  normalization=o['normalize'][1],uvrange=uvrange,smooth=True,
                  apply_weights=o['apply_weights'][2],catcher=catcher,RMSFactorInitHMP=1.,options=o,
                  **ddf_kw)

        if o['method'] is not None:
            ddf_shift(ThisImageName,facet_offset_file,options=o,catcher=catcher,dicomodel=CurrentBaseDicoModelName+'.DicoModel')

        update_frequencies(ThisImageName,freq)
            
    return mslistnames # cache for these can be added to the tidy up
                       # at pipeline end
