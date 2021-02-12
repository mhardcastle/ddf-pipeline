from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import zip
from past.utils import old_div
import numpy as np
from pipeline import *
from pyrap.tables import table
from auxcodes import getpos,MSList
from astropy import units as u
from astropy.coordinates import Angle
import os

# Utility functions for making the cube

def get_freqs_hdus(filenames):
    hdus=[]
    freqs=[]
    g=glob.glob(filenames)
    for f in g:
        hdus.append(fits.open(f))
        header=hdus[-1][0].header
        freqs.append(header['CRVAL4'])

    freqs,hdus = (list(x) for x in zip(*sorted(zip(freqs, hdus), key=lambda pair: pair[0])))
    return freqs,hdus
    
def make_cube(freqs,hdus,outfile):

    chans=[]
    for h in hdus:
        ch,stokes,y,x=h[0].data.shape
        chans.append(ch)
        
    newdata=np.zeros((np.sum(chans),stokes,y,x),dtype=np.float32)
    print('Output file shape is', newdata.shape)
    for i,h in enumerate(hdus):
        if i==0:
            chb=0
        else:
            chb=sum(chans[:i])
        newdata[chb:chb+chans[i],:,:,:]=h[0].data

    ohdu=hdus[0]
    ohdu[0].data=newdata
    ohdu[0].header['NAXIS4']=np.sum(chans)
    hdus[0].writeto(outfile,overwrite=True)

def make_cube_stokes(freqs,hdus,outfile,stokesparam):

    chans=[]
    for h in hdus:
        ch,stokes,y,x=h[0].data.shape
        chans.append(ch)

    report('Original shape %s %s %s %s'%(np.sum(chans),stokes,y,x))
    newdata=np.zeros((np.sum(chans),y,x),dtype=np.float32)
    print('Output file shape is', newdata.shape)
    for i,h in enumerate(hdus):
        report('Doing planes from small cube %s for stokes %s'%(i,stokesparam))
        if i==0:
            chb=0
        else:
            chb=sum(chans[:i])
        newdata[chb:chb+chans[i],:,:]=h[0].data[:,stokesparam,:,:]

    ohdu=hdus[0]
    ohdu[0].data=newdata
    ohdu[0].header['NAXIS3']=np.sum(chans)
    hdus[0].writeto(outfile,overwrite=True)


def do_polcubes(colname,
                CurrentDDkMSSolName,
		uvrange,imageoutname,
		ddf_kw,
                beamsize,imsize,cellsize,robust,
                options,catcher):

    o=options

    m=MSList(o['full_mslist'])
    ufreqs=sorted(set(m.freqs))
    for i,freq in enumerate(ufreqs):
        print('Image %i: channel map for frequency %.3f MHz' % (i,old_div(freq,1e6)))
        # iterate over frequencies, finding all MS with the same values
        fmslist=[]
        for ms,f,chan in zip(m.mss,m.freqs,m.channels):
            if f==freq:
                fmslist.append(ms)
        mslistname='stokes-mslist-%i.txt' % i
        with open(mslistname,'w') as file:
            for ms in fmslist:
                file.write(ms+'\n')
        channels=len(chan)
                
        ThisImageName = '%s_QU_Cube%s'%(imageoutname,i)

        ddf_image(ThisImageName,mslistname,
                  cleanmode='SSD',ddsols=CurrentDDkMSSolName,
                  applysols=o['apply_sols'][6],
                  polcubemode=True,
		  AllowNegativeInitHMP=True,
                  majorcycles=0,robust=robust,
                  colname=colname,use_dicomodel=False,
                  uvrange=uvrange,beamsize=beamsize,
                  imsize=imsize,cellsize=cellsize,peakfactor=0.001,
                  smooth=True,automask=True,automask_threshold=5,normalization=o['normalize'][2],channels=channels,
                  startchan=0,endchan=channels,options=o,
                  catcher=catcher)

    if o['split_polcubes']:

        outfile='%s_StokesQ.cube.dirty.fits'%imageoutname
        if os.path.isfile(outfile) or os.path.isfile(outfile+'.fz'):
            warn('Uncorrected cube file already exists, not making it')
        else:
            report('Making uncorrected cube')
            freqs,hdus=get_freqs_hdus('%s_QU_Cube*.cube.dirty.fits'%imageoutname)
            outfile = '%s_StokesQ.cube.dirty.fits'%imageoutname
            make_cube_stokes(freqs,hdus,outfile,0)
            
        outfile='%s_StokesU.cube.dirty.fits'%imageoutname
        if os.path.isfile(outfile) or os.path.isfile(outfile+'.fz'):
            warn('Uncorrected cube file already exists, not making it')
        else:
            report('Making uncorrected cube')
            freqs,hdus=get_freqs_hdus('%s_QU_Cube*.cube.dirty.fits'%imageoutname) 
            outfile = '%s_StokesU.cube.dirty.fits'%imageoutname
            make_cube_stokes(freqs,hdus,outfile,1)

        outfile='%s_StokesQ.cube.dirty.corr.fits'%imageoutname
        if os.path.isfile(outfile) or os.path.isfile(outfile+'.fz'):
            warn('Uncorrected cube file already exists, not making it')
        else:
            report('Making uncorrected cube')
            freqs,hdus=get_freqs_hdus('%s_QU_Cube*.cube.dirty.corr.fits'%imageoutname)
            outfile = '%s_StokesQ.cube.dirty.corr.fits'%imageoutname
            make_cube_stokes(freqs,hdus,outfile,0)

        outfile='%s_StokesU.cube.dirty.corr.fits'%imageoutname
        if os.path.isfile(outfile) or os.path.isfile(outfile+'.fz'):
            warn('Uncorrected cube file already exists, not making it')
        else:
            report('Making uncorrected cube')
            freqs,hdus=get_freqs_hdus('%s_QU_Cube*.cube.dirty.corr.fits'%imageoutname)
            outfile = '%s_StokesU.cube.dirty.corr.fits'%imageoutname
            make_cube_stokes(freqs,hdus,outfile,1)
    else:

        outfile='%s_QU.cube.dirty.fits'%imageoutname
        if os.path.isfile(outfile) or os.path.isfile(outfile+'.fz'):
            warn('Uncorrected cube file already exists, not making it')
        else:
            report('Making uncorrected cube')
            freqs,hdus=get_freqs_hdus('%s_QU_Cube*.cube.dirty.fits'%imageoutname)
            make_cube(freqs,hdus,outfile)

        outfile='%s_QU.cube.dirty.corr.fits'%imageoutname
        if os.path.isfile(outfile) or os.path.isfile(outfile+'.fz'):
            warn('Corrected cube file already exists, not making it')
        else:
            freqs,hdus=get_freqs_hdus('%s_QU_Cube*.cube.dirty.corr.fits'%imageoutname)
            report('Making corrected cube')
            make_cube(freqs,hdus,outfile)
            
