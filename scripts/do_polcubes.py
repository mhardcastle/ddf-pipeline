import numpy as np
from pipeline import *
from pyrap.tables import table
from auxcodes import getpos
from astropy import units as u
from astropy.coordinates import Angle

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
    print 'Output file shape is', newdata.shape
    for i,h in enumerate(hdus):
        if i==0:
            chb=0
        else:
            chb=sum(chans[:i])
        newdata[chb:chb+chans[i],:,:,:]=h[0].data

    ohdu=hdus[0]
    ohdu[0].data=newdata
    ohdu[0].header['NAXIS4']=np.sum(chans)
    hdus[0].writeto(outfile,clobber=True)

def do_polcubes(colname,
                CurrentDDkMSSolName,
                uvrange,
                ddf_kw,
                options=None,catcher=None):

    o=options

    mslist=o['full_mslist']
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    for i in range(0,len(filenames)):
        filename = filenames[i]
        freqs=[]
        t = pt.table(filename+'/SPECTRAL_WINDOW', readonly=True, ack=False)
        chanfreqs=t[0]['CHAN_FREQ']
        for freq in chanfreqs:
            if freq not in freqs:
                freqs.append(freq)
        channels=len(freqs)

        ThisImageName = 'image_full_low_QU_Cube%s'%i

        ddf_image(ThisImageName,filename,
                  cleanmode='SSD',ddsols=CurrentDDkMSSolName,
                  applysols='AP',
                  polcubemode=True,
		  AllowNegativeInitHMP=True,
                  majorcycles=0,robust=o['low_robust'],
                  colname=colname,use_dicomodel=False,
                  uvrange=uvrange,beamsize=o['low_psf_arcsec'],
                  imsize=2500,cellsize=o['low_cell'],peakfactor=0.001,
                  smooth=True,automask=True,automask_threshold=5,normalization=o['normalize'][2],channels=channels,
                  startchan=0,endchan=channels,options=o,
                  catcher=catcher)

    outfile='image_full_low_QU.cube.dirty.fits'
    if os.path.isfile(outfile):
        warn('Uncorrected cube file already exists, not making it')
    else:
        report('Making uncorrected cube')
        freqs,hdus=get_freqs_hdus('image_full_low_QU_Cube*.cube.dirty.fits')
        make_cube(freqs,hdus,outfile)

    outfile='image_full_low_QU.cube.dirty.corr.fits'
    if os.path.isfile(outfile):
        warn('Corrected cube file already exists, not making it')
    else:
        freqs,hdus=get_freqs_hdus('image_full_low_QU_Cube*.cube.dirty.corr.fits')
        report('Making corrected cube')
        make_cube(freqs,hdus,outfile)
            

