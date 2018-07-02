
import numpy as np
from pipeline import *
from pyrap.tables import table
from auxcodes import getpos
from astropy import units as u
from astropy.coordinates import Angle

def do_polcubes(colname,
                         CurrentDDkMSSolName,
                         uvrange,changroups,
                         ddf_kw,
                         options=None,catcher=None):
    o=options

    mslist=o['full_mslist']
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    freqs=[]
    for filename in filenames:
        t = pt.table(filename+'/SPECTRAL_WINDOW', readonly=True, ack=False)
        chanfreqs=t[0]['CHAN_FREQ']
        for freq in chanfreqs:
            if freq not in freqs:
                freqs.append(freq)
    channels=len(freqs)
    chanranges = np.linspace(0,channels,changroups+1)

    # Make 4 cubes each with a quarter of the bandwidth.
    for i in range(0,len(chanranges)-1):
        ThisImageName = 'image_full_low_QU_Cube%s'%i

        ddf_image(ThisImageName,o['full_mslist'],
                  cleanmode='SSD',ddsols=CurrentDDkMSSolName,
                  applysols='AP',
                  polcubemode=True,
		  AllowNegativeInitHMP=True,
                  majorcycles=0,robust=o['low_robust'],
                  colname=colname,use_dicomodel=False,
                  uvrange=uvrange,beamsize=o['low_psf_arcsec'],
                  imsize=2500,cellsize=o['low_cell'],peakfactor=0.001,
                  smooth=True,automask=True,automask_threshold=5,normalization=o['normalize'][2],channels=int(chanranges[i+1]-chanranges[i]),
                  startchan=int(chanranges[i]),endchan=int(chanranges[i+1]),options=o,
                  catcher=catcher)



