
import numpy as np
from pipeline import *
from pyrap.tables import table
from auxcodes import getpos

def do_polcubes(colname,
                         CurrentDDkMSSolName,
                         uvrange,
                         ddf_kw,
                         options=None,catcher=None):
    o=options

    mslist=o['full_mslist']
    filenames=[l.strip() for l in open(mslist,'r').readlines()]

    name,ra,dec = getpos(filenames[0])
    
    # Make 4 cubes of 1250 4.5arcsec pixels (~1.5 degrees) which overlap 
    cuberas = [ra-0.75,ra-0.75,ra+0.75,ra+0.75]
    cubedecs = [dec+0.75,dec-0.75,dec+0.75,dec-0.75]

    for i in range(0,len(cuberas)):
        cubera,cubedec = cuberas[i],cubedecs[i]
        ThisImageName = 'image_full_low_QU_Cube%s'%i

        ddf_image(ThisImageName,o['full_mslist'],
                  cleanmode='SSD',ddsols=CurrentDDkMSSolName,
                  applysols='AP',
                  polcubemode=True,
		  AllowNegativeInitHMP=True,
                  majorcycles=0,robust=o['low_robust'],
                  colname=colname,use_dicomodel=False,
                  uvrange=uvrange,beamsize=o['low_psf_arcsec'],
                  imsize=1250,cellsize=o['low_cell'],peakfactor=0.001,
                  smooth=True,automask=True,automask_threshold=5,normalization=o['normalize'][2],
                  phasecenter=[cubera,cubedec],options=o,
                  catcher=catcher)

