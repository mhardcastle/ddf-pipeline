from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
# hexagonal tiling of a fits image with fixed-radius circles
# return the central position list

from builtins import range
from past.utils import old_div
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from auxcodes import flatten
import numpy as np

def hextile(image,radius):

    pos=[]
    hs=radius*np.sqrt(3)
    hdus=fits.open(image)
    hdu=flatten(hdus)
    maxy,maxx=hdu.data.shape
    w=WCS(hdu.header)
    print('Hex tiling image')
    # co-ords of bottom left of image
    ra_c,dec_c=w.wcs_pix2world(old_div(maxx,2),old_div(maxy,2),0)
    ra_factor=np.cos(dec_c*np.pi/180.0)
    ra_ll,dec_ll=w.wcs_pix2world(0,0,0)
    ra_lr,dec_lr=w.wcs_pix2world(maxx,0,0)
    ra_ul,dec_ul=w.wcs_pix2world(0,maxy,0)
    c_c=SkyCoord(ra_c*u.degree,dec_c*u.degree,frame='icrs')
    c_ll=SkyCoord(ra_ll*u.degree,dec_ll*u.degree,frame='icrs')
    c_lr=SkyCoord(ra_lr*u.degree,dec_lr*u.degree,frame='icrs')
    dra,ddec=[v.value for v in c_c.spherical_offsets_to(c_ll)]
    nha=old_div(dra*2,hs)
    print('Number of hexes across',nha)
    c_ul=SkyCoord(ra_ul*u.degree,dec_ul*u.degree,frame='icrs')
    dra,ddec=[v.value for v in c_c.spherical_offsets_to(c_ul)]
    nhu=old_div(2*ddec,hs)
    print('Number of hexes up',nhu)
    nha=int(0.5+nha)
    nhu=int(0.5+nhu)
    for j in range(nhu):
        for i in range(nha):
            xc=old_div((1.0*maxx*(i+(j % 2)*0.5)),nha)
            yc=old_div((maxy*(j+0.5)),nhu)
            ra_p,dec_p=w.wcs_pix2world(xc,yc,0)
            pos.append((float(ra_p),float(dec_p)))
    return ra_factor,pos

def plotcircle(ra,dec,xsize,ysize,color):
    circle1=Ellipse((ra,dec),width=2*xsize,height=2*ysize,angle=0,color=color,alpha=0.2)
    plt.gcf().gca().add_artist(circle1)
    plt.scatter(ra,dec)

if __name__=='__main__':
    import os
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    from astropy.table import Table
    ra_factor,pos=hextile('image_ampphase1_di.int.restored.fits',0.5)
    if ra_factor<0.1: ra_factor=0.5
    for p in pos:
        ra=p[0]
        dec=p[1]
        plotcircle(ra,dec,0.5/ra_factor,0.5,'blue')
    plt.gca().invert_xaxis()
    plt.xlabel('RA')
    plt.ylabel('Dec')
    catfile='image_full_ampphase_di_m.NS.offset_cat.fits'
    if os.path.isfile(catfile):
        t=Table.read(catfile)
        plt.scatter(t['RA'],t['DEC'],alpha=0.1,marker='.')
    plt.show()
