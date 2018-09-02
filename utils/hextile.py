# hexagonal tiling of a fits image with fixed-radius circles
# return the central position list

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
    print 'Hex tiling image'
    # co-ords of bottom left of image
    ra_c,dec_c=w.wcs_pix2world(maxx/2,maxy/2,0)
    ra_factor=np.cos(dec_c*np.pi/180.0)
    ra_ll,dec_ll=w.wcs_pix2world(0,0,0)
    ra_lr,dec_lr=w.wcs_pix2world(maxx,0,0)
    ra_ul,dec_ul=w.wcs_pix2world(0,maxy,0)
    c_ll=SkyCoord(ra_ll*u.degree,dec_ll*u.degree,frame='icrs')
    c_lr=SkyCoord(ra_lr*u.degree,dec_lr*u.degree,frame='icrs')
    dra,ddec=[v.value for v in c_lr.spherical_offsets_to(c_ll)]
    nha=dra/hs
    print 'Number of hexes across',nha
    c_ul=SkyCoord(ra_ul*u.degree,dec_ul*u.degree,frame='icrs')
    dra,ddec=[v.value for v in c_ll.spherical_offsets_to(c_ul)]
    nhu=ddec/hs
    print 'Number of hexes up',nhu
    nha=int(0.5+nha)
    nhu=int(0.5+nhu)
    for j in range(nhu):
        for i in range(nha):
            ra_p=ra_lr+(i*hs+(j % 2)*0.5*hs)/ra_factor
            dec_p=dec_ll+((j+0.25)*hs)
            c_pos=SkyCoord(ra_p*u.degree,dec_p*u.degree,frame='icrs')
            pos.append((c_pos.ra.value,c_pos.dec.value))
    return ra_factor,pos

def plotcircle(ra,dec,xsize,ysize,color):
    circle1=Ellipse((ra,dec),width=2*xsize,height=2*ysize,angle=0,color=color,alpha=0.2)
    plt.gcf().gca().add_artist(circle1)
    plt.scatter(ra,dec)

if __name__=='__main__':
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    ra_factor,pos=hextile('image_ampphase1_di.int.restored.fits',0.5)
    for p in pos:
        ra=p[0]
        dec=p[1]
        plotcircle(ra,dec,0.5/ra_factor,0.5,'blue')
    plt.gca().invert_xaxis()
    plt.xlabel('RA')
    plt.ylabel('Dec')

    plt.show()
