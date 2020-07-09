from __future__ import print_function
from __future__ import division
# Code to determine what facet a source is in

from builtins import range
from builtins import object
from past.utils import old_div
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import re

def point_inside_polygon(x,y,poly):
    # code from http://www.ariel.com.au/a/python-point-int-poly.html
    n = len(poly)
    inside =False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = old_div((y-p1y)*(p2x-p1x),(p2y-p1y))+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

class RegPoly(object):
    ''' Code for manipulating a region file as a list of polygons '''
    def coordconv(self,ra,dec):
        ''' return the original co-ords and the spherical offset
        The former is useful to convert co-ordinates to a standard form '''
        c=SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='icrs')
        return (c.ra.value,c.dec.value),[v.value for v in self.cref.spherical_offsets_to(c)]

    def __init__(self,inreg,cra,cdec):
        '''Convert a region to a polygon using a co-ordinate system relative
        to a specified reference point, ideally the image centre. Use
        astropy co-ordinates to avoid issues with co-ordinate
        singularities
        '''
        
        self.inreg=inreg
        self.cra=cra
        self.cdec=cdec
        self.cref=SkyCoord(ra=cra*u.degree, dec=cdec*u.degree,frame='icrs')
        
        lines=open(inreg).readlines()
        inpoly=False
        # clist and llist are in the offset co-ordinates
        # oclist and ollist are in original ra/dec co-ordinates, e.g. for plotting, but after standard form conversion
        clist=[]
        oclist=[]
        llist=[]
        ollist=[]
        for l in lines:
            if l[0:4]=='line':
                if not inpoly:
                    coords=[]
                    ocoords=[]
                    inpoly=True
                bits=l[5:].split(',')
                ra=float(bits[0])
                dec=float(bits[1])
                result=self.coordconv(ra,dec)
                ra,dec=result[0]
                dra,ddec=result[1]
                ocoords.append((ra,dec))
                coords.append((dra,ddec))
            else:
                if inpoly:
                    clist.append(coords)
                    oclist.append(ocoords)
                    inpoly=False
                if l[0:5]=='point':
                    m=re.match(r"point\((.*),(.*)\).*text=\{\[(.*)\]\}",l)
                    if m is None:
                        raise RuntimeError('Failed to parse region label %s' % l)
                    ra=float(m.group(1))
                    dec=float(m.group(2))
                    result=self.coordconv(ra,dec)
                    ra,dec=result[0]
                    dra,ddec=result[1]
                    ollist.append((ra,dec,m.group(3)))
                    llist.append((dra,ddec,m.group(3)))

        if inpoly:
            clist.append(coords)
            oclist.append(ocoords)

        bbox=[]
        for poly in clist:
            a=np.array(poly)
            xmin=np.min(a[:,0])
            xmax=np.max(a[:,0])
            ymin=np.min(a[:,1])
            ymax=np.max(a[:,1])
            bbox.append([xmin,xmax,ymin,ymax])
        self.bbox=np.array(bbox)
            
        # make the lists
        self.clist=clist
        self.oclist=oclist
        self.llist=llist
        self.ollist=ollist
        self.label_poly()
        self.labels_to_integers()
        
    def label_poly(self):
        # labels end up scrambled wrt their polygons, so we unscramble them
        llist=self.llist
        plab=[None]*len(llist)
        for dra,ddec,label in llist:
            i=self.which_poly(dra,ddec,convert=False)
            if i is None:
                raise RuntimeError('Failed to locate a label! (%s)' % label)
            else:
                plab[i]=label
        self.plab=plab

    def which_poly(self,ra,dec,convert=True):
        if convert:
            dra,ddec=self.coordconv(ra,dec)[1]
        else:
            dra=ra
            ddec=dec
        check=((dra>=self.bbox[:,0]) & (dra<=self.bbox[:,1]) & (ddec>=self.bbox[:,2]) & (ddec<=self.bbox[:,3]))
        for i,poly in enumerate(self.clist):
            if not check[i]:
                continue
            if point_inside_polygon(dra,ddec,poly):
                return i
        return None

    def labels_to_integers(self):
        plab_int=[]
        for p in self.plab:
            m=re.match(r".*_S(.*)",p)
            if m is None:
                raise RuntimeError('Label no regexp match: %s' % p)
            plab_int.append(int(m.group(1)))
        self.plab_int=plab_int

    def add_facet_labels(self,t):
        ''' Add integer labels to an astropy table t '''
        facets=[]
        dra,ddec=self.coordconv(np.array(t['RA']),np.array(t['DEC']))[1] # strip units
        
        for i in range(len(t)):
            poly=self.which_poly(dra[i],ddec[i],convert=False)
            if poly is not None:
                facets.append(self.plab_int[poly])
            else:
                facets.append(-1)
        t['Facet']=facets
        return t
        
def plot_offsets(t,poly,color):
    import matplotlib.pyplot as plt
    basesize=10
    rarange=(np.min(t['RA']),np.max(t['RA']))
    decrange=(np.min(t['DEC']),np.max(t['DEC']))
    mdec=np.mean(decrange)
    xstrue=(rarange[1]-rarange[0])*np.cos(mdec*np.pi/180.0)
    ystrue=decrange[1]-decrange[0]
    plt.figure(figsize=(old_div(basesize*xstrue,ystrue), basesize))
    plt.xlim(rarange)
    plt.ylim(decrange)
    plt.xlabel('RA')
    plt.ylabel('DEC')

    print(np.mean(t['FIRST_dRA']),np.mean(t['FIRST_dDEC']))
    
    for p in poly:
        x=[pt[0] for pt in p]
        y=[pt[1] for pt in p]
        plt.plot(x,y,color='black',ls=':')

    mra=[]
    mdec=[]
    mdra=[]
    mddec=[]
    maxreg=np.max(t['Facet'])
    for i in range(0,int(maxreg)+1):
        ts=t[t['Facet']==i]
        if len(ts)>5:
            mra.append(np.mean(ts['RA']))
            mdec.append(np.mean(ts['DEC']))
            mdra.append(np.mean(ts['FIRST_dRA']))
            mddec.append(np.mean(ts['FIRST_dDEC']))
            plt.quiver([mra[-1]]*len(ts),[mdec[-1]]*len(ts),ts['FIRST_dRA'],ts['FIRST_dDEC'],units = 'xy', angles='xy', scale=1.0,color=color,alpha=0.03)
            print(i,len(ts),mra[-1],mdec[-1],mdra[-1],mddec[-1])

    plt.gca().invert_xaxis()
    plt.quiver(mra,mdec,mdra,mddec,units = 'xy', angles='xy', scale=1.0,color=color)
    plt.quiver(np.mean(t['RA']),np.mean(t['DEC']),1.0,0.0,units = 'xy', angles='xy', scale=1.0,color='green')
    plt.text(np.mean(t['RA']),np.mean(t['DEC']),'1 arcsec',color='green')

def label_table(t,regfile,cra,cdec):
    ''' convenience function to label a fits table based on a region file '''
    r=RegPoly(regfile,cra,cdec)
    r.add_facet_labels(t)
    return t

#if __name__=='__main__':
#
#    do_plot_facet_offsets('image_full_ampphase1m.cat.fits_FIRST_match_filtered.fits','image_full_ampphase1m.tessel.reg')
#    plt.show()

