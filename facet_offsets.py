# Code to determine what facet a source is in

from astropy.io import fits
from astropy.table import Table
import numpy as np
import re
import matplotlib.pyplot as plt

def region_to_poly(inreg):
    lines=open(inreg).readlines()

    inpoly=False
    clist=[]
    coords=[]
    llist=[]
    for l in lines:
        if l[0:4]=='line':
            if not inpoly:
                coords=[]
                inpoly=True
            bits=l[5:].split(',')
            coords.append((float(bits[0]),float(bits[1])))
    #        coords.append((bits[2],bits[3].split(')')[0]))
        else:
            if inpoly:
                clist.append(coords)
                inpoly=False
            if l[0:5]=='point':
                m=re.match(r"point\((.*),(.*)\).*text=\{\[(.*)\]\}",l)
                if m is None:
                    print 'Failed to match:',l
                llist.append((float(m.group(1)),float(m.group(2)),m.group(3)))

    if inpoly:
        clist.append(coords)

    return clist,llist
    
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
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def which_poly(x,y,plist):
    for i,poly in enumerate(plist):
        if point_inside_polygon(x,y,poly):
            return i
    return None

def assign_labels_to_poly(plist,llist):
    # labels end up scrambled wrt their polygons, so we unscramble them
    plabel=[None]*len(plist)
    for ra,dec,label in llist:
        i=which_poly(ra,dec,plist)
        if i is None:
            print 'Failed to locate a label!'
            stop
        else:
            plabel[i]=label
    return plabel

def labels_to_integers(plab):
    plab_int=[]
    for p in plab:
        m=re.match(r".*_S(.*)",p)
        if m is None:
            print 'No match:',p
        plab_int.append(int(m.group(1)))
    return plab_int

def add_facet_labels(t,polys,pli):
    facets=[]
    for r in t:
        poly=which_poly(r['RA'],r['DEC'],polys)
        if poly is not None:
            facets.append(pli[poly])
        else:
            facets.append(-1)
    t['Facet']=facets

def plot_offsets(t,poly,color):

    basesize=10
    rarange=(np.min(t['RA']),np.max(t['RA']))
    decrange=(np.min(t['DEC']),np.max(t['DEC']))
    mdec=np.mean(decrange)
    xstrue=(rarange[1]-rarange[0])*np.cos(mdec*np.pi/180.0)
    ystrue=decrange[1]-decrange[0]
    plt.figure(figsize=(basesize*xstrue/ystrue, basesize))
    plt.xlim(rarange)
    plt.ylim(decrange)
    plt.xlabel('RA')
    plt.ylabel('DEC')

    print np.mean(t['FIRST_dRA']),np.mean(t['FIRST_dDEC'])
    
    for p in poly:
        x=[pt[0] for pt in p]
        y=[pt[1] for pt in p]
        plt.plot(x,y,color='black',ls=':')

#    lines=open(polyfile).readlines()
#    i=0
#    for l in lines:
#        if l[0:7]=='polygon':
#            i+=1
#            ts=t[t['Region']==i]
#            if len(ts)>0:
#                bits=l[8:-2].split(', ')
#                x=[float(b) for b in bits[::2]]
#                y=[float(b) for b in bits[1::2]]
#                plt.plot(x,y,color='black',ls=':')

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
            print i,len(ts),mra[-1],mdec[-1],mdra[-1],mddec[-1]

    plt.gca().invert_xaxis()
    plt.quiver(mra,mdec,mdra,mddec,units = 'xy', angles='xy', scale=1.0,color=color)
    plt.quiver(np.mean(t['RA']),np.mean(t['DEC']),1.0,0.0,units = 'xy', angles='xy', scale=1.0,color='green')
    plt.text(np.mean(t['RA']),np.mean(t['DEC']),'1 arcsec',color='green')

def do_plot_facet_offsets(t,regfile,savefig=None):

    if isinstance(t,str):
        t=Table.read(t)
    polys,labels=region_to_poly(regfile)
    plab=assign_labels_to_poly(polys,labels)
    pli=labels_to_integers(plab)

    add_facet_labels(t,polys,pli)
    plot_offsets(t,polys,'red')
    if savefig is not None:
        plt.savefig(savefig)

if __name__=='__main__':

    do_plot_facet_offsets('image_full_ampphase1m.cat.fits_FIRST_match_filtered.fits','image_full_ampphase1m.tessel.reg')
    plt.show()

