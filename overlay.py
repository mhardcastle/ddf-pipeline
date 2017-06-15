# overlay code ripped out of show_overlay's various incarnations

# the idea is to separate the plotting of the overlay from the
# determination of what maps to use

# lofarhdu must be a flattened (2D) LOFAR map. opthdu can be any type
# of optical map. rms_use must be the loc

import aplpy
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.table import Table

def find_noise(a):
    b=a.flatten()
    for i in range(10):
        m=np.nanmean(b)
        s=np.nanstd(b)
        b=b[b<(m+5.0*s)]
    return m,s

def find_noise_area(hdu,ra,dec,size):
    # ra, dec, size in degrees
    ysize,xsize=hdu[0].data.shape
    w=WCS(hdu[0].header)
    ras=[ra-size,ra-size,ra+size,ra+size]
    decs=[dec-size,dec+size,dec-size,dec+size]
    xv=[]
    yv=[]
    for r,d in zip(ras,decs):
        x,y=w.wcs_world2pix(r,d,0)
        xv.append(x)
        yv.append(y)
    xmin=int(min(xv))
    if xmin<0: xmin=0
    xmax=int(max(xv))
    if xmax>=xsize: xmax=xsize-1
    ymin=int(min(yv))
    if ymin<0: ymin=0
    ymax=int(max(yv))
    if ymax>=ysize: ymax=ysize-1
    print xmin,xmax,ymin,ymax
    subim=hdu[0].data[ymin:ymax,xmin:xmax]
    mean,noise=find_noise(subim)
    for i in range(5,0,-1):
        vmax=np.nanmedian(subim[subim>(mean+i*noise)])
        if not(np.isnan(vmax)):
            break
    return mean,noise,vmax
    

def show_overlay(lofarhdu,opthdu,ra,dec,size,firsthdu=None,rms_use=None,bmaj=None,bmin=None,bpa=None,title=None,save_name=None,plotpos=None,block=True,interactive=False,plot_coords=True,overlay_cat=None,lw=1.0,show_lofar=True,no_labels=False,show_grid=True,overlay_region=None,overlay_scale=1.0,circle_radius=None,coords_color='white',coords_lw=1,coords_ra=None,coords_dec=None,marker_ra=None,marker_dec=None,marker_color='white',marker_lw=3,noisethresh=1):

    try:
        from matplotlib.cbook import MatplotlibDeprecationWarning
        import warnings
        warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
    except:
        print 'Cannot hide warnings'

    print '========== Doing overlay at',ra,dec,'with size',size,'==========='
    print '========== Title is',title,'=========='

    if coords_ra is None:
        coords_ra=ra
        coords_dec=dec

    if show_lofar:
        lofarmax=np.nanmax(lofarhdu[0].data)
        if rms_use is None:
            rms_use=find_noise_area(lofarhdu,ra,dec,size)[1]
            print 'Using LOFAR rms',rms_use
        drlimit=1000
        print lofarmax/drlimit,rms_use*2.0
        minlevel=max([lofarmax/drlimit,rms_use*2.0])
        levels=minlevel*2.0**np.linspace(0,14,30)

    hdu=opthdu
    mean,noise,vmax=find_noise_area(hdu,ra,dec,size)
    print 'Optical parameters are',mean,noise,vmax
    f = aplpy.FITSFigure(hdu,north=True)
    print 'centring on',ra,dec,size
    f.recenter(ra,dec,width=size,height=size)
    f.show_colorscale(vmin=mean+noisethresh*noise,vmax=vmax,stretch='log')
    #f.show_colorscale(vmin=0,vmax=1e-3)
    #f.show_colorscale(vmin=0,vmax=1.0)
    if bmaj is not None:
        f._header['BMAJ']=bmaj
        f._header['BMIN']=bmin
        f._header['BPA']=bpa

    if show_lofar: f.show_contour(lofarhdu,colors='yellow',linewidths=lw, levels=levels)

    if firsthdu is not None:
        firstrms=find_noise_area(firsthdu,ra,dec,size)[1]
        print 'Using FIRST rms',firstrms
        firstlevels=firstrms*3*2.0**np.linspace(0,14,30)
        f.show_contour(firsthdu,colors='lightgreen',linewidths=lw, levels=firstlevels)

    if bmaj is not None:
        f.add_beam()
        f.beam.set_corner('bottom left')
        f.beam.set_edgecolor('red')
        f.beam.set_facecolor('white')

    if plot_coords:
        f.show_markers(coords_ra,coords_dec,marker='+',facecolor=coords_color,edgecolor=coords_color,linewidth=coords_lw,s=1500,zorder=100)

    if marker_ra is not None:
        f.show_markers(marker_ra,marker_dec,marker='x',facecolor=marker_color,edgecolor=marker_color,linewidth=marker_lw,s=1500,zorder=100)
        

    if plotpos is not None:
        if not isinstance(plotpos,list):
            plotpos=[(plotpos,'x'),]
        for t,marker in plotpos:
            if len(t)>0:
                f.show_markers(t['ra'],t['dec'],marker=marker,facecolor='white',edgecolor='white',linewidth=2,s=750)

    if circle_radius is not None:
        f.show_circles([ra,],[dec,],[circle_radius,],facecolor='none',edgecolor='cyan',linewidth=5)

    if overlay_cat is not None:
        t=overlay_cat
        f.show_ellipses(t['RA'],t['DEC'],t['Maj']*2/overlay_scale,t['Min']*2/overlay_scale,angle=90+t['PA'],edgecolor='red',linewidth=3,zorder=100)

    if overlay_region is not None:
        f.show_regions(overlay_region)

    if no_labels:
        f.axis_labels.hide()
        f.tick_labels.hide()
    if show_grid:
        f.add_grid()
        f.grid.show()
        f.grid.set_xspacing(1.0/60.0)
        f.grid.set_yspacing(1.0/60.0)

    if title is not None:
        plt.title(title)
    plt.tight_layout()
    if save_name is None:
        plt.show(block=block)
    else:
        plt.savefig(save_name)

    if interactive:
        def onclick(event):
            xp=event.xdata
            yp=event.ydata
            xw,yw=f.pixel2world(xp,yp)
            if event.button==2:
                print title,xw,yw
        fig=plt.gcf()

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
    return f
