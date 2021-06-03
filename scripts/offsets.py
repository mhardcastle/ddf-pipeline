#!/usr/bin/env python

# Compute optical offsets and do a shift
# Steps are:
# -1) download optical (or other) data -- do in parallel if possible, so should be done by the time this is run as a script
# 0) merge the downloads to a single catalogue
# 1) run pybdsm on the full-res image, filter?
# 2) label the catalogue
# 3) find offsets by matching 
# 4) fit to offset histograms
# 5) apply shift -- to be done by ddf-pipeline

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import matplotlib
matplotlib.use('Agg')
from auxcodes import report,run,warn,die,get_centpos,sepn
import requests
import os
from get_cat import get_cat
import glob
from astropy.table import Table, vstack, unique
import numpy as np
from scipy.optimize import curve_fit
import sys
from scipy.special import gammaln
from facet_offsets import RegPoly
from astropy.io import fits
from astropy.wcs import WCS
import pickle
try:
    import bdsf as bdsm
except ImportError:
    import lofar.bdsm as bdsm
degtorad=np.pi/180.0

import sys

def model(x,norm,sigma,offset,bl,radius=60):
    return old_div(bl*np.sqrt(radius**2.0-x**2.0),radius)+norm*np.exp(old_div(-(x-offset)**2.0,(2*sigma**2.0)))

class Offsets(object):
    def __init__(self,prefix,n=45,cellsize=1.5,imroot=None,fitmethod='mcmc'):
        self.prefix=prefix
        self.n=n
        self.chains=[]
        self.cellsize=cellsize
        self.imroot=imroot
        self.fitmethod=fitmethod
        if imroot is not None:
            self.read_regfile(imroot+'.tessel.reg')

    def read_regfile(self,regfile):
        cra,cdec=get_centpos()
        self.r=RegPoly(regfile,cra,cdec)
        self.polys=self.r.oclist
        self.labels=self.r.ollist
        self.plab=self.r.plab
        self.pli=self.r.plab_int

    def find_offsets(self,tf,ot,sep=1.0):
        self.dral=[]
        self.ddecl=[]
        self.lofar_table=tf
        for f in range(self.n):
            t=tf[tf['Facet']==f]
            if len(t)==0:
                print('No sources in facet',f)
                self.dral.append(None)
                self.ddecl.append(None)
                continue
            minra=np.min(t['RA'])
            maxra=np.max(t['RA'])
            mindec=np.min(t['DEC'])
            maxdec=np.max(t['DEC'])
            otf=ot[(ot['ra']>=(minra-sep/60.0)) & (ot['ra']<=(maxra+sep/60.0)) &
                   (ot['dec']>=(mindec-sep/60.0)) & (ot['dec']<=(maxdec+sep/60.0))]
            print('Facet %2i has %4i LOFAR sources and %6i comparison sources' % (f,len(t),len(otf)))

            dral=[]
            ddecl=[]

            for r in t:
                ra=r['RA']
                dec=r['DEC']
                dra=3600.0*(ra-otf['ra'])*np.cos(dec*np.pi/180.0)
                ddec=3600.0*(dec-otf['dec'])
                d2d=np.sqrt(dra**2.0+ddec**2.0)

                d2dmask = d2d<sep*60.0

                dral+=list(dra[d2dmask])
                ddecl+=list(ddec[d2dmask])

            self.dral.append((np.array(dral)).flatten())
            self.ddecl.append((np.array(ddecl)).flatten())

    def save_offsets(self):
        for i in range(self.n):
            np.save(self.prefix+'/dra-'+str(f)+'.npy',self.dral[i])
            np.save(self.prefix+'/ddec-'+str(f)+'.npy',self.ddecl[i])

    def load_offsets(self):
        self.dral=[]
        self.ddecl=[]
        for i in range(self.n):
            self.dral.append(np.load(self.prefix+'/dra-'+str(i)+'.npy'))
            self.ddecl.append(np.load(self.prefix+'/ddec-'+str(i)+'.npy'))
            
    def fit_chi2(self,h):
        height=np.median(h)
        norm=np.max(h)-height
        peak=self.bcenter[np.argmax(h)]
        popt,pcov=curve_fit(model,self.bcenter,h,[norm,0.5,peak,height],1.0+np.sqrt(h+0.75))
        return popt,np.sqrt(np.diagonal(pcov))

    def lnlike(self,X,h):
        if X[0]<0 or X[3]<0 or X[1]<0:
            return -np.inf
        mv=model(self.bcenter,*X)
        # Eq A3 of 3C305 paper; mv is mu, h is n
        lv=h*np.log(mv)-mv-gammaln(h+1)
        return np.sum(lv)

    def lnpost(self,parms,h):
        return self.lnprior(parms)+self.lnlike(parms,h)

    def lnprior(self,X):
        # gaussian norm, sigma, offset; baseline norm
        if X[0]<0 or X[3]<0 or X[1]<0:
            return -np.inf
        if X[1]>5:
            return -np.inf
        if np.abs(X[2])>5:
            return -np.inf
        #return -np.log(X[0])-np.log(X[3])
        return 0

    def fit_emcee(self,h):
        import emcee
        height=np.median(h)
        norm=np.max(h)-height
        peak=self.bcenter[np.argmax(h)]
        if np.abs(peak)>3.0:
            peak=0.0
        nwalkers=24
        ndim=4
        parms=[]
        for i in range(nwalkers):
            parms.append([norm+np.random.normal(0,0.5),0.5+np.random.normal(0,0.05),peak+np.random.normal(0,0.1),height+np.random.normal(0,2.0)])
        parms=np.array(parms)
        for i in (0,1,3):
            parms[:,i]=np.abs(parms[:,i])
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnpost, args=(h,))
        
        sampler.run_mcmc(parms,1000)
        chain=sampler.chain
        # find initial errors
        samples=chain[:, 200:, :].reshape((-1, ndim))
        samplest=samples.transpose()
        prange=np.percentile(samplest,(10,90),axis=1)
        wanted=[]
        for k in range(nwalkers):
            wparms=np.mean(chain[k,200:,:],axis=0)
            if np.all(wparms>prange[0]) and np.all(wparms<prange[1]):
                wanted.append(k)
        
        # now use only the walkers that didn't get lost
        chain=chain[wanted, :, :]
        samples=chain[:, 200:, :].reshape((-1, ndim))
        samplest=samples.transpose()

        self.chains.append(chain)
        means=np.mean(samplest,axis=1)
        errors=np.percentile(samplest,(16,84),axis=1)-means
        err=(errors[1]-errors[0])/2.0
        return means,err

    def fit_offsets(self,minv=-40,maxv=40,nbins=150):
        if self.fitmethod=='mcmc':
            fitfn=self.fit_emcee
        elif self.fitmethod=='chi2':
            fitfn=self.fit_chi2
        else:
            raise NotImplementedError('Fit method '+self.fitmethod)
        self.bins=np.linspace(minv,maxv,nbins+1)
        self.bcenter=0.5*(self.bins[:-1]+self.bins[1:])
        self.rar=[]
        self.decr=[]
        self.rae=[]
        self.dece=[]
        self.rah=[]
        self.dech=[]
        for i in range(self.n):
            print('Facet',i)
            if self.dral[i] is None:
                print('Not fitting, no data')
                self.rar.append([0,0,0,0])
                self.rae.append([100,100,100,100])
                self.decr.append([0,0,0,0])
                self.dece.append([100,100,100,100])
                pass
            else:
                h,_=np.histogram(self.dral[i],self.bins)
                self.rah.append(h)
                p,perr=fitfn(h)
                print('RA Offset is ',p[2],'+/-',perr[2])
                self.rar.append(p)
                self.rae.append(perr)
                h,_=np.histogram(self.ddecl[i],self.bins)
                self.dech.append(h)
                p,perr=fitfn(h)
                self.decr.append(p)
                self.dece.append(perr)
                print('DEC Offset is ',p[2],'+/-',perr[2])
        self.rar=np.array(self.rar)
        self.rae=np.array(self.rae)
        self.decr=np.array(self.decr)
        self.dece=np.array(self.dece)

    def save_fits(self):
        np.save(self.prefix+'-facet_offsets.npy',np.array([self.rar[:,2],self.decr[:,2],self.rae[:,2],self.dece[:,2]]).T)

    def plot_fits(self,pdffile):
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt
        with PdfPages(pdffile) as pdf:
            for i in range(self.n):
                plt.subplot(2,1,1)
                plt.plot(self.bcenter,self.rah[i])
                plt.plot(self.bcenter,model(self.bcenter,*self.rar[i]))
                plt.subplot(2,1,2)
                plt.plot(self.bcenter,self.dech[i])
                plt.plot(self.bcenter,model(self.bcenter,*self.decr[i]))
                plt.suptitle('Facet '+str(i))
                pdf.savefig()
                plt.close()

    def plot_chains(self,pdffile):
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt
        with PdfPages(pdffile) as pdf:
            for j,c in enumerate(self.chains):
                labels=['norm','sigma','offset','bline']
                ndim=len(labels)
                for i in range(len(labels)):
                    plt.subplot(ndim,1,i+1)
                    plt.plot(c[:,:,i].transpose())
                    plt.ylabel(labels[i])
                    plt.xlabel('Samples')
                facet=old_div(j,2)
                chain=j%2
                plt.suptitle('Facet %i chain %i' % (facet,chain))
                pdf.savefig()
                plt.close()

    def plot_offsets(self,lofar_table=None):
        import matplotlib.pyplot as plt
        if lofar_table is not None:
            self.lofar_table=Table.read(lofar_table)
        tf=self.lofar_table

        poly,labels=self.polys,self.labels

        basesize=10
        rarange=(np.min(tf['RA']),np.max(tf['RA']))
        decrange=(np.min(tf['DEC']),np.max(tf['DEC']))
        mdec=np.mean(decrange)
        xstrue=(rarange[1]-rarange[0])*np.cos(mdec*np.pi/180.0)
        ystrue=decrange[1]-decrange[0]
        plt.figure(figsize=(old_div(basesize*xstrue,ystrue), basesize))
        plt.xlim(rarange)
        plt.ylim(decrange)
        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.title('Offsets with method '+self.prefix)

        for p in poly:
            x=[pt[0] for pt in p]
            y=[pt[1] for pt in p]
            plt.plot(x,y,color='black',ls=':')

        mra=[]
        mdec=[]
        mdra=[]
        mddec=[]
        for f in range(self.n):
            t=tf[tf['Facet']==f]
            if len(t)>0:
                mra.append(np.mean(t['RA']))
                mdec.append(np.mean(t['DEC']))
                plt.text(mra[-1],mdec[-1],str(f),color='blue')
                mdra.append(self.rar[f,2])
                mddec.append(self.decr[f,2])
                print(f,len(t),mra[-1],mdec[-1],mdra[-1],mddec[-1])

                plt.gca().invert_xaxis()
                plt.quiver(mra,mdec,mdra,mddec,units = 'xy', angles='xy', scale=1.0,color='red')
                plt.quiver(np.mean(tf['RA']),np.mean(tf['DEC']),1.0,0.0,units = 'xy', angles='xy', scale=1.0,color='green')
                plt.text(np.mean(tf['RA']),np.mean(tf['DEC']),'1 arcsec',color='green')

        plt.savefig('SKO-'+self.prefix+'.png')

    def offsets_to_facetshift(self,filename):

        cellsize=self.cellsize
        outfile=open(filename,'w')
        lines=open('%s.facetCoord.txt'%self.imroot).readlines()
        for l in lines:
            bits=[b.strip() for b in l.split(',')]
            rar=float(bits[2])
            ra=old_div(rar,degtorad)
            decr=float(bits[3])
            dec=old_div(decr,degtorad)
            number=self.r.which_poly(ra,dec)
            #print 'Direction',pli[number]
            direction=self.pli[number]
            print(rar,decr,old_div(-self.rar[direction,2],cellsize),old_div(self.decr[direction,2],cellsize), file=outfile)
        outfile.close()

    def make_astrometry_map(self,outname,factor):
        # factor tells us how much bigger than cellsize the pixels will be
        hdus=fits.open(self.imroot+'.app.restored.fits')
        _,_,yd,xd=hdus[0].data.shape
        xd=int(xd/factor)
        yd=int(yd/factor)
        hdus[0].header['CDELT1']*=factor
        hdus[0].header['CDELT2']*=factor
        hdus[0].header['CRPIX1']/=factor
        hdus[0].header['CRPIX2']/=factor
        w=WCS(hdus[0].header)
        rmap=np.ones((1,1,yd,xd))*np.nan
        # this would be faster with use of e.g. PIL
        for y in range(yd):
            print('.', end=' ')
            sys.stdout.flush()
            xv=np.arange(xd)
            yv=y*np.ones_like(xv)
            ra,dec,_,_=w.wcs_pix2world(xv,yv,0,0,0)
            dra,ddec=self.r.coordconv(ra,dec)[1]
            for i,x in enumerate(xv):
                number=self.r.which_poly(dra[i],ddec[i],convert=False)
                if number is not None:
                    direction=self.pli[number]
                    rmap[0,0,y,x]=np.sqrt(self.rae[direction,2]**2.0+self.dece[direction,2]**2.0)
        print()
        hdus[0].data=rmap
        hdus.writeto(outname,overwrite=True)

    def save(self,filename):
        f = open(filename, 'wb')
        pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    @staticmethod
    def load(filename):
        with open(filename, 'rb') as f:
            return pickle.load(f)

def merge_cat(rootname,rastr='ra',decstr='dec'):
    g=glob.glob(rootname+'/*.vo')
    tlist=[]
    for f in g:
        try:
            t=Table.read(f)
        except:
            print('Error reading table',f)
            raise
        t2=Table()
        t2['ra']=t[rastr]
        t2['dec']=t[decstr]
        tlist.append(t2)

    t=vstack(tlist)
    t2=unique(t,keys=['ra','dec'])
    t2.write(rootname+'.fits',overwrite=True)
    return t2

def do_offsets(o):
    # o is the options file

    if o['mode']!='normal' and  o['mode']!='test':
        raise NotImplementedError('Offsets called with mode '+o['mode'])

    method=o['method']
    image_root='image_full_ampphase_di_m.NS'
    
    report('Determining astrometric offsets with method '+method+' in mode '+o['mode'])
    report('Merging downloaded catalogues')
    if os.path.isfile(method+'.fits'):
        warn('Merged file exists, reading from disk instead')
        data=Table.read(method+'.fits')
    else:
        if method=='pslocal':
            data=Table.read(method+'/'+method+'.txt',format='ascii')
            data['RA'].name='ra'
            data['DEC'].name='dec'
            data.write(method+'.fits')
        else:    
            kwargs={}
            if 'panstarrs' in method:
                kwargs['rastr']='ramean'
                kwargs['decstr']='decmean'
            data=merge_cat(method,**kwargs)

    if o['mode']=='test':
        image_root+='_shift'
        method+='-test'

    report('Running PyBDSM on LOFAR image, please wait...')
    catfile=image_root+'.offset_cat.fits'
    gaulfile=catfile.replace('cat','gaul')
    if os.path.isfile(catfile):
        warn('Catalogue already exists, skipping pybdsf run')
    else:
        if o['mode']=='test':
            suffix='facetRestored'
        else:
            suffix='restored'
        pbimage=image_root+'.int.'+suffix+'.fits'
        nonpbimage=image_root+'.app.'+suffix+'.fits'
        img = bdsm.process_image(pbimage, detection_image=nonpbimage, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0,output_opts=True, output_all=True, atrous_do=False, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None)
        img.write_catalog(outfile=catfile,catalog_type='srl',format='fits',correct_proj='True')
        img.write_catalog(outfile=gaulfile,catalog_type='gaul',format='fits',correct_proj='True')

    lofar=Table.read(catfile)
    print(len(lofar),'LOFAR sources before filtering')
    filter=(lofar['E_RA']*3600.0)<2.0
    filter&=(lofar['E_DEC']*3600.0)<2.0
    filter&=(lofar['Maj']*3600.0)<10
    lofar=lofar[filter]
    print(len(lofar),'LOFAR sources after filtering')
    regfile=image_root+'.tessel.reg'
    cra,cdec=get_centpos()
    report('Set up structure')
    if o['clusterfile'] is None:
        clusterfile="image_dirin_SSD_m.npy.ClusterCat.npy"
    else:
        clusterfile=o['clusterfile']
    NDir=np.load(clusterfile).shape[0]
    oo=Offsets(method,n=NDir,imroot=image_root,cellsize=o['cellsize'],fitmethod=o['fit'])
    report('Label table')
    lofar_l=oo.r.add_facet_labels(lofar)
    report('Finding offsets')
    oo.find_offsets(lofar_l,data)
    report('Fitting offsets')
    oo.fit_offsets()
    report('Making plots and saving output')
    #oo.plot_fits(method+'-fits.pdf')
    oo.save_fits()
    oo.plot_offsets()
    if 'test' not in o['mode']:
        oo.save(method+'-fit_state.pickle')
        report('Making astrometry error map, please wait')
        oo.make_astrometry_map('astromap.fits',20)
        oo.offsets_to_facetshift('facet-offset.txt')

if __name__=='__main__':
    from options import options
    from parset import option_list
    o=options(sys.argv[1:],option_list)
    do_offsets(o)
