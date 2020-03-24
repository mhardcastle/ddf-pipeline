#!/usr/bin/env python

from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import optparse
import pickle
import numpy as np
import os
from itertools import product as ItP
#from DDFacet.Other import MyLogger
#log=MyLogger.getLogger("ClassSmooth")

SaveName="last_Smooth.obj"

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = list(range(order+1))
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

    
def read_options():
    desc="""Questions and suggestions: cyril.tasse@obspm.fr"""
    global options
    opt = optparse.OptionParser(usage='Usage: %prog --ms=somename.MS <options>',version='%prog version 1.0',description=desc)

    data = optparse.OptionGroup(opt, "* Data-related options", "Won't work if not specified.")
    data.add_option('--MSName',help='MSName [no default]',default='')
    data.add_option('--SolsFile',help='Solfile [no default]',default='')
    opt.add_option_group(data)
    group = optparse.OptionGroup(opt, "* Smoothing options", "Defaults may be sensible")
    group.add_option('--WSize',help='Smooth window size',default=53)
    group.add_option('--Order',help='Smooth order',default=2)
    opt.add_option_group(group)
    group = optparse.OptionGroup(opt, "* Misc options", "Defaults may be sensible")
    group.add_option('--Plot',help='Enable plotting',default=False)
    opt.add_option_group(group)


    options, arguments = opt.parse_args()
    f = open(SaveName,"wb")
    pickle.dump(options,f)


def NormMatrices(G):
    nt,nch,na,_,_=G.shape

    for iChan,it in ItP(list(range(nch)),list(range(nt))):
        Gt=G[it,iChan,:,:]
        u,s,v=np.linalg.svd(Gt[0])
        # #J0/=np.linalg.det(J0)
        # J0=Gt[0]
        # JJ=np.dot(J0.T.conj(),J0)
        # sqJJ=ModLinAlg.sqrtSVD(JJ)
        # sqJJinv=ModLinAlg.invSVD(JJ)
        # U=np.dot(J0,sqJJinv)
        U=np.dot(u,v)
        for iAnt in range(0,na):
            Gt[iAnt,:,:]=np.dot(U.T.conj(),Gt[iAnt,:,:])
            #Gt[iAnt,:,:]=np.dot(np.dot(u,Gt[iAnt,:,:]),v.T.conj())
            #Gt[iAnt,:,:]=np.dot(Gt[iAnt,:,:],J0)
    return G


class ClassSmooth(object):
    def __init__(self,MSName,SolsName,Type="linear",WSize=53,Order=2,PolMode="Full",OutName=None,doplot=False):

        SolsName="killMS.%s.sols.npz"%SolsName
        self.FileName="/".join([os.path.abspath(MSName),SolsName])
        self.OutName=OutName
        print("Smoothing from %s"%self.FileName)
        self.DicoFile=dict(np.load(self.FileName))
        self.Sols=self.DicoFile["Sols"]
        self.Sols=self.Sols.view(np.recarray)
        self.WSize=WSize
        self.Order=Order
        self.PolMode=PolMode
        self.doplot=doplot
        self.NormAllDirs()

    def NormAllDirs(self):
        print("  Normalising Jones matrices ....")
        nt,nch,na,nd,_,_=self.Sols.G.shape
        for iDir in range(nd):
            G=self.Sols.G[:,:,:,iDir,:,:]
            self.Sols.G[:,:,:,iDir,:,:]=NormMatrices(G)

    def Smooth(self):
        Sols0=self.Sols
        nt0,nch,na,nd,_,_=Sols0.G.shape
        G0=Sols0.G.reshape((nt0,nch,na,nd,4))

        Sols1=np.zeros(nt0,dtype=Sols0.dtype)
        Sols1=Sols1.view(np.recarray)
        #nt1,nch,na,nd,_,_=Sols1.G.shape

        G1=Sols1.G.reshape((nt0,na,nd,4))
        if self.PolMode=="Full":
            Pols=list(range(4))

        Sols1.t0=Sols0.t0
        Sols1.t1=Sols0.t1
#        Sols1.tm=Sols0.tm

        print("  Smoothing")
        for iDir in range(nd):
            for iAnt in range(na):
                for ipol in Pols:
                    # Amplitude
                    yp=np.abs(G0[:,0,iAnt,iDir,ipol])
                    # Do some smoothing on the amps
                    yp = np.copy(savitzky_golay(yp, self.WSize, self.Order))
                    G1[:,iAnt,iDir,ipol]=yp[:]
                    # Phase
                    yp=np.angle(G0[:,0,iAnt,iDir,ipol])
                    G1[:,iAnt,iDir,ipol]*=np.exp(1j*yp[:])
                if self.doplot:
                    import matplotlib.pyplot as plt
                    xp=(Sols0.t0+Sols0.t1)/2.
                    op0=np.abs
                    op1=np.angle
                    plt.clf()
                    plt.suptitle('Direction = %i Antenna = %i' % (iDir,iAnt))
                    for ipol in Pols:
                        plt.subplot(2,1,1)
                        plt.scatter(xp,op0(G0[:,0,iAnt,iDir,ipol]))
                        plt.plot(xp,op0(G1[:,iAnt,iDir,ipol]),marker=".",ls="",label=str(ipol))
                        plt.legend(loc=0)

                        plt.subplot(2,1,2)
                        plt.scatter(xp,op1(G0[:,0,iAnt,iDir,ipol]))
                        plt.plot(xp,op1(G1[:,iAnt,iDir,ipol]))

                    plt.draw()
                    plt.show(False)
                    plt.pause(0.1)

        G1=G1.reshape((nt0,nch,na,nd,2,2))
        Sols1.G=G1
        self.Sols1=Sols1
     
    def Save(self):
        self.DicoFile["Sols"]=self.Sols1
        OutName=self.OutName
        if OutName==None:
            FileName=self.FileName.split("/")[-1]
            Path="/".join(self.FileName.split("/")[0:-1])+"/"
            Name=".".join(FileName.split(".")[1:-2])
            OutName="%skillMS.%s.Smooth.sols.npz"%(Path,Name)
        print("  Saving smoothed solutions in: %s"%OutName)
        np.savez(OutName,**self.DicoFile)

        
def test():
    FileName="killMS.KAFCA.sols.npz"
    CI=ClassSmooth(FileName)
    CI.Smooth()
    CI.Save()

def main(options=None):
    if options is None:
        f = open(SaveName,'rb')
        options = pickle.load(f)
    #FileName="killMS.KAFCA.sols.npz"

    SolsFile=options.SolsFile

    MSName=options.MSName
    if ".txt" in MSName:
        f=open(MSName)
        Ls=f.readlines()
        f.close()
        MSName=[]
        for l in Ls:
            ll=l.replace("\n","")
            MSName.append(ll)
        lMS=MSName
        print("In batch mode, running Smooth on the following MS:")
        for MS in lMS:
            print("  %s"%MS)
    elif "*" in options.MSName:
        Patern=options.MSName
        lMS=sorted(glob.glob(Patern))
        print("In batch mode, running Smooth on the following MS:")
        for MS in lMS:
            print("  %s"%MS)
    else:
        lMS=[options.MSName]



    for MSName in lMS:
        CI=ClassSmooth(MSName,SolsFile,WSize=options.WSize,Order=options.Order,doplot=options.Plot)
        CI.Smooth()
        CI.Save()


if __name__=="__main__":
    read_options()
    f = open(SaveName,'rb')
    options = pickle.load(f)


    main(options=options)
