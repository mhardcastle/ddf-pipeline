import numpy as np
from scipy.optimize import leastsq
import scipy
import pyfits
import os,sys
from pipeline_logging import run_log
from subprocess import call

class bcolors:
   HEADER = '\033[95m'
   OKBLUE = '\033[94m'
   OKGREEN = '\033[92m'
   WARNING = '\033[93m'
   FAIL = '\033[91m'
   ENDC = '\033[0m'

def die(s):
   print bcolors.FAIL+s+bcolors.ENDC
   raise Exception(s)

def report(s):
   print bcolors.OKGREEN+s+bcolors.ENDC

def warn(s):
   print bcolors.OKBLUE+s+bcolors.ENDC

def run(s,proceed=False,dryrun=False,log=None,quiet=False):
   report('Running: '+s)
   if not dryrun:
      if log is None:
         #print 'Log is none, using "call" to run '+s
         retval=call(s,shell=True)
      else:
         retval=run_log(s,log,quiet)
      if not(proceed) and retval!=0:
         os.system('CleanSHM.py')
         die('FAILED to run '+s+': return value is '+str(retval))

      return retval
   else:
      warn('Dry run, skipping this step')

def close_higher_bsmooth_number(n,B):
    """Finds the closest higher B-smooth number (good numbers for FFTs) http://en.wikipedia.org/wi$
        """
    fivesmooths = np.array([1])
    i = np.max([n-300,5])
    while np.max(fivesmooths) < n:
        if prime_factors(i) <=B and (i/2.0 == int(i)/2): #Image must have an even number of pixels
            fivesmooths = np.append(fivesmooths,i)
        i+=1
    print fivesmooths
    fivesmooths=np.sort(fivesmooths)
    closest_smooth_index = np.where(abs(fivesmooths-n)==np.min(abs(fivesmooths-n)))
    if len(fivesmooths[closest_smooth_index[0]]) == 2:
	closest_smooth = fivesmooths[closest_smooth_index[0][1]]
        return closest_smooth
    print closest_smooth_index
    if fivesmooths[closest_smooth_index[0]] < n:
        closest_smooth = fivesmooths[closest_smooth_index[0]+1]
    else:
        closest_smooth = fivesmooths[closest_smooth_index[0]]
    return closest_smooth[0]


#-----------------------------------------------------------

def find_imagenoise(workingimage,estnoise):
    f = pyfits.open('%s'%workingimage)
    noisearray = f[0].data.flatten()
    maxpixel = np.max(noisearray)
    noisearray = np.random.permutation(noisearray)[:10000]
    noisepix = np.array(filter(lambda x: abs(x) > 10E-8,noisearray)) # Filter out the edge pixels which have values around 1E-10
    noisepix = np.array(filter(lambda x: abs(x)<50.0*estnoise,noisepix))
    f.close()
    rms = fit_gaussian_histogram(noisepix,'n')
    return rms

#-----------------------------------------------------------

def prime_factors(n):
    """Find the prime factors of a number
    """
    factors=[]
    d=2
    while(d*d<n):
        while(n>1):            
            while n%d==0:
                factors.append(d)
                n=n/d
            d+=1
    return factors[-1]

#-----------------------------------------------------------
    
def close_higher_bsmooth_number(n,B):
    """Finds the closest higher B-smooth number (good numbers for FFTs) http://en.wikipedia.org/wi$
        """
    fivesmooths = np.array([1])
    i = np.max([n-300,5])
    while np.max(fivesmooths) < n:
        if prime_factors(i) <=B and (i/2.0 == int(i)/2): #Image must have an even number of pixels
            fivesmooths = np.append(fivesmooths,i)
        i+=1
    print fivesmooths
    fivesmooths=np.sort(fivesmooths)
    closest_smooth_index = np.where(abs(fivesmooths-n)==np.min(abs(fivesmooths-n)))
    if len(fivesmooths[closest_smooth_index[0]]) == 2:
	closest_smooth = fivesmooths[closest_smooth_index[0][1]]
        return closest_smooth
    print closest_smooth_index
    if fivesmooths[closest_smooth_index[0]] < n:
        closest_smooth = fivesmooths[closest_smooth_index[0]+1]
    else:
        closest_smooth = fivesmooths[closest_smooth_index[0]]
    return closest_smooth[0]

#-----------------------------------------------------------
#------------------------------------------------------------

def model_gaussian(t, coeffs):
    return coeffs[0] + coeffs[1] * np.exp( - ((t-coeffs[2])**2.0/(2*coeffs[3]**2)))

#------------------------------------------------------------
def residuals_gaussian(coeffs, y, t):
    return y - model_gaussian(t, coeffs)

#------------------------------------------------------------
def fit_gaussian_histogram(pixelvals,plotting):

    fitnumbers,cellsizes = np.histogram(pixelvals,100)
    sigmaguess = np.std(pixelvals)/(abs(cellsizes[1]-cellsizes[0]))
    x0 = [0.0,max(fitnumbers),np.where(fitnumbers==max(fitnumbers))[0][0],sigmaguess] #Offset amp, amp, x-offset, sigma
    t = np.arange(len(fitnumbers))
    x, flag = scipy.optimize.leastsq(residuals_gaussian, x0, args=(fitnumbers, t))

    if plotting == 'y':
        pylab.plot(fitnumbers)
        pylab.plot(t,fitnumbers,t,model_gaussian(t,x))
        pylab.show()
        pylab.close()
        pylab.cla()

    #print 'Sigma is %s'%(x[3]*abs(cellsizes[1]-cellsizes[0]))
    
    return (x[3]*abs(cellsizes[1]-cellsizes[0]))
