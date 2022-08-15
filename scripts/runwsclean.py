#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import map
from builtins import str
from builtins import range
from past.utils import old_div
import logging
logging.basicConfig(filename='selfcal.log', format='%(levelname)s:%(asctime)s ---- %(message)s', datefmt='%m/%d/%Y %I:%M:%S', level=logging.DEBUG)
from auxcodes import flatten
import matplotlib
matplotlib.use('Agg')
import os, sys
import numpy as np
import losoto
import losoto.lib_operations
import glob
from astropy.io import fits
import pyrap.tables as pt
import os.path
from losoto import h5parm
import bdsf
import pyregion
import argparse
import pickle
import aplpy
from lofar.stationresponse import stationresponse
try:
  from getcpus import getcpus
  getcpuworks = True
except:
  getcpuworks = False 

def removenans(parmdb, soltab):
   H5 = h5parm.h5parm(parmdb, readonly=False)
   vals =H5.getSolset('sol000').getSoltab(soltab).getValues()[0]
   weights = H5.getSolset('sol000').getSoltab(soltab).getValues(weight=True)[0]   
   
   idxnan  = np.where((~np.isfinite(vals))) 
   
   print(idxnan)
   print('Found some NaNs', vals[idxnan])

   if H5.getSolset('sol000').getSoltab(soltab).getType() == 'phase':
       vals[idxnan] = 0.0
   if H5.getSolset('sol000').getSoltab(soltab).getType() == 'amplitude':
       vals[idxnan] = 1.0            
   
   weights[idxnan] = 0.0

   H5.getSolset('sol000').getSoltab(soltab).setValues(weights,weight=True)
   H5.getSolset('sol000').getSoltab(soltab).setValues(vals)
   H5.close()
   return


def losotolofarbeam(parmdb, soltabname, ms, inverse=False, useElementResponse=True, useArrayFactor=True, useChanFreq=True):
    """
    Do the beam correction via this imported losoto operation
    """

    H5 = h5parm.h5parm(parmdb, readonly=False)
    soltab = H5.getSolset('sol000').getSoltab(soltabname)

    #t = pt.table(ms)
    sr = stationresponse(ms, inverse, useElementResponse, useArrayFactor, useChanFreq)

    numants = pt.taql('select gcount(*) as numants from '+ms+'::ANTENNA').getcol('numants')[0]
    times = soltab.getAxisValues('time')

    for vals, coord, selection in soltab.getValuesIter(returnAxes=['ant','time','pol','freq'], weight=False):
        vals = losoto.lib_operations.reorderAxes( vals, soltab.getAxesNames(), ['ant','time','freq','pol'] )

        for stationnum in range(numants):
            logging.debug('Working on station number %i' % stationnum)
            for itime, time in enumerate(times):
                beam = sr.evaluateStation(time=time, station=stationnum)
                # Reshape from [nfreq, 2, 2] to [nfreq, 4]
                beam = beam.reshape(beam.shape[0], 4)

                if soltab.getAxisLen('pol') == 2:
                    beam = beam[:,[0,3]] # get only XX and YY

                if soltab.getType() == 'amplitude':
                    vals[stationnum, itime, :, :] = np.abs(beam)
                elif soltab.getType() == 'phase':
                    vals[stationnum, itime, :, :] = np.angle(beam)
                else:
                    logging.error('Beam prediction works only for amplitude/phase solution tables.')
                    return 1

        vals = losoto.lib_operations.reorderAxes( vals, ['ant','time','freq','pol'], [ax for ax in soltab.getAxesNames() if ax in ['ant','time','freq','pol']] )
        soltab.setValues(vals, selection)
    
    H5.close()
    return
 

#losotolofarbeam('P214+55_PSZ2G098.44+56.59.dysco.sub.shift.avg.weights.ms.archive_templatejones.h5', 'amplitude000', 'P214+55_PSZ2G098.44+56.59.dysco.sub.shift.avg.weights.ms.archive', inverse=False, useElementResponse=False, useArrayFactor=True, useChanFreq=True)

#sys.exit()

def cleanup(mslist):
    
    for ms in mslist:
        os.system('rm -rf ' + ms)
        
    os.system('rm -f *first-residual.fits')    
    os.system('rm -f *psf.fits') 
    os.system('rm -f *-00*-*.fits')
    os.system('rm -f *dirty.fits')
    os.system('rm -f solintimage*model.fits')
    os.system('rm -f solintimage*residual.fits')
    os.system('rm -f *pybdsf.log')
    return


def flagms_startend(ms, tecsolsfile, tecsolint):
    

    taql = 'taql'
       
    msout = ms + '.cut'
    
    H5 =h5parm.h5parm(tecsolsfile)
    tec = H5.getSolset('sol000').getSoltab('tec000').getValues() 
    tecvals = tec[0]
    
    axis_names = H5.getSolset('sol000').getSoltab('tec000').getAxesNames()
    time_ind = axis_names.index('time')
    ant_ind = axis_names.index('ant')
    
    #['time', 'ant', 'dir', 'freq']
    reftec = tecvals[:,0,0,0] 
    
    #print np.shape( tecvals[:,:,0,0]), np.shape( reftec[:,None]), 
    tecvals = tecvals[:,:,0,0] - reftec[:,None] # reference to zero
    
    times   = tec[1]['time']
    
    #print tecvals[:,0]
    
    goodtimesvec = []
    
    for timeid, time in enumerate(times):
    
      tecvals[timeid,:]
      
      #print timeid, np.count_nonzero( tecvals[timeid,:])
      goodtimesvec.append(np.count_nonzero( tecvals[timeid,:]))



    goodstartid = np.argmax (np.array(goodtimesvec) > 0)
    goodendid   = len(goodtimesvec) - np.argmax (np.array(goodtimesvec[::-1]) > 0)
    
    print('First good solutionslot,', goodstartid, ' out of', len(goodtimesvec))
    print('Last good solutionslot,', goodendid, ' out of', len(goodtimesvec))    
    H5.close()
    
    cmd = taql + " ' select from " + ms + " where TIME in (select distinct TIME from " + ms 
    cmd+= " offset " + str(goodstartid*np.int(tecsolint)) 
    cmd+= " limit " + str((goodendid-goodstartid)*np.int(tecsolint)) +") giving " 
    cmd+= msout + " as plain'"
    
    print(cmd)
    os.system(cmd)
    
    os.system('rm -rf ' + ms)
    os.system('mv ' + msout + ' ' + ms)
    return


#flagms_startend('P215+50_PSZ2G089.52+62.34.dysco.sub.shift.avg.weights.ms.archive','phaseonlyP215+50_PSZ2G089.52+62.34.dysco.sub.shift.avg.weights.ms.archivesolsgrid_9.h5', 2)
#sys.exit()

def removestartendms(ms, starttime=None, endtime=None):

    # chdeck if output is already there and remove    
    if os.path.isdir(ms + '.cut'):
          os.system('rm -rf ' + ms + '.cut')  
    if os.path.isdir(ms + '.cuttmp'):
          os.system('rm -rf ' + ms + '.cuttmp')  

        
    cmd = 'DP3 msin=' + ms + ' ' + 'msout.storagemanager=dysco msout=' + ms + '.cut '
    cmd+=  'msin.weightcolumn=WEIGHT_SPECTRUM steps=[] ' 
    if starttime is not None:
      cmd+= 'msin.starttime=' + starttime + ' '
    if endtime is not None:  
      cmd+= 'msin.endtime=' + endtime   + ' '   
    print(cmd)  
    os.system(cmd)
    
    cmd = 'DP3 msin=' + ms + ' ' + 'msout.storagemanager=dysco msout=' + ms + '.cuttmp '
    cmd+= 'msin.weightcolumn=WEIGHT_SPECTRUM_SOLVE steps=[] '  
    if starttime is not None:
      cmd+= 'msin.starttime=' + starttime + ' '
    if endtime is not None:  
      cmd+= 'msin.endtime=' + endtime   + ' '
    print(cmd)
    os.system(cmd)    


    # Make a WEIGHT_SPECTRUM from WEIGHT_SPECTRUM_SOLVE
    t  = pt.table(ms + '.cut' , readonly=False)

    print('Adding WEIGHT_SPECTRUM_SOLVE') 
    desc = t.getcoldesc('WEIGHT_SPECTRUM')
    desc['name']='WEIGHT_SPECTRUM_SOLVE'
    t.addcols(desc)

    t2 = pt.table(ms + '.cuttmp' , readonly=True)
    imweights = t2.getcol('WEIGHT_SPECTRUM')
    t.putcol('WEIGHT_SPECTRUM_SOLVE', imweights)

    # Fill WEIGHT_SPECTRUM with WEIGHT_SPECTRUM from second ms
    t2.close()
    t.close() 

    # clean up
    os.system('rm -rf ' + ms + '.cuttmp')



  
    return




def which(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None



def plotimage(fitsimagename, outplotname, mask=None, rmsnoiseimage=None):

   #logging.basicConfig(level=logging.ERROR)   # to block astropy/aplpy warnings about fits headers
   #image noise for plotting
   if rmsnoiseimage == None:
      hdulist = fits.open(fitsimagename)
   else:
      hdulist = fits.open(rmsnoiseimage)   
   imagenoise = findrms(np.ndarray.flatten(hdulist[0].data))
   hdulist.close() 
  
   #image noise info
   hdulist = fits.open(fitsimagename) 
   imagenoiseinfo = findrms(np.ndarray.flatten(hdulist[0].data))
   ffits = flatten(hdulist)
  
   f = aplpy.FITSFigure(ffits)
   f.show_colorscale(vmax=16*imagenoise, vmin=-6*imagenoise, cmap='bone')
   f.set_title(fitsimagename+' (noise = {} mJy/beam)'.format(round(imagenoiseinfo*1e3, 3)))
   f.add_beam()
   f.beam.set_frame(True)
   f.beam.set_color('white')
   f.beam.set_edgecolor('black')
   f.beam.set_linewidth(1.)
   f.add_grid()
   f.grid.set_color('white')
   f.grid.set_alpha(0.5)
   f.grid.set_linewidth(0.2)
   f.add_colorbar()
   f.colorbar.set_axis_label_text('Flux (Jy beam$^{-1}$)')
   if mask is not None:
      maskhdu=fits.open(mask)
      fmask=flatten(maskhdu)
      f.show_contour(fmask, colors='red', levels=[0.1*imagenoise], filled=False, smooth=1, alpha=0.6, linewidths=1)
   f.save(outplotname, dpi=120, format='png')
   #logging.basicConfig(level=logging.DEBUG)
   logging.info(fitsimagename + ' RMS noise: ' + str(imagenoiseinfo))
   hdulist.close()   

def archive(mslist, outtarname, regionfile, fitsmask, imagename, ncpu=None):
  for ms in mslist:
    msout = ms + '.calibrated'
    if os.path.isdir(msout):
      os.system('rm -rf ' + msout)
    cmd  = 'DP3 '
    if ncpu is not None:
       cmd +='numthreads=%i' % ncpu 
    cmd += 'msin=' + ms + ' msout=' + msout + ' '
    cmd +='msin.datacolumn=CORRECTED_DATA msout.storagemanager=dysco steps=[]'
    os.system(cmd)

  msliststring = ' '.join(map(str, glob.glob('*.calibrated') ))
  cmd = 'tar -zcf ' + outtarname + ' ' + msliststring + ' selfcal.log ' + regionfile + ' ' + \
        fitsmask + ' ' + imagename + ' '
  if os.path.isfile(outtarname):
      os.system('rm -f ' + outtarname)
  logging.info('Creating archived calibrated tarball: ' + outtarname)    
  os.system(cmd)
  
  for ms in mslist:
    msout = ms + '.calibrated'   
    os.system('rm -rf ' + msout)
  return


def reweight(mslist, pixsize, imsize, channelsout, niter, robust, multiscale=False, fitsmask=None):
   """
   determine the solution time and frequency intervals based on the amount of compact source flux
   """
   
   rmslist = []

   logging.info('Adjusting weights')

   for ms in mslist:
          imageout =  'rmsimage' + ms.split('.ms')[0] 
          makeimage([ms], imageout, pixsize, imsize, channelsout, np.int(old_div(niter,(len(mslist)**(1./3.)))), robust, multiscale=multiscale, predict=False,fitsmask=fitsmask)
          
          hdulist = fits.open(imageout + '-MFS-image.fits')
          imagenoise = findrms(np.ndarray.flatten(hdulist[0].data))
          hdulist.close() 
          rmslist.append(imagenoise)
          
   weightslist = []       
   return 

def determinesolints(mslist, pixsize, imsize, channelsout, niter, robust, TEC, multiscale=False,ncpu=None):
   """
   determine the solution time and frequency intervals based on the amount of compact source flux
   """
   if os.path.isfile('nchan_phase.p') and os.path.isfile('solint_phase.p') and \
      os.path.isfile('solint_ap.p') and os.path.isfile('nchan_ap.p'):
    
      with open('nchan_phase.p', 'rb') as f:
         nchan_phase_F = pickle.load(f)        
  
      with open('solint_phase.p', 'rb') as f:
         solint_phase_F = pickle.load(f)        

      with open('solint_ap.p', 'rb') as f:
         solint_ap_F = pickle.load(f)
  
      with open('nchan_ap.p', 'rb') as f:
         nchan_ap_F = pickle.load(f)        
  
   else:
      decl  = getdeclinationms(mslist[0])
      declf = declination_sensivity_factor(decl)
      
      logging.info('Noise inceases by a factor of: ' + str(declf) + ' for this observations at delclination [deg] ' + str(decl))

      nchan_phase_F  = []
      solint_phase_F = []
      solint_ap_F    = []
      nchan_ap_F     = [] 

      for ms in mslist:
          imageout =  'solintimage' + ms.split('.ms')[0] 
          makeimage([ms], imageout, pixsize, imsize, channelsout, np.int(old_div(niter,2)), robust, multiscale=multiscale, predict=False, ncpu=ncpu)
          csf = determine_compactsource_flux(imageout + '-MFS-image.fits') 
        
          nchan_phase, solint_phase, solint_ap = calculate_solintnchan(old_div(csf,declf) )
          solint_phase = 2*solint_phase
          nchan_ap = nchan_phase # keep the same
          logging.info('MS, NCHAN, SOLINT_PHASE, SOLINT_AP: ' + str(ms) + ' ' + str(nchan_phase) + ' ' + str (solint_phase) + ' ' + str(solint_ap))
          if TEC:
            nchan_phase  = 5.
        
          nchan_phase_F.append(nchan_phase)
          solint_phase_F.append(solint_phase)
          solint_ap_F.append(solint_ap)
          nchan_ap_F.append(nchan_ap)

      with open('nchan_phase.p', 'wb') as f:
         pickle.dump(nchan_phase_F,f)
  
      with open('solint_phase.p', 'wb') as f:
         pickle.dump(solint_phase_F,f)

      with open('solint_ap.p', 'wb') as f:
         pickle.dump(solint_ap_F,f)
  
      with open('nchan_ap.p', 'wb') as f:
         pickle.dump(nchan_ap_F,f)

   return nchan_phase_F, solint_phase_F, solint_ap_F, nchan_ap_F

def create_beamcortemplate(ms,ncpu=None):
  """
  create a DP3 gain H5 template solutution file that can be filled with losoto
  """
  H5name = ms + '_templatejones.h5'   

  cmd = 'DP3 '
  if ncpu is not None:
     cmd +='numthreads=%i ' % ncpu
  cmd += 'msin=' + ms + ' msin.datacolumn=DATA msout=. '
  cmd += 'msin.modelcolumn=DATA '
  cmd += 'steps=[ddecal] ddecal.type=ddecal '
  cmd += 'ddecal.maxiter=1 ddecal.usemodelcolumn=True ddecal.nchan=1 '
  cmd += 'ddecal.mode=complexgain ddecal.h5parm=' + H5name  + ' '
  cmd += 'ddecal.solint=10'

  print(cmd)
  os.system(cmd)

  return H5name

def create_losoto_beamcorparset(ms):
    """
    Create a losoto parset to fill the beam correction values'.
    """
    parset = 'losotobeam.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')

    f.write('pol = [XX,YY]\n')
    f.write('soltab = [sol000/*]\n\n\n')
    
    #f.write('[beam]\n')
    #f.write('ms = %s\n' % ms)
    #f.write('operation = LOFARBEAM\n')
    #f.write('useElementResponse = False\n')
    #f.write('useArrayFactor = True\n')
    #f.write('useChanFreq = True\n\n\n')

    f.write('[plotphase]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-0.5,0.5]\n')
    f.write('prefix = plotlosoto%s/phases_beam\n' % ms)
    f.write('refAnt = CS001HBA0\n\n\n')

    f.write('[plotamp]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/amplitude000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [0.2,1]\n')
    f.write('prefix = plotlosoto%s/amplitudes_beam\n' %ms)

    f.close()
    return parset

def create_losoto_tecandphaseparset(ms):
    parset = 'losoto_plotfasttecandphase.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')
  
    f.write('pol = [XX, YY]\n')
    f.write('Ncpu = 0\n\n\n')

    f.write('[plottecandphase]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('axesInPlot = [time]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-3.14,3.14]\n')
    f.write('soltabsToAdd = tec000\n')
    f.write('prefix = plotlosoto%s/fasttecandphase\n' % ms)
    f.write('refAnt = CS003HBA0\n')
  
    f.close()
    return parset


def create_losoto_fastphaseparset(ms):
    parset = 'losoto_plotfastphase.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')

    f.write('pol = [XX, YY]\n')
    f.write('soltab = [sol000/*]\n')
    f.write('Ncpu = 0\n\n\n')

    f.write('[plotphase]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-3.14,3.14]\n')
    f.write('prefix = plotlosoto%s/fastphase\n' % ms)
    f.write('refAnt = CS001HBA0\n')

    f.close()
    return parset


def create_losoto_flag_apgridparset(ms):

    parset= 'losoto_flag_apgrid.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')

    f.write('pol = [XX, YY]\n')
    f.write('soltab = [sol000/*]\n')
    f.write('Ncpu = 0\n\n\n')
   
    f.write('[plotamp]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/amplitude000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [0,2.5]\n')
    f.write('prefix = plotlosoto%s/slowamp\n\n\n' % ms)

    f.write('[plotphase]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-3.14,3.14]\n')
    f.write('prefix = plotlosoto%s/slowphase\n' % ms)
    f.write('refAnt = CS003HBA0\n\n\n')

    f.write('[flagamp]\n')
    f.write('soltab = [sol000/amplitude000]\n')
    f.write('operation = FLAG\n')
    f.write('axesToFlag = [time,freq]\n')
    f.write('mode = smooth\n')
    f.write('maxCycles = 3\n')
    f.write('windowNoise = 7\n')
    f.write('maxRms = 7.\n')
    f.write('order  = [5,5]\n\n\n')
  
    f.write('[flagphase]\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('operation = FLAG\n')
    f.write('axesToFlag = [time,freq]\n')
    f.write('mode = smooth\n')
    f.write('maxCycles = 3\n')
    f.write('windowNoise = 7\n')
    f.write('maxRms = 7.\n')
    f.write('order  = [5,5]\n\n\n')

    f.write('[plotampafter]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/amplitude000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [0,2.5]\n')
    f.write('prefix = plotlosoto%s/ampsm\n\n\n' % ms)

    f.write('[plotphase_after]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-3.14,3.14]\n')
    f.write('prefix = plotlosoto%s/phasesm\n' % ms)
    f.write('refAnt = CS003HBA0\n')
  
    f.close()
    return parset

def create_losoto_mediumsmoothparset(ms, boxsize):
    parset= 'losoto_mediansmooth.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')

    f.write('pol = [XX, YY]\n')
    f.write('soltab = [sol000/*]\n')
    f.write('Ncpu = 0\n\n\n')

    f.write('[smoothphase]\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('operation= SMOOTH\n')
    f.write('axesToSmooth = [freq,time]\n')
    f.write('size = [%s,%s]\n' % (boxsize, boxsize))
    f.write('mode = runningmedian\n\n\n')

    f.write('[smoothamp]\n')
    f.write('soltab = [sol000/amplitude000]\n')
    f.write('operation= SMOOTH\n')
    f.write('axesToSmooth = [freq,time]\n')
    f.write('size = [%s,%s]\n' % (boxsize, boxsize))
    f.write('mode = runningmedian\n\n\n')

    f.write('[plotamp_after]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/amplitude000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [0,2.5]\n')
    f.write('prefix = plotlosoto%s/amps_smoothed\n\n\n' % ms)

    f.write('[plotphase_after]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-3.14,3.14]\n')
    f.write('prefix = plotlosoto%s/phases_smoothed\n\n\n' % ms)
    f.write('refAnt = CS003HBA0\n')


    f.write('[plotphase_after1rad]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-1,1]\n')
    f.write('prefix = plotlosoto%s/phases_smoothed1rad\n' % ms)
    f.write('refAnt = CS003HBA0\n')

    f.close()
    return parset


def beamcor(ms,ncpu=None):
    """
    correct a ms for the beam in the phase center (array_factor only)
    """
    H5name = create_beamcortemplate(ms)
    parset = create_losoto_beamcorparset(ms)

    losotolofarbeam(H5name, 'phase000', ms, useElementResponse=False, useArrayFactor=True, useChanFreq=True)
    losotolofarbeam(H5name, 'amplitude000', ms, useElementResponse=False, useArrayFactor=True, useChanFreq=True)   

    losoto = 'losoto'
    taql = 'taql'


    
    cmdlosoto = losoto + ' ' + H5name + ' ' + parset
    print(cmdlosoto)
    os.system(cmdlosoto)
    
    cmd = 'DP3 '
    if ncpu is not None:
       cmd +='numthreads=%i ' % ncpu
    cmd += 'msin=' + ms + ' msin.datacolumn=DATA msout=. '
    cmd += 'msin.weightcolumn=WEIGHT_SPECTRUM '
    cmd += 'msout.datacolumn=CORRECTED_DATA steps=[ac1,ac2] msout.storagemanager=dysco '
    cmd += 'ac1.parmdb='+H5name + ' ac2.parmdb='+H5name + ' '
    cmd += 'ac1.type=applycal ac2.type=applycal '
    cmd += 'ac1.correction=phase000 ac2.correction=amplitude000 ac2.updateweights=True ' 
    print('DP3 applycal:', cmd)
    os.system(cmd)
   
    
    os.system(taql + " 'update " + ms + " set DATA=CORRECTED_DATA'")
    return


def findrms(mIn,maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    m=mIn[np.abs(mIn)>maskSup]
    rmsold=np.std(m)
    diff=1e-1
    cut=3.
    bins=np.arange(np.min(m),np.max(m),(np.max(m)-np.min(m))/30.)
    med=np.median(m)
    for i in range(10):
        ind=np.where(np.abs(m-med)<rmsold*cut)[0]
        rms=np.std(m[ind])
        if np.abs(old_div((rms-rmsold),rmsold))<diff: break
        rmsold=rms
    return rms


def findamplitudenoise(parmdb):
      """
      find the 'amplitude noise' in a parmdb, return non-clipped rms value
      """
      H5 = h5parm.h5parm(parmdb, readonly=True) 
      amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
      weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
      H5.close()

      idx = np.where(weights != 0.0)
      
      amps = amps[idx]
      amps = amps[np.isfinite(amps)]
      amps = np.log10(np.ndarray.flatten(amps))
      
      
      noise = findrms(amps)
      
      logging.info('Noise and clipped noise' + str(parmdb) + ' ' + str(np.std(amps)) + ' ' + str(noise))

      # do not return clipped noise, we are intersted in finding data with high outliers
      return np.std(amps)


def getimsize(boxfile, cellsize=1.5):
   """
   find imsize need to image a DS9 boxfile region
   """
   r = pyregion.open(boxfile)
   
   xs = np.ceil(old_div((r[0].coord_list[2])*1.6*3600.,cellsize))
   ys = np.ceil(old_div((r[0].coord_list[3])*1.6*3600.,cellsize))

   imsize = np.ceil(xs) # // Round up decimals to an integer
   if(imsize % 2 == 1): 
       imsize = imsize + 1
   return np.int(imsize)


def smoothsols(parmdb, ms):
    

    losoto = 'losoto'    
    
    cmdlosoto = losoto + ' ' + parmdb + ' '
    noise = findamplitudenoise(parmdb)
    smooth = False
    if noise >= 0.1:
      cmdlosoto += create_losoto_mediumsmoothparset(ms, '9') 
      smooth = True      
    if noise < 0.1 and noise >= 0.08:
      cmdlosoto += create_losoto_mediumsmoothparset(ms, '7') 
      smooth = True    
    if noise < 0.08 and noise >= 0.07:
      cmdlosoto += create_losoto_mediumsmoothparset(ms, '5')
      smooth = True
    if noise < 0.07 and noise >= 0.04:
      cmdlosoto += create_losoto_mediumsmoothparset(ms, '3')
      smooth = True
    print(cmdlosoto)
    if smooth:
       os.system(cmdlosoto)
    return


def applycal(ms, parmdb, soltype, preapplyphase, TEC=False, weight_spectrum='WEIGHT_SPECTRUM_SOLVE',ncpu=None):

    # APPLYCAL CASE I (rare)
    if (soltype == 'complexgain') and (preapplyphase == False):
      cmd = 'DP3 '
      if ncpu is not None:
         cmd += 'numthreads=%i ' % ncpu
      cmd += 'msin=' + ms + ' msin.datacolumn=DATA msout=. '
      cmd += 'msin.weightcolumn='+weight_spectrum + ' '
      cmd += 'msout.datacolumn=CORRECTED_DATA steps=[ac1,ac2] msout.storagemanager=dysco '
      cmd += 'ac1.parmdb='+parmdb + ' ac2.parmdb='+parmdb + ' '
      cmd += 'ac1.type=applycal ac2.type=applycal '
      cmd += 'ac1.correction=phase000 ac2.correction=amplitude000 ' 
      print('DP3 applycal:', cmd)
      os.system(cmd)

    # APPLYCAL CASE II
    if soltype == 'scalarphase' and TEC == False:
      cmd = 'DP3 '
      if ncpu is not None:
         cmd += 'numthreads=%i ' % ncpu

      cmd += 'msin=' + ms + ' msin.datacolumn=DATA msout=. '
      cmd += 'msin.weightcolumn='+weight_spectrum + ' '
      cmd += 'msout.datacolumn=CORRECTED_DATA steps=[ac1] msout.storagemanager=dysco '
      cmd += 'ac1.parmdb=phaseonly'+parmdb + ' ac1.type=applycal '
      cmd += 'ac1.correction=phase000 '
      print('DP3 applycal:', cmd)
      os.system(cmd)
    
    # APPLYCAL CASE III  
    if soltype == 'scalarphase' and TEC == True:
      cmd = 'DP3 '
      if ncpu is not None:
         cmd += 'numthreads=%i ' % ncpu
      
      cmd += 'msin=' + ms + ' msin.datacolumn=DATA msout=. '
      cmd += 'msin.weightcolumn='+weight_spectrum + ' '
      cmd += 'msout.datacolumn=CORRECTED_DATA steps=[ac1,ac2] msout.storagemanager=dysco '
      cmd += 'ac1.parmdb=phaseonly'+parmdb + ' ac1.type=applycal '
      cmd += 'ac1.correction=phase000 '
      cmd += 'ac2.parmdb=phaseonly'+parmdb + ' ac2.type=applycal '
      cmd += 'ac2.correction=tec000 '
      print('DP3 applycal:', cmd)
      os.system(cmd)

    # APPLYCAL CASE IV      
    if (soltype == 'complexgain') and (preapplyphase == True):
       
      cmd = 'DP3 '
      if ncpu is not None:
         cmd += 'numthreads=%i ' % ncpu
     
      cmd += 'msin=' + ms + ' msin.datacolumn=DATA msout=. '
      cmd += 'msin.weightcolumn='+weight_spectrum + ' msout.storagemanager=dysco '
      if TEC == False:
        cmd += 'msout.datacolumn=CORRECTED_DATA steps=[ac0,ac1,ac2] '
        cmd += 'ac0.parmdb=phaseonly'+parmdb + ' ac1.parmdb='+parmdb + ' ac2.parmdb='+parmdb + ' '
        cmd += 'ac0.type=applycal ac1.type=applycal ac2.type=applycal '
        cmd += 'ac0.correction=phase000 ac1.correction=phase000 ac2.correction=amplitude000 ' 
      if TEC == True:
        cmd += 'msout.datacolumn=CORRECTED_DATA steps=[ac0,ac1,ac2,ac3] '
        cmd += 'ac0.parmdb=phaseonly'+parmdb + ' ac0.type=applycal ac0.correction=phase000 '
        cmd += 'ac1.parmdb=phaseonly'+parmdb + ' ac1.type=applycal ac1.correction=tec000 '
        cmd += 'ac2.parmdb='+parmdb + ' ac2.type=applycal ac2.correction=phase000 '
        cmd += 'ac3.parmdb='+parmdb + ' ac3.type=applycal ac3.correction=amplitude000 '
                    
      print('DP3 applycal:', cmd)
      os.system(cmd) 

    return


def change_refant(parmdb, soltab):
    '''
    Changes the reference antenna, if needed, for phase
    '''
    H5     = h5parm.h5parm(parmdb, readonly=False) 
    phases = H5.getSolset('sol000').getSoltab(soltab).getValues()[0]
    weights= H5.getSolset('sol000').getSoltab(soltab).getValues(weight=True)[0]
    
    #print 'SHAPE', np.shape(weights)#, np.size(weights[:,:,0,:,:])

    
    antennas = list(H5.getSolset('sol000').getSoltab(soltab).getValues()[1]['ant'])
    #print antennas
    
    idx0    = np.where((weights[:,:,0,:,:] == 0.0))[0]
    idxnan  = np.where((~np.isfinite(phases[:,:,0,:,:])))[0]
    
    #print idx0

    refant = ' '    
    if ((old_div(np.float(len(idx0)),np.float(np.size(weights[:,:,0,:,:])))) > 0.5) or ((old_div(np.float(len(idxnan)),np.float(np.size(weights[:,:,0,:,:])))) > 0.5):
      logging.info('Trying to changing reference anntena')
    

      for antennaid,antenna in enumerate(antennas[1::]):
            print(antenna)
            idx0    = np.where((weights[:,:,antennaid+1,:,:] == 0.0))[0]

            idxnan  = np.where((~np.isfinite(phases[:,:,antennaid+1,:,:])))[0]
            
            print(idx0, idxnan, ((old_div(np.float(len(idx0)),np.float(np.size(weights[:,:,antennaid+1,:,:]))))))
            
            if  ((old_div(np.float(len(idx0)),np.float(np.size(weights[:,:,antennaid+1,:,:])))) < 0.5) and ((old_div(np.float(len(idxnan)),np.float(np.size(weights[:,:,antennaid+1,:,:])))) < 0.5):
              logging.info('Found new reference anntena,' + str(antenna))
              refant = antenna
              break
    
    
    if refant != ' ':
        for antennaid,antenna in enumerate(antennas):
            phases[:,:,antennaid,:,:] = phases[:,:,antennaid,:,:] - phases[:,:,antennas.index(refant),:,:]
            
        H5.getSolset('sol000').getSoltab(soltab).setValues(phases)     

    H5.close()
    return


def calculate_solintnchan(compactflux):
    
    if compactflux >= 3.5:
        nchan = 5.
        solint_phase = 1.
        
    if compactflux <= 3.5:
        nchan = 5.
        solint_phase = 1.
  
    if compactflux <= 1.0:
        nchan= 10.
        solint_phase = 2
 
    if compactflux <= 0.75:
        nchan= 15.
        solint_phase = 3.

 
    #solint_ap = 100. / np.sqrt(compactflux)
    solint_ap = 120. /(compactflux**(1./3.)) # do third power-scaling
    #print solint_ap
    if solint_ap < 60.:
        solint_ap = 60.  # shortest solint_ap allowed
    if solint_ap > 180.:
        solint_ap = 180.  # longest solint_ap allowed
 
    if compactflux <= 0.4:
        nchan= 15.
        solint_ap = 180.
 
    return np.int(nchan), np.int(solint_phase), np.int(solint_ap)




def determine_compactsource_flux(fitsimage):
    
    hdul = fits.open(fitsimage)
    bmaj = hdul[0].header['BMAJ']
    bmin = hdul[0].header['BMIN']
    avgbeam = 3600.*0.5*(bmaj + bmin)
    pixsize = 3600.*(hdul[0].header['CDELT2'])
    rmsbox1 = np.int(7.*avgbeam/pixsize)
    rmsbox2 = np.int((rmsbox1/10.) + 1.)
    
    img = bdsf.process_image(fitsimage,mean_map='zero', rms_map=True, rms_box = (rmsbox1,rmsbox2))
    hdul.close()

    
    return img.total_flux_gaus
#print determine_compactsource_flux('imtry1_0-MFS-image.fits')
#sys.exit()


def getdeclinationms(ms):
    '''
    return approximate declination of pointing center of the ms
    input: a ms
    output: declination in degrees
    '''
    t = pt.table(ms +'/FIELD', readonly=True)
    direction = np.squeeze ( t.getcol('PHASE_DIR') )
    t.close()
    return 360.*direction[1]/(2.*np.pi)

#print getdeclinationms('1E216.dysco.sub.shift.avg.weights.set0.ms')
#sys.exit()

def declination_sensivity_factor(declination):
    '''
    compute sensitivy factor lofar data, reduced by delclination, eq. from G. Heald.
    input declination is units of degrees
    '''
    factor = 1./(np.cos(2.*np.pi*(declination - 52.9)/360.)**2)

    return factor

#print declination_sensivity_factor(-3.7)
#sys.exit()

def flaglowamps(parmdb, lowampval=0.1):
    '''
    flag bad amplitudes in H5 parmdb, those with values < lowampval
    '''
    H5 = h5parm.h5parm(parmdb, readonly=False) 
    amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
    idx = np.where(amps < lowampval)
    #print idx
    #print amps[idx]
    weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
    
    #print 'W', np.sum(weights[idx])
    weights[idx] = 0.0
    amps[idx] = 1.0
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(weights,weight=True)
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(amps)

    #also put phases weights and phases to zero
    phases =H5.getSolset('sol000').getSoltab('phase000').getValues()[0]
    weights_p = H5.getSolset('sol000').getSoltab('phase000').getValues(weight=True)[0]
    weights_p[idx] = 0.0
    phases[idx] = 0.0
    #print idx
    H5.getSolset('sol000').getSoltab('phase000').setValues(weights_p,weight=True)
    H5.getSolset('sol000').getSoltab('phase000').setValues(phases)
    
    H5.close()
    return


def flaghighgamps(parmdb, highampval=10.):
    '''
    flag bad amplitudes in H5 parmdb, those with values > highampval
    '''
    H5 = h5parm.h5parm(parmdb, readonly=False) 
    amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
    idx = np.where(amps > highampval)
    #print idx
    #print amps[idx]
    weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
    
    #print 'W', np.sum(weights[idx])
    weights[idx] = 0.0
    amps[idx] = 1.0
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(weights,weight=True)
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(amps)

    #also put phases weights and phases to zero
    phases =H5.getSolset('sol000').getSoltab('phase000').getValues()[0]
    weights_p = H5.getSolset('sol000').getSoltab('phase000').getValues(weight=True)[0]
    weights_p[idx] = 0.0
    phases[idx] = 0.0
    #print idx
    H5.getSolset('sol000').getSoltab('phase000').setValues(weights_p,weight=True)
    H5.getSolset('sol000').getSoltab('phase000').setValues(phases)
    
    H5.close()
    return

def flagbadamps(parmdb):
    '''
    flag bad amplitudes in H5 parmdb, those with amplitude==1.0
    '''
    H5 = h5parm.h5parm(parmdb, readonly=False) 
    amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
    idx = np.where(amps == 1.0)
    #print idx
    weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
    
    weights[idx] = 0.0
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(weights,weight=True)

    #also put phases weights and phases to zero
    phases =H5.getSolset('sol000').getSoltab('phase000').getValues()[0]
    weights_p = H5.getSolset('sol000').getSoltab('phase000').getValues(weight=True)[0]
    weights_p[idx] = 0.0
    phases[idx] = 0.0

    H5.getSolset('sol000').getSoltab('phase000').setValues(weights_p,weight=True)
    H5.getSolset('sol000').getSoltab('phase000').setValues(phases)
    
    H5.close()
    return



def normamps(parmdb):
    '''
    normalize amplitude solutions to one
    '''
    
    if len(parmdb) == 1:
      H5 = h5parm.h5parm(parmdb[0], readonly=False) 
      amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
      weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
      idx = np.where(weights != 0.0)
    
      amps = np.log10(amps)
      logging.info('Mean amplitudes before normalization: ' + str(10**(np.nanmean(amps[idx]))))
      amps = amps - (np.nanmean(amps[idx]))
      logging.info('Mean amplitudes after normalization: ' + str(10**(np.nanmean(amps[idx]))))
      amps = 10**(amps)
      #print np.nanmean(amps[idx])

      H5.getSolset('sol000').getSoltab('amplitude000').setValues(amps) 
      H5.close()

    else:
      #amps = []  
      for i, parmdbi in enumerate(parmdb):
          H5 = h5parm.h5parm(parmdbi, readonly=True) 
          ampsi = np.copy(H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0])
          weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
          idx = np.where(weights != 0.0)
          logging.info(parmdbi + '  Normfactor: '+ str(10**(np.nanmean(np.log10(ampsi[idx])))))
          if i == 0:
            amps = np.ndarray.flatten(ampsi[idx])
          else:
            amps = np.concatenate((amps, np.ndarray.flatten(ampsi[idx])),axis=0)

          #print np.shape(amps), parmdbi
          H5.close()
      normmin = (np.nanmean(np.log10(amps))) 
      logging.info('Global normfactor: ' + str(10**normmin))
      # now write the new H5 files
      for parmdbi in parmdb:  
         H5   = h5parm.h5parm(parmdbi, readonly=False) 
         ampsi = np.copy(H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0])
         ampsi = (np.log10(ampsi)) - normmin
         ampsi = 10**ampsi
         H5.getSolset('sol000').getSoltab('amplitude000').setValues(ampsi) 
         H5.close()
    return


#flagbadamps('test.h5')  
#flaglowamps('test.h5')
#flaghighgamps('test.h5')
#normamps(['P10.h5', 'P12.h5','P14.h5','P15.h5'])
#sys.exit()

def flaglowamps(parmdb, lowampval=0.1):
    '''
    flag bad amplitudes in H5 parmdb, those with values < lowampval
    '''
    H5 = h5parm.h5parm(parmdb, readonly=False) 
    amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
    idx = np.where(amps < lowampval)
    #print idx
    #print amps[idx]
    weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
    
    #print 'W', np.sum(weights[idx])
    weights[idx] = 0.0
    amps[idx] = 1.0
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(weights,weight=True)
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(amps)

    #also put phases weights and phases to zero
    phases =H5.getSolset('sol000').getSoltab('phase000').getValues()[0]
    weights_p = H5.getSolset('sol000').getSoltab('phase000').getValues(weight=True)[0]
    weights_p[idx] = 0.0
    phases[idx] = 0.0
    #print idx
    H5.getSolset('sol000').getSoltab('phase000').setValues(weights_p,weight=True)
    H5.getSolset('sol000').getSoltab('phase000').setValues(phases)
    
    #H5.getSolset('sol000').getSoltab('phase000').flush()
    #H5.getSolset('sol000').getSoltab('amplitude000').flush()
    H5.close()
    return


def flaghighgamps(parmdb, highampval=10.):
    '''
    flag bad amplitudes in H5 parmdb, those with values > highampval
    '''
    H5 = h5parm.h5parm(parmdb, readonly=False) 
    amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
    idx = np.where(amps > highampval)
    #print idx
    #print amps[idx]
    weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
    
    #print 'W', np.sum(weights[idx])
    weights[idx] = 0.0
    amps[idx] = 1.0
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(weights,weight=True)
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(amps)

    #also put phases weights and phases to zero
    phases =H5.getSolset('sol000').getSoltab('phase000').getValues()[0]
    weights_p = H5.getSolset('sol000').getSoltab('phase000').getValues(weight=True)[0]
    weights_p[idx] = 0.0
    phases[idx] = 0.0
    print(idx)
    H5.getSolset('sol000').getSoltab('phase000').setValues(weights_p,weight=True)
    H5.getSolset('sol000').getSoltab('phase000').setValues(phases)
    
    #H5.getSolset('sol000').getSoltab('phase000').flush()
    #H5.getSolset('sol000').getSoltab('amplitude000').flush()
    H5.close()
    return


def removenegativefrommodel(imagenames):
    '''
    replace negative pixel values in WSCLEAN model images with zeros
    '''
    for image in imagenames:
        print('remove negatives from model: ', image)
        hdul = fits.open(image)
        data = hdul[0].data
        
        data[np.where(data < 0.0)] = 0.0
        hdul[0].data = data
        hdul.writeto(image, overwrite=True)
        hdul.close()
    return

def makeimage(mslist, imageout, pixsize, imsize, channelsout, niter, robust, uvtaper=False, multiscale=True, predict=True, uvmin=' ', fitsmask=None, idg=False, deepmultiscale=False, ncpu=None):

    msliststring = ' '.join(map(str, mslist))
    os.system('rm -f ' + imageout + '-*.fits')
    imcol = 'CORRECTED_DATA'
    t = pt.table(mslist[0],readonly=True) # just test for first ms in mslist
    colnames =t.colnames()
    if 'CORRECTED_DATA' not in colnames: # for first imaging run
      imcol = 'DATA' 
    t.close()
    
    #baselineav = str (old_div(2.5e3*60000.*2.*np.pi *np.float(pixsize),(24.*60.*60*np.float(imsize))) )
    baselineav = str (old_div(1.87e3*60000.*2.*np.pi *np.float(1.5),(24.*60.*60*np.float(imsize))) )
    #baselineav = str (5e4*60000.*2.*np.pi *np.float(pixsize)/(24.*60.*60*np.float(imsize)) ) # this would mean 4sec for 20000 px image
    #imcol = 'DATA' 
    
    
    wsclean = 'wsclean'
       
       
    cmd = wsclean + ' '
    if ncpu is not None:
       cmd += '-j %i ' % ncpu
    #if not deepmultiscale:
    cmd += '-no-update-model-required -minuv-l 80 '
    cmd += '-size ' + imsize + ' ' + imsize + ' -reorder '
    cmd += '-weight briggs ' + robust + ' -weighting-rank-filter 3 -clean-border 1 '
    cmd += '-mgain 0.8 -fit-beam -data-column ' + imcol +' -join-channels -channels-out '
    cmd += channelsout + ' -padding 1.4 -auto-mask 2.5 -auto-threshold 1.0 '
    #cmd += '-parallel-deconvolution ' + str(np.int(imsize)/2) + ' ' 
    if multiscale:
       cmd += '-multiscale '+' -multiscale-scales 0,4,8,16,32,64 '
    if fitsmask != None:
      if os.path.isfile(fitsmask): 
        cmd += '-fits-mask '+ fitsmask + ' '
        #cmd += '-beam-shape 6arcsec 6arcsec 0deg ' 
      else:
        print('fitsmask: ', fitsmask, 'does not exist')
        sys.exit()
    if uvtaper:
       cmd += '-taper-gaussian 15arcsec '
    
    
    if idg:
      cmd += '-use-idg -grid-with-beam -use-differential-lofar-beam -idg-mode cpu '
      cmd += '-beam-aterm-update 800 '
      cmd += '-pol iquv -link-polarizations i '
    else:
      cmd += '-fit-spectral-pol 3 '# -beam-shape 6arcsec 6arcsec 0deg '   
      cmd += '-pol i '
      cmd += '-baseline-averaging ' + baselineav + ' '
      
    
    cmd += '-name ' + imageout + ' -scale ' + pixsize + 'arcsec ' 

    print('WSCLEAN: ', cmd + '-niter ' + str(niter) + ' ' + msliststring)
    logging.info(cmd + '-niter ' + str(niter) + ' ' + msliststring)
    os.system(cmd + '-niter ' + str(niter) + ' ' + msliststring)

    if deepmultiscale:
        
      # predict first to fill MODEL_DATA so we can continue with clean
      cmdp = wsclean + ' -size ' 
      cmdp += imsize + ' ' + imsize + ' -channels-out ' + channelsout + ' -padding 1.4 -predict ' 
      if idg:
        cmdp += '-use-idg -grid-with-beam -use-differential-lofar-beam -idg-mode cpu '
        cmdp += '-beam-aterm-update 800 '
        cmdp += '-pol iquv '
      
      cmdp += '-name ' + imageout + ' -scale ' + pixsize + 'arcsec ' + msliststring
      print('PREDICT STEP for continue: ', cmdp)
      os.system(cmdp)
       
      # NOW continue cleaning  
      cmd += '-niter ' + str(old_div(niter,5)) + ' -multiscale -continue ' + msliststring
      print('WSCLEAN continue: ', cmd)
      os.system(cmd)

    # REMOVE nagetive model components, these are artifacts (only for Stokes I)
    if idg:
      removenegativefrommodel(sorted(glob.glob(imageout +'-????-I-model*.fits')))  # only Stokes I
    else:    
      removenegativefrommodel(sorted(glob.glob(imageout + '-????-model.fits')))

    if predict:
      cmd = wsclean + ' -size ' 
      cmd += imsize + ' ' + imsize + ' -channels-out ' + channelsout + ' -padding 1.4 -predict ' 
      if idg:
        cmd += '-use-idg -grid-with-beam -use-differential-lofar-beam -idg-mode cpu '
        cmd += '-beam-aterm-update 800 '
        cmd += '-pol iquv '
      
      cmd += '-name ' + imageout + ' -scale ' + pixsize + 'arcsec ' + msliststring
      print('PREDICT STEP: ', cmd)
      os.system(cmd)

 
def runDP3(ms, solint_ap, solint_phaseonly, nchan_phase, nchan_ap, parmdb, soltype, \
             preapplyphase, weight_spectrum='WEIGHT_SPECTRUM_SOLVE',uvmin=0, TEC=False,\
             FOVedge=False, smoothcal=True, ncpu=None):
    
    losotoparset = create_losoto_flag_apgridparset(ms)
    losotoparset_phase = create_losoto_fastphaseparset(ms)
    losotoparset_tecandphase = create_losoto_tecandphaseparset(ms)


    losoto = 'losoto'
       

    
    if os.path.isfile(parmdb):
      print('H5 file exists  ', parmdb)
      #os.system('rm -f ' + parmdb)
    if os.path.isfile('phaseonly' + parmdb) and preapplyphase == True:
      print('H5 file exists  ', 'phaseonly' + parmdb)
      #os.system('rm -f ' + parmdb)
    
    if soltype == 'scalarphase' and preapplyphase == True:
        print('Not supported')
        sys.exit()
    
 
    cmd = 'DP3 '
    if ncpu is not None:
       cmd += 'numthreads=%i ' % ncpu

    cmd += 'msin=' + ms + ' msin.datacolumn=DATA msout=. msin.modelcolumn=MODEL_DATA '
    cmd += 'msin.weightcolumn='+weight_spectrum + ' '
    cmd += 'steps=[ddecal] ' + 'msout.storagemanager=dysco ddecal.type=ddecal '
    cmd += 'ddecal.maxiter=100 ddecal.propagatesolutions=True '
    cmd += 'ddecal.usemodelcolumn=True '
    if uvmin != 0:
        cmd += 'ddecal.uvlambdamin=' + str(uvmin) + ' '      
   
    # CASE I   
    if soltype == 'complexgain' and preapplyphase == True:
        if TEC == False:
           cmd += 'ddecal.mode=scalarphase ' # scalarphase, phaseonly, complexgain, tec, tecandphase
           cmd += 'ddecal.tolerance=1.e-5 '
        if TEC == True:
           cmd += 'ddecal.mode=tecandphase ' # scalarphase, phaseonly, complexgain, tec, tecandphase
           cmd += 'ddecal.approximatetec=True '
           cmd += 'ddecal.stepsize=0.2 '
           cmd += 'ddecal.maxapproxiter=45 '
           cmd += 'ddecal.tolerance=1e-4 '
           cmd += 'ddecal.approxtolerance=6e-3 '
        
        cmd += 'ddecal.solint=' + str(solint_phaseonly) + ' '
        cmd += 'ddecal.nchan=' + str(nchan_phase) + ' '
        cmd += 'ddecal.h5parm=phaseonly' + parmdb + ' '

    # CASE II (not really that useful)  
    if soltype == 'complexgain' and preapplyphase == False:
        cmd += 'ddecal.mode=complexgain ' # scalarphase, phaseonly
        cmd += 'ddecal.solint=' + str(solint_ap) + ' '
        cmd += 'ddecal.nchan=' + str(nchan_ap) + ' '
        cmd += 'ddecal.h5parm=' + parmdb + ' '
        cmd += 'ddecal.tolerance=1.e-5 '
        if FOVedge:
          cmd += 'ddecal.smoothnessconstraint=8e6 '
 
    # CASE III      
    if soltype == 'scalarphase' and preapplyphase == False:
        if TEC == False:
           cmd += 'ddecal.mode=scalarphase ' # scalarphase, phaseonly
           cmd += 'ddecal.tolerance=1.e-5 '
        if TEC == True:
           cmd += 'ddecal.mode=tecandphase ' #  tec, tecandphase
           cmd += 'ddecal.approximatetec=True '
           cmd += 'ddecal.stepsize=0.2 '
           cmd += 'ddecal.maxapproxiter=45 '
           cmd += 'ddecal.tolerance=1.e-4 '
           cmd += 'ddecal.approxtolerance=6e-3 '
        
        cmd += 'ddecal.solint=' + str(solint_phaseonly) + ' '
        cmd += 'ddecal.nchan=' + str(nchan_phase) + ' '
        cmd += 'ddecal.h5parm=phaseonly' + parmdb + ' ' 
        
 
    print('DP3 solve:', cmd)
    os.system(cmd)
    #sys.exit()  
    if preapplyphase: # APPLY FIRST 
        cmd = 'DP3 '
        if ncpu is not None:
           cmd += 'numthreads=%i ' % ncpu
        cmd += 'msin=' + ms + ' msin.datacolumn=DATA msout=. '
        cmd += 'msin.weightcolumn='+weight_spectrum + ' msout.storagemanager=dysco '
        if TEC == False:
          cmd += 'msout.datacolumn=CORRECTED_DATA_PHASE steps=[ac1] '
          cmd += 'ac1.parmdb=phaseonly'+parmdb + ' ac1.type=applycal '
          cmd += 'ac1.correction=phase000 '
          print('DP3 PRE-APPLY PHASE-ONLY:', cmd)
          os.system(cmd)
          cmdlosotophase = losoto + ' '
          cmdlosotophase += 'phaseonly' + parmdb + ' ' + losotoparset_phase
          os.system(cmdlosotophase)
        if TEC == True:
          cmd += 'msout.datacolumn=CORRECTED_DATA_PHASE steps=[ac1,ac2] '
          cmd += 'ac1.parmdb=phaseonly'+parmdb + ' ac1.type=applycal '
          cmd += 'ac1.correction=phase000 '
          cmd += 'ac2.parmdb=phaseonly'+parmdb + ' ac2.type=applycal '
          cmd += 'ac2.correction=tec000 '
          print('DP3 PRE-APPLY TECANDPHASE:', cmd)
          os.system(cmd)
          cmdlosotophase = losoto + ' '
          cmdlosotophase += 'phaseonly' + parmdb + ' ' + losotoparset_tecandphase
          os.system(cmdlosotophase)

        # RUN DP3 again
        cmd = 'DP3 '
        if ncpu is not None:
           cmd += 'numthreads=%i ' % ncpu

        cmd += 'msin=' + ms + ' msin.datacolumn=CORRECTED_DATA_PHASE msout=. '
        cmd += 'msin.weightcolumn='+weight_spectrum + ' '
        cmd += 'msin.modelcolumn=MODEL_DATA '
        cmd += 'steps=[ddecal] ' + 'msout.storagemanager=dysco ddecal.type=ddecal '
        cmd += 'ddecal.maxiter=100 ddecal.tolerance=1.e-5 ddecal.propagatesolutions=True '
        cmd += 'ddecal.usemodelcolumn=True '
        cmd += 'ddecal.nchan=' + str(nchan_ap) + ' '
        cmd += 'ddecal.mode=' + soltype + ' ' # scalarphase, phaseonly, complexgain, tec, tecandphase
        cmd += 'ddecal.h5parm=' + parmdb + ' ' 
        cmd += 'ddecal.solint=' + str(solint_ap) + ' ' 
        if FOVedge:
          cmd += 'ddecal.smoothnessconstraint=8e6 '

        if uvmin != 0:
           cmd += 'ddecal.uvlambdamin=' + str(uvmin) + ' '
        print('DP3 SLOW GAIN solve:', cmd)
        os.system(cmd)
        os.system('cp ' + parmdb + ' ' + parmdb + '.backup')

    #  ---- PROCESS SOLUTIONS, FILTERING AND PLOTTING ----
    # MAKE losoto command
    cmdlosoto = losoto + ' '
    
    cmdlosoto += parmdb + ' ' + losotoparset


    
    #  CASE I (rare)
    if (soltype == 'complexgain') and (preapplyphase == False):
      flagbadamps(parmdb)  
      flaglowamps(parmdb)
      flaghighgamps(parmdb)
      change_refant(parmdb,'phase000')
      removenans(parmdb, 'amplitude000')
      removenans(parmdb, 'phase000')
      # FLAG/SMOOTH solutions
      os.system(cmdlosoto)
      
      if smoothcal:
        smoothsols(parmdb, ms)
      

    #  CASE II
    if soltype == 'scalarphase' and TEC == False:
      cmdlosotophase = losoto + ' '
      cmdlosotophase += 'phaseonly' + parmdb + ' ' + losotoparset_phase
      os.system(cmdlosotophase)
    
    # CASE III  
    if soltype == 'scalarphase' and TEC == True:
      cmdlosotophase = losoto + ' '
      cmdlosotophase += 'phaseonly' + parmdb + ' ' + losotoparset_tecandphase
      os.system(cmdlosotophase)

    #  CASE IV      
    if (soltype == 'complexgain') and (preapplyphase == True):
      flagbadamps(parmdb)
      if FOVedge == True:
        flaglowamps(parmdb,lowampval=0.1)
        flaghighgamps(parmdb, highampval=10.)  
      else:
        flaglowamps(parmdb,lowampval=0.1)
        flaghighgamps(parmdb, highampval=10.)
      
      # FLAG/SMOOTH solutions
      removenans(parmdb, 'amplitude000')
      removenans(parmdb, 'phase000')      
      change_refant(parmdb,'phase000')
      os.system(cmdlosoto)
      
      if smoothcal:
         smoothsols(parmdb, ms)
      
      #if os.path.isdir('plotlosoto' + ms):
      #  os.system('rm -rf plotlosoto' + ms)
      #os.system('mv plotlosoto plotlosoto' + ms) 
       
#plotimage('imselfcal_0-MFS-image.fits', 'test.png')
#sys.exit()

# ------------- MAIN LOOP ------------------

if __name__=='__main__':

   switchtogaincycle = 3 # number of phase-only selfcal cycles


   # ---- INPUT -----#

   parser = argparse.ArgumentParser(description='Calibrate and image DR2 cutout')
   parser.add_argument('-b','--boxfile', help='boxfile, required argument', type=str)
   parser.add_argument('--fitsmask', help='fitsmask for deconvolution, if not provided use automasking', type=str)
   parser.add_argument('--H5sols', help='prefix name for H5 solution file, default=solsgrid', default='solsgrid', type=str)
   parser.add_argument('--imsize', help='image size, required if boxfile is not used', type=int)
   parser.add_argument('-n', '--niter', help='niter, default=15000', default=15000, type=int)
   parser.add_argument('--robust', help='Briggs robust paramter, default=-0.5', default=-0.5, type=float)
   parser.add_argument('--channelsout', help='channelsout, default=6', default=6, type=int)
   parser.add_argument('-u', '--uvmin', help='inner uv-cut for calibration in lambda, default=350', default=350., type=float)
   parser.add_argument('--no-tec', help='do not use TEC fitting', action='store_false')
   parser.add_argument('--multiscale', help='use multiscale deconvolution, not recommended', action='store_true')
   parser.add_argument('--pixelscale', help='pixels size in arcsec, deafult=1.5', default=1.5, type=float)
   parser.add_argument('--idg', help='use the Image Domain gridder', action='store_true')
   parser.add_argument('--no-beamcor', help='do not correct the visilbities for the array factor', action='store_false')
   parser.add_argument('-i','--imagename', help='imagename, default=image', required=True, type=str)
   parser.add_argument('--start', help='start selfcal cycle at this iteration, default=0', default=0, type=int)
   parser.add_argument('--stop', help='stop selfcal cycle at this iteration, default=10', default=10, type=int)
   parser.add_argument('--no-smoothcal', help='median smooth amplitudes', action='store_false')
   parser.add_argument('--maskthreshold', help='threshold for MakeMask.py, default=5', default=5, type=int)
   if getcpuworks:
      parser.add_argument('--ncpu', help='number of cpu to use, default=%i' % getcpus(), default=getcpus(), type=int)
   else:
      parser.add_argument('--ncpu', help='number of cpu to use, default=%i' % os.cpu_count(), default= os.cpu_count(), type=int)

   parser.add_argument('ms', nargs='*', help='msfile(s)')  

   args = vars(parser.parse_args())
   #print args

   if which('DP3') == None:
     print('Cannot find DP3, forgot to source lofarinit.[c]sh?')
     sys.exit()



   if args['boxfile'] == None and args['imsize'] == None:
     print('Incomplete input detected, either boxfile or imsize is required')
     sys.exit()
   if args['boxfile'] != None and args['imsize'] != None:
     print('Wrong input detected, both boxfile and imsize are set')
     sys.exit()


   mslist = sorted(args['ms'])
   if args['boxfile'] != None:
     imsize   = str(getimsize(args['boxfile'], args['pixelscale']))
   if args['imsize'] != None:
     imsize = str(args['imsize']) 
   TEC = args['no_tec']
   idg = args['idg']
   ncpu = args['ncpu']
   multiscale = args['multiscale']
   imageout  = args['imagename'] + '_'
   dobeamcor = args['no_beamcor']
   if args['fitsmask'] != None:
     fitsmask = args['fitsmask']
   else:
     fitsmask = None
   uvmin = args['uvmin']  
   robust = str(args['robust'])
   channelsout = str(args['channelsout'])  
   niter = args['niter']
   parmdb = args['H5sols']  + '_'
   pixsize = str(args['pixelscale'])  

   #print args['no_smoothcal']
   if args['boxfile'] != None:
     outtarname = (args['boxfile'].split('/')[-1]).split('.reg')[0] + '.tar.gz'
   else:
     outtarname = 'calibrateddata' + '.tar.gz' 

   logging.info('Imsize:                    ' + str(imsize))
   logging.info('Pixelscale:                ' + str(pixsize))
   logging.info('Niter:                     ' + str(niter))
   logging.info('Uvmin:                     ' + str(uvmin))
   logging.info('Multiscale:                ' + str(multiscale))
   logging.info('Beam correction:           ' + str(dobeamcor))
   logging.info('IDG:                       ' + str(idg))
   logging.info('TEC:                       ' + str(TEC))
   if args['boxfile'] != None:
     logging.info('Bobxfile:                  ' + args['boxfile'])
   logging.info('Mslist:                    ' + ' '.join(map(str,mslist)))
   logging.info('User specified clean mask: ' + str(fitsmask))
   logging.info('Threshold for MakeMask:    ' + str(args['maskthreshold']))
   logging.info('Briggs robust:             ' + str(robust))
   logging.info('Imagename prefix:          ' + imageout)
   logging.info('Solution file prefix:      ' + parmdb)
   logging.info('Output file will be:       ' + outtarname)


   makemask = 'MakeMask.py'

   for ms in mslist:
     if not os.path.isdir(ms):
       print(ms, ' does not exist')
       sys.exit()

   if beamcor and idg:
     print('beamcor=True and IDG=True is not possible')
     sys.exit()


   deepmultiscale = False

   # GET SOLUTION TIMSCALES
   nchan_phase,solint_phase,solint_ap,nchan_ap = determinesolints(mslist, \
                                                 pixsize, imsize, channelsout, \
                                                 np.int(niter), robust, TEC)

   # ----- START SELFCAL LOOP -----
   for i in range(args['start'],args['stop']):

     # AUTOMATICALLY PICKUP PREVIOUS MASK (in case of a restart)
     if (i > 0) and (args['fitsmask'] == None):
       if idg:  
         if os.path.isfile(imageout + str(i-1) + '-MFS-I-image.fits.mask.fits'):
             fitsmask = imageout + str(i-1) + '-MFS-I-image.fits.mask.fits'
       else:
         if os.path.isfile(imageout + str(i-1) + '-MFS-image.fits.mask.fits'):
             fitsmask = imageout + str(i-1) + '-MFS-image.fits.mask.fits'


     # BEAM CORRECTION
     if dobeamcor and i == 0:
         for ms in mslist:
           beamcor(ms,ncpu=ncpu)

     # do an additional clean run with "-continue" using multiscale
     #if i>= 6 and not multiscale:
     #    deepmultiscale=True


     # IMAGE
     makeimage(mslist, imageout + str(i), pixsize, imsize, channelsout, np.int(niter), robust, uvtaper=False, multiscale=multiscale, idg=idg, fitsmask=fitsmask, deepmultiscale=deepmultiscale, ncpu=ncpu)

     # MAKE FIGURE WITH APLPY
     if idg:
       plotimage(imageout + str(i) +'-MFS-I-image.fits',imageout + str(i) + '.png' , \
                 mask=fitsmask, rmsnoiseimage=imageout + str(0) +'-MFS-I-image.fits')
     else:
       plotimage(imageout + str(i) +'-MFS-image.fits',imageout + str(i) + '.png' , \
                 mask=fitsmask, rmsnoiseimage=imageout + str(0) +'-MFS-image.fits')


     # SOLVE
     for msnumber, ms in enumerate(mslist):
       if i < switchtogaincycle:
         runDP3(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
                  np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
                  ms + parmdb + str(i) + '.h5' ,'scalarphase', False, uvmin=uvmin, TEC=TEC, FOVedge=False, ncpu=ncpu)
       else:
         runDP3(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
                  np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
                  ms + parmdb + str(i) + '.h5'  ,'complexgain', True, uvmin=uvmin, TEC=TEC, FOVedge=False, smoothcal=args['no_smoothcal'], ncpu=ncpu)

     # NORMALIZE GLOBAL GAIN (done in log-space)
     if i >= switchtogaincycle:
        print('Doing global gain normalization')  
        parmdblist = []  
        for msnumber, ms in enumerate(mslist):
          parmdblist.append(ms + parmdb + str(i) + '.h5')  
        normamps(parmdblist)

     # APPLYCAL
     for msnumber, ms in enumerate(mslist):
       if i < switchtogaincycle:
         applycal(ms, ms + parmdb + str(i) +'.h5' ,'scalarphase', False, TEC=TEC, ncpu=ncpu)
       else:
         applycal(ms, ms + parmdb + str(i) +'.h5' ,'complexgain', True, TEC=TEC, ncpu=ncpu)   

     # MAKE MASK
     if args['fitsmask'] == None:
       if idg:  
         imagename  = imageout + str(i) + '-MFS-I-image.fits'
       else:
         imagename  = imageout + str(i) + '-MFS-image.fits'
       cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
       os.system(cmdm)
       fitsmask = imagename + '.mask.fits'

     # CUT FLAGGED DATA FROM MS AT START AND END
     if (i == 0) or (i == switchtogaincycle) or (i == switchtogaincycle + 1) or (i == switchtogaincycle + 2) \
       or (i == switchtogaincycle + 3) or (i == switchtogaincycle + 4):
        for msnumber, ms in enumerate(mslist):  
          flagms_startend(ms, 'phaseonly' + ms + parmdb + str(i) + '.h5', np.int(solint_phase[msnumber]))



   archive(mslist, outtarname, args['boxfile'], fitsmask, imagename)    
   cleanup(mslist)

