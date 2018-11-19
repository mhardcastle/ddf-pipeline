import pyrap.tables as pt
import os,sys
import numpy as np
import argparse
import pyregion
from astropy.io import fits
from astropy.wcs import WCS
from astropy.io import ascii
import glob

# copy SOLSDIR, big-mslist.txt, 'image_full_ampphase_di_m.NS.mask01.fits', 'image_full_ampphase_di_m.NS.DicoModel', image_dirin_SSD_m.npy.ClusterCat.npy


# NOTE, applybeam NDPPP step does not work on phase-shifted data, do not use it.


def fixsymlinks():
  # need to make it....
  os.path.islink('filecheck')
  return



def columnchecker(mslist, colname):
    
    for ms in mslist:
      t = pt.table(ms, ack=False) 
      if colname not in t.colnames():
          print colname, ' not present in ', ms
          sys.exit()
      t.close()

def filechecker():
  '''
  Check if files are present to avoid errors to avoid crashes
  '''
  if not os.path.isfile('image_full_ampphase_di_m.NS.DicoModel'):
    raise IOError('image_full_ampphase_di_m.NS.DicoModel does not exist')
  if not os.path.isfile('image_full_ampphase_di_m.NS.mask01.fits'):
    raise IOError('image_full_ampphase_di_m.NS.mask01.fits does not exist')   
  if not os.path.isfile('image_dirin_SSD_m.npy.ClusterCat.npy'):
    raise IOError('image_dirin_SSD_m.npy.ClusterCat.npy does not exist')   
  if not os.path.isdir('SOLSDIR'):
   raise IOError('SOLSDIR directory does not exist')
  return


def striparchivename():
  mslist = glob.glob('L*_SB*.ms.archive')
  for ms in mslist:
      outname = ms.rstrip('.archive')
      cmd = 'mv ' + ms + ' ' + outname
      print (cmd)
      os.system(cmd)

  return


def addextraweights(msfiles):
   '''
   Adds the column WEIGHT_SPECTRUM_FROM_IMAGING_WEIGHT from IMAGING_WEIGHT from DR2
   Input msfiles (list of ms)
   '''

   for ms in msfiles:
     ts  = pt.table(ms, readonly=False)
     colnames = ts.colnames()
     if 'WEIGHT_SPECTRUM_FROM_IMAGING_WEIGHT' not in colnames:
       desc = ts.getcoldesc('WEIGHT_SPECTRUM')
       desc['name']='WEIGHT_SPECTRUM_FROM_IMAGING_WEIGHT'
       ts.addcols(desc)
       ts.close() # to write results

     else:
         print 'WEIGHT_SPECTRUM_FROM_IMAGING_WEIGHT already exists'
         ts.close()
         
     ts  = pt.table(ms, readonly=False)    

     if 'IMAGING_WEIGHT' in colnames:
       iw = ts.getcol('IMAGING_WEIGHT')
       ws_tmp = ts.getcol('WEIGHT_SPECTRUM_FROM_IMAGING_WEIGHT')
       n, nfreq, npol = np.shape(ws_tmp)

       for i in range(npol):
         print 'Copying over correlation ', i, ms
         ws_tmp[:,:,i] = iw
         ts.putcol('WEIGHT_SPECTRUM_FROM_IMAGING_WEIGHT',ws_tmp)
     else:
         print 'IMAGING_WEIGHT column is not present in:', ms
     ts.close()

   return


def mask_region(infilename,ds9region,outfilename):

    hdu=fits.open(infilename)
    hduflat = flatten(hdu)
    map=hdu[0].data

    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)
    rmsval = np.mean(hdu[0].data[0][0][np.where(manualmask == True)])
    hdu[0].data[0][0][np.where(manualmask == True)] = 0.0
    hdu.writeto(outfilename,overwrite=True)

    return

def mask_except_region(infilename,ds9region,outfilename):

    hdu=fits.open(infilename)
    hduflat = flatten(hdu)
    map=hdu[0].data

    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)
    hdu[0].data[0][0][np.where(manualmask == False)] = 0.0
    hdu.writeto(outfilename,overwrite=True)

    return

def flatten(f):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        return fits.PrimaryHDU(header=f[0].header,data=f[0].data)

    w = WCS(f[0].header)
    wn=WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2
    copy=('EQUINOX','EPOCH','BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
    for k in copy:
        r=f[0].header.get(k)
        if r is not None:
            header[k]=r

    slice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            slice.append(np.s_[:],)
        else:
            slice.append(0)
        
    hdu = fits.PrimaryHDU(header=header,data=f[0].data[tuple(slice)])
    return hdu

def removecolumn(msfile,colname):
     t = pt.table(msfile,readonly=False)
     colnames =t.colnames()
     if colname in colnames:
        print 'Removing ',  colname, 'from ', msfile
        t.removecols(colname)
     t.close()
     return


def getregionboxcenter(regionfile):
    """
    Extract box center of a DS9 box region. 
    Input is regionfile Return NDPPP compatible string for phasecenter shifting
    """
    r = pyregion.open(regionfile)
    
    if len(r[:]) > 1:
      print 'Only one region can be specified, your file contains', len(r[:])
      sys.exit() 
    
    if r[0].name != 'box':
      print 'Only box region supported'
      sys.exit()
    
    ra  = r[0].coord_list[0]
    dec = r[0].coord_list[1]
    boxsizex = r[0].coord_list[2]
    boxsizey = r[0].coord_list[3]
    angle = r[0].coord_list[4]
    if boxsizex != boxsizey:
      print 'Only a sqaure box region supported, you have these sizes:', boxsizex, boxsizey
      sys.exit()
    if np.abs(angle) > 1:
      print 'Only nomrally oriented sqaure boxes are supported, you region is oriented under angle:', angle
      sys.exit()   
    
    regioncenter =  ('{:12.8f}'.format(ra) + 'deg,' + '{:12.8f}'.format(dec) + 'deg').replace(' ', '')
    return regioncenter

def mscolexist(ms, colname):
    """ Check if a colname exists in the measurement set ms, returns either True or False """
    t = pt.table(ms,readonly=True)
    colnames =t.colnames()
    if colname in colnames: # check if the column is in the list
      exist = True
    else:
      exist = False  
    t.close()
    return exist




# Use MakeMask to remove sources within the cluster from the mask.
fullmask    = 'image_full_ampphase_di_m.NS.mask01.fits'
indico      = 'image_full_ampphase_di_m.NS.DicoModel'
outdico     = 'image_full_ampphase_di_m_SUB.NS.DicoModel'

parser = argparse.ArgumentParser(description='Keep soures insize box region, subtract everything else and create new ms')
parser.add_argument('-b','--boxfile', help='boxfile, required argument', required=True, type=str)
parser.add_argument('-m','--mslist', help='DR2 mslist file, default=big-mslist.txt', default='big-mslist.txt', type=str)
parser.add_argument('-c','--column', help='Input column for the ms, default=DATA', default='DATA', type=str) #DATA_DI_CORRECTED
parser.add_argument('-f','--freqavg', help='channel averaging, default=4', default=4, type=int)
parser.add_argument('-t','--timeavg', help='timesample averaging, default=2', default=2, type=int)
parser.add_argument('-n','--ncpu', help='number of cpu to use, default=34', default=34, type=int)
parser.add_argument('-p','--prefixname', help='prefixname for output ms, default=object', default='object', type=str)
parser.add_argument('--nodysco', help='Do not dysco compress output', action='store_false')
parser.add_argument('--split', help='Do not concat but keep 10 SB blocks', action='store_true')
parser.add_argument('--aoflaggerbefore', help='Do an extra round of AOflagger on input data', action='store_true')
parser.add_argument('--aoflaggerafter', help='Do an extra round of AOflagger on averaged output data', action='store_true')
parser.add_argument('--maxamplitude', help='flag amplitudes above this number, default=1e6', default=1.e6, type=float)
#parser.add_argument('--takeoutbeam', help='Correct for the beam on the phase-shifted target data', action='store_true')

args = vars(parser.parse_args())

striparchivename()

if not os.path.isfile(args['mslist']):
    # try to make it
    os.system('/net/rijn/data2/rvweeren/LoTSS_ClusterCAL/make_mslists.py')
    if not os.path.isfile(args['mslist']):
      raise IOError('File', args['mslist'], 'does not exist and coult not be created')

filechecker()

solsfile = glob.glob('DDS3_full*smoothed.npz')
if len(solsfile) != 1:
     print 'Cannot find the correct solution file'
     sys.exit()
solsfile = str(solsfile[0])


boxfile     = args['boxfile']
ncpu        = args['ncpu']
timestepavg = args['timeavg']
freqstepavg = args['freqavg']
obsid       = args['prefixname']


dopredict   = True
dosubtract  = True
doconcat    = True
dokmscal     = True
doflagafter = args['aoflaggerafter']
amplmax     = args['maxamplitude']  # flag data with amplitues above this number
takeoutbeam = False # not supported by NDPPP #args['takeoutbeam']


aoflagger   = args['aoflaggerbefore']
dysco       = args['nodysco']
split       = args['split']  # ouput seperate ms for DDF pipeline

colname = 'DATA_SUB'

#print doflagafter, takeoutbeam, aoflagger, dysco, split


msfiles   = ascii.read(args['mslist'],data_start=0)
msfiles   = list(msfiles[:][msfiles.colnames[0]]) # convert to normal list of strings

t = pt.table(msfiles[0] + '/OBSERVATION')
fieldname = t.getcol('LOFAR_TARGET')['array'][0]
t.close()

msoutconcat   = fieldname + '_' + obsid + '.dysco.sub.shift.avg.weights.ms.archive'
phasecenter = '[' + getregionboxcenter(boxfile) + ']'
print phasecenter

outmask = 'cutoutmask.fits' #boxfile.split('.reg')[0] + 'mask' # just a name, can be anything


if os.path.isdir(msoutconcat):
  print 'MS exists:', msoutconcat
  sys.exit()
if os.path.isdir(msoutconcat+'.tmpweight'):
  print 'MS exists:', msoutconcat+'.tmpweight'
  sys.exit()
  


columnchecker(msfiles, args['column'])


if dopredict:
    
    # Apparently it can be dangerous to remove a column for tiled storagemanagers, comment out!
    #for ms in msfiles:
    #  removecolumn(ms, 'PREDICT_SUB')
    #  removecolumn(ms, 'DATA_SUB')
    
    os.system('rm -f ' + outdico)  # clean up
    os.system('rm -f ' + outmask)
    mask_region(fullmask,boxfile,outmask)
    
    os.system("MaskDicoModel.py --MaskName=%s --InDicoModel=%s --OutDicoModel=%s"%(outmask,indico,outdico))

    os.system("DDF.py --Output-Name=image_full_ampphase_di_m.NS_SUB --Data-MS=" + args['mslist'] + " --Deconv-PeakFactor 0.001000 --Data-ColName " + args['column'] + " --Parallel-NCPU="+str(ncpu) + " --Facets-CatNodes=image_dirin_SSD_m.npy.ClusterCat.npy --Beam-CenterNorm=1 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust -0.500000 --Image-NPix=20000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 1.500000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --Cache-Weight=reset --Output-Mode=Predict --Output-RestoringBeam 6.000000 --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=5.00 --Mask-External=" + outmask + " --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=" + solsfile + " --Predict-InitDicoModel=" + outdico + " --Selection-UVRangeKm=[0.100000,1000.000000] --GAClean-MinSizeInit=10 --Cache-Reset 1 --Beam-Smooth=1 --Predict-ColName='PREDICT_SUB'")


# Subtract the columns
if dosubtract:
    for ms in msfiles:
        t = pt.table(ms,readonly=False, ack=True)
        colnames =t.colnames()

        if ('PREDICT_SUB' in colnames) and (args['column'] in colnames):
            print 'Reading', 'PREDICT_SUB'
            f=t.getcol('PREDICT_SUB')
            print 'Reading', args['column']
            d=t.getcol(args['column'])


            if colname not in colnames:
               # Append new column containing all sources
               desc = t.getcoldesc(args['column'])
               desc['name']= colname
               t.addcols(desc)
            print 'Writing %s'%colname
            t.putcol(colname,d-f)
        else:
            print 'Warning, ', msfile, ' does not contain PREDICT_SUB and/or ' + args['column'] +', skipping.....'
        
        t.close()

    addextraweights(msfiles)

if dokmscal:
  outmask_target = 'inregionmask.fits'
  outdico_target =   'image_full_ampphase_di_m_TAR.NS.DicoModel'

  mask_except_region(fullmask,boxfile,outmask_target)
    
  os.system("MaskDicoModel.py --MaskName=%s --InDicoModel=%s --OutDicoModel=%s"%(outmask_target,indico,outdico_target))

  os.system("DDF.py --Output-Name=image_full_ampphase_di_m.NS_TAR --Data-MS=" + args['mslist'] + " --Deconv-PeakFactor 0.001000 --Data-ColName " + args['column'] + " --Parallel-NCPU="+str(ncpu) + " --Facets-CatNodes=image_dirin_SSD_m.npy.ClusterCat.npy --Beam-CenterNorm=1 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust -0.500000 --Image-NPix=20000 --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell 1.500000 --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --Cache-Weight=reset --Output-Mode=Predict --Output-RestoringBeam 6.000000 --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=5.00 --Mask-External=" + outmask + " --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=" + solsfile + " --Predict-InitDicoModel=" + outdico + " --Selection-UVRangeKm=[0.100000,1000.000000] --GAClean-MinSizeInit=10 --Cache-Reset 1 --Beam-Smooth=1 --Predict-ColName='PREDICT_TAR'")
  
  for ms in msfiles:



    os.system('kMS.py --MSName %s --SolverType KAFCA --PolMode Scalar --BaseImageName Predict_DDT --dt 0.536871 --NIterKF 6 --CovQ 0.100000 --LambdaKF=0.500000 --NCPU 32 --OutSolsName DIT --PowerSmooth=0.000000 --InCol DATA_SUB --Weighting Natural --UVMinMax=0.100000,1000.000000 --SolsDir=SOLSDIR --SolverType CohJones --PolMode Scalar --SkyModelCol PREDICT_TAR --OutCol DATA_SUB_CORRECTED --ApplyToDir 0 --dt 1.0 --NChanSols 1'%(ms))
  colname="DATA_SUB_CORRECTED"



# can manually update mslist for other selection 
#msfiles   = ascii.read('big-mslist.txt',data_start=0)
#msfiles   = list(msfiles[:][msfiles.colnames[0]]) # convert to normal list of strings

#msoutconcat   = obsid + '.dysco.sub.shift.avg.weights.ms.set0.archive'
    
if doconcat:    
    msfilesconcat = []

    #remove ms from the list where column DATA_SUB does not exist (to prevent NDPPP crash)
    for msnumber, ms in enumerate(msfiles):
      if os.path.isdir(ms):
        if mscolexist(ms,colname):
          msfilesconcat.append(ms)
        else:
          msfilesconcat.append('missing' + str(msnumber))
      else:  
        msfilesconcat.append('missing' + str(msnumber))    
        
        
    # FIRST CONCAT WITH WEIGHT_SPECTRUM_FROM_IMAGING_WEIGHT
    cmd =  'NDPPP msin="' + str(msfilesconcat) + '" msin.orderms=False '
    if aoflagger:
      if takeoutbeam:  
        cmd += 'steps=[phaseshift,aoflagger,applybeam,average] '
      else:
        cmd += 'steps=[phaseshift,aoflagger,average] '  
    else:
      if takeoutbeam:    
        cmd += 'steps=[phaseshift,applybeam,average] '
      else: 
        cmd += 'steps=[phaseshift,average] '   
    cmd += 'msin.datacolumn=%s msin.missingdata=True '%colname
    if dysco:
      cmd += 'msout.storagemanager=dysco '
    cmd += 'msin.weightcolumn=WEIGHT_SPECTRUM_FROM_IMAGING_WEIGHT '  
    cmd += 'msout=' + msoutconcat + ' '
    cmd += 'phaseshift.type=phaseshift phaseshift.phasecenter=' + phasecenter + ' '
    cmd += 'aoflagger.type=aoflagger '
    cmd += 'average.type=averager '
    cmd += 'average.timestep=' + str(timestepavg) + ' average.freqstep=' + str(freqstepavg) + ' ' 
    cmd += 'applybeam.type=applybeam applybeam.usechannelfreq=True '
    cmd += 'applybeam.beammode=array_factor ' # do no update weights from beam because in this case we just want IMAGING_WEIGHT
    
    print cmd
    os.system(cmd)



    # SECOND CONCAT WITH WEIGHT_SPECTRUM
    cmd =  'NDPPP msin="' + str(msfilesconcat) + '" msin.orderms=False '
    if aoflagger:
      if takeoutbeam:  
        cmd += 'steps=[phaseshift,aoflagger,applybeam,average] '
      else:
        cmd += 'steps=[phaseshift,aoflagger,average] '  
    else:
      if takeoutbeam:    
        cmd += 'steps=[phaseshift,applybeam,average] '
      else: 
        cmd += 'steps=[phaseshift,average] '   
    cmd += 'msin.datacolumn=%s msin.missingdata=True '%colname
    if dysco:
      cmd += 'msout.storagemanager=dysco '
 
    cmd += 'msin.weightcolumn=WEIGHT_SPECTRUM '
    cmd += 'msout=' + msoutconcat + '.tmpweight ' 
    cmd += 'phaseshift.type=phaseshift phaseshift.phasecenter=' + phasecenter + ' '
    cmd += 'aoflagger.type=aoflagger '
    cmd += 'average.type=averager '
    cmd += 'average.timestep=' + str(timestepavg) + ' average.freqstep=' + str(freqstepavg) + ' ' 
    cmd += 'applybeam.type=applybeam applybeam.usechannelfreq=True '
    cmd += 'applybeam.beammode=array_factor applybeam.updateweights=True '
    print cmd
    os.system(cmd)

    # Make a WEIGHT_SPECTRUM from WEIGHT_SPECTRUM_SOLVE
    t  = pt.table(msoutconcat, readonly=False)

    print 'Adding WEIGHT_SPECTRUM_SOLVE' 
    desc = t.getcoldesc('WEIGHT_SPECTRUM')
    desc['name']='WEIGHT_SPECTRUM_SOLVE'
    t.addcols(desc)

    t2 = pt.table(msoutconcat + '.tmpweight', readonly=True)
    imweights = t2.getcol('WEIGHT_SPECTRUM')
    t.putcol('WEIGHT_SPECTRUM_SOLVE', imweights)

    # Fill WEIGHT_SPECTRUM with WEIGHT_SPECTRUM from second ms
    t2.close()
    t.close() 

    # clean up
    os.system('rm -rf ' + msoutconcat + '.tmpweight')

    print ' '
    print ' '
    print 'Ouput column WEIGHT_SPECTRUM used for imaging (contains IMAGING_WEIGHT from DR2)'
    print 'Ouput column WEIGHT_SPECTRUM_SOLVE used for calibration (contains WEIGHT_SPECTRUM from DR2)'


    #cmd += 'msin.starttime=12May2015/19:23:22.0 msin.endtime=13May2015/01:43:00.0 '

if doflagafter:
     cmd = 'NDPPP msin=' + msoutconcat + ' msout=. msin.datacolumn=DATA ' 
     cmd += 'steps=[aoflagger,preflag] aoflagger.type=aoflagger preflag.type=preflagger '
     cmd += 'preflag.amplmax=' + str(amplmax) + ' '
     os.system(cmd)

if split:

 nchanperblock = np.int(20/freqstepavg)
 t = pt.table(msoutconcat + '/SPECTRAL_WINDOW', readonly=True)
 nchan = t.getcol('NUM_CHAN')[0]
 t.close()
 
 
 for chan in range(0,nchan,nchanperblock):
   msout = obsid + '_chan' + str(chan) + '-' + str(chan+nchanperblock-1) + '.ms'

   cmd  = 'NDPPP msin=' + msoutconcat + ' msout='+ msout + ' msin.datacolumn=DATA ' 
   if dysco:
     cmd += 'msout.storagemanager=dysco '
   cmd += 'msin.weightcolumn=WEIGHT_SPECTRUM_SOLVE '
   cmd += 'steps=[] ' + 'msin.startchan=' + str(chan) + ' '
   cmd += 'msin.nchan=' + str(nchanperblock) + ' ' + 'msout=' + msout + ' '
   print cmd
   os.system(cmd)

    
    
    
