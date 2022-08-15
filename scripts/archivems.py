from __future__ import print_function
# MS archiving script contributed by R. van Weeren
# LOFAR software must be on path and dysco on LD_LIBRARY_PATH

import os, sys
import glob
import argparse
import casacore.tables as pt
from casacore.tables.tableutil import makescacoldesc, makearrcoldesc, \
    makecoldesc, maketabdesc
import numpy as np

def make_imaging_weight_column(ms):
    
    t = pt.table(ms, readonly=False, ack=False)
    cnames = t.colnames()
  
    try:
        cdesc = t.getcoldesc('DATA')
    except:
        raise ValueError('Column DATA does not exist')
    # Determine if the DATA storage specification is tiled.
    hasTiled = False
    try:
        dminfo = t.getdminfo("DATA")
        if dminfo['TYPE'][:5] == 'Tiled':
            hasTiled = True
    except:
        hasTiled = False
    # Use TiledShapeStMan if needed.
    if not hasTiled:
        dminfo = {'TYPE': 'TiledShapeStMan',
                  'SPEC': {'DEFAULTTILESHAPE': [4, 32, 128]}}    

    if 'IMAGING_WEIGHT' in cnames:
        six.print_("Column IMAGING_WEIGHT not added; it already exists")
    else:
        # Add IMAGING_WEIGHT which is 1-dim and has type float.
        # It needs a shape, otherwise the CASA imager complains.
        shp = []
        if 'shape' in cdesc:
            shp = cdesc['shape']
        if len(shp) > 0:
            shp = [shp[0]]  # use nchan from shape
        else:
            shp = [t.getcell('DATA', 0).shape[0]]  # use nchan from actual data
        cd = makearrcoldesc('IMAGING_WEIGHT', 0, ndim=1, shape=shp,
                            valuetype='float')
        dminfo = {'TYPE': 'TiledShapeStMan',
                  'SPEC': {'DEFAULTTILESHAPE': [32, 128]}}
        dminfo['NAME'] = 'imagingweight'
        t.addcols(maketabdesc(cd), dminfo)

        print("added column IMAGING_WEIGHT")
     
    t.flush()
    t.close()

    return


parser = argparse.ArgumentParser(description='Compress and copy column of a list of ms files for archiving')
parser.add_argument('-m','--mslist', help='DR2 mslist file, default=big-mslist.txt', default='big-mslist.txt', type=str)
parser.add_argument('-c','--column', help='Column that is copied from the MS and compressed, default=DATA_DI_CORRECTED', default='DATA_DI_CORRECTED', type=str)
parser.add_argument('--inmsdir', help='Forces the ouput MS to be in the same directory as the mslist', action='store_true')
parser.add_argument('--skipimweights', help='Do not copy over IMAGING_WEIGHT column', action='store_true')
parser.add_argument('--preserve', help='Preserve existing files', action='store_true')
args = vars(parser.parse_args())

dirname = os.path.dirname(args['mslist'])
msfiles = [l.rstrip() for l in open(args['mslist']).readlines()]

for ms in msfiles:
    
    if dirname !='': # user gave a mslist that is not in the current directory
      msin = dirname + '/' + ms # add the full path to the input file
    else:
      msin = ms  
  
    if args['inmsdir']:
      msout = dirname + '/' + ms + '.archive'
    else:
      msout = ms + '.archive' # write ms in current working directory

    if not args['preserve']:
        cmd  = 'rm -r '+msout+'; '
    else:
        cmd = ''
    cmd += 'DP3 msin=' + msin + ' msin.datacolumn=' + args['column'] + ' '
    cmd += 'msout.storagemanager=dysco msout=' + msout  + ' steps=[] '
    cmd += 'msin.weightcolumn=WEIGHT_SPECTRUM '

    if not os.path.isdir(msout):
      print(cmd)  
      result=os.system(cmd)
      if result!=0:
          os.system('rm -r '+msout)
          raise RuntimeError('DP3 call failed')
      
      if not args['skipimweights']:
        tin = pt.table(msin, ack=False)
        if 'IMAGING_WEIGHT' in tin.colnames():
           make_imaging_weight_column(msout)
           tout = pt.table(msout, readonly=False, ack=False)
           imw = tin.getcol('IMAGING_WEIGHT')
           tout.putcol('IMAGING_WEIGHT', imw)
           tout.close()
        else:
           print('Warning: IMAGING_WEIGHT does not exist in input ms', ms)
      
        tin.close()
    else:
        print('Skipping',msout,'as it already exists')
        
