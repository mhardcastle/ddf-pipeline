#!/bin/bash

### change these lines to suit your local setup

### set DDF_DIR to the installation directory
export DDF_DIR=/home/mjh/git/public

### a minimal PYTHONPATH is desirable but it must include all of the
### prerequisites. Only add the things you need
export PYTHONPATH=/home/mjh/python_modules/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/usr/lib/python2.7/site-packages
export PYTHONPATH=/soft/lofar-130916/lib64/python2.7/site-packages:$PYTHONPATH
export PYTHONPATH=/soft/pyrap-220216/usr/lib64/python2.7/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=/soft/lofar-130916/lib64
export LD_LIBRARY_PATH=/soft/casacore-220216/lib:$LD_LIBRARY_PATH
### don't change anything below here

echo -e Source directory for DDF software: $DDF_DIR

export KILLMS_DIR=$DDF_DIR
export PYTHONPATH=$DDF_DIR/ddf-pipeline/scripts:$DDF_DIR/ddf-pipeline/utils:$PYTHONPATH
export PYTHONPATH=$DDF_DIR:$PYTHONPATH
export PYTHONPATH=$DDF_DIR/DDFacet:$PYTHONPATH

export LD_LIBRARY_PATH=$DDF_DIR/DDFacet/DDFacet/cbuild:$LD_LIBRARY_PATH

export PATH=$DDF_DIR/killMS2:$PATH
export PATH=$DDF_DIR/DDFacet/DDFacet:$PATH
export PATH=$DDF_DIR/SkyModel:$PATH
export PATH=$DDF_DIR/ddf-pipeline/scripts:$PATH

# force numpy to stay single-threaded -- maybe
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
