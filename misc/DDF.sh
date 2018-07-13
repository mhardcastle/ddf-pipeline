#!/bin/bash

### set DDF_DIR to the installation directory
export DDF_DIR=INSTALLDIR
### don't modify anything below here

echo -e Source directory for DDF software: $DDF_DIR

export KILLMS_DIR=$DDF_DIR
export PYTHONPATH=$DDF_DIR/ddf-pipeline/scripts:$DDF_DIR/ddf-pipeline/utils:$PYTHONPATH
export PYTHONPATH=$DDF_DIR:$PYTHONPATH
export PYTHONPATH=$DDF_DIR/DDFacet/build/lib:$DDF_DIR/DDFacet/build/lib/DDFacet:$PYTHONPATH

export LD_LIBRARY_PATH=$DDF_DIR/DDFacet/DDFacet/cbuild:$LD_LIBRARY_PATH

export PATH=$DDF_DIR/killMS:$PATH
export PATH=$DDF_DIR/DDFacet/build/scripts-2.7:$PATH
export PATH=$DDF_DIR/SkyModel:$PATH
export PATH=$DDF_DIR/DynSpecMS:$PATH
export PATH=$DDF_DIR/ddf-pipeline/scripts:$PATH

# force numpy to stay single-threaded -- maybe
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
