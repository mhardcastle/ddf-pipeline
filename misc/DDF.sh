#!/bin/bash

### set DDF_DIR to the installation directory
export DDF_DIR=INSTALLDIR
export KILLMS_DIR=$DDF_DIR
### don't modify anything below here

#echo -e Source directory for DDF software: $DDF_DIR

export PYTHONPATH=$DDF_DIR/ddf-pipeline/scripts:$DDF_DIR/ddf-pipeline/utils:$PYTHONPATH

export PATH=$DDF_DIR/DynSpecMS:$PATH
export PATH=$DDF_DIR/ddf-pipeline/scripts:$PATH
export PATH=$DDF_DIR/lotss-hba-survey:$PATH
export PATH=$DDF_DIR/killMS/killMS:$PATH

export PYTHONPATH=$DDF_DIR/lotss-query:$PYTHONPATH
export PYTHONPATH=$DDF_DIR/lotss-hba-survey:$PYTHONPATH
export PYTHONPATH=$DDF_DIR:$PYTHONPATH

# force numpy to stay single-threaded -- maybe
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
