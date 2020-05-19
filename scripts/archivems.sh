#!/bin/bash

# Set up the LOFAR path and run the archivems python script
# We wrap the script up to avoid having to have the LOFAR software
# on path for all the pipeline operations

# parameter is the directory to work in
# second (optical) parameter is column to archive

source $LOFARSOFT

dir=${2-DATA_DI_CORRECTED}
cd ${1}
python $DDF_DIR/ddf-pipeline/scripts/archivems.py -c ${dir}
