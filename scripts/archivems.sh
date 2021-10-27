#!/bin/bash

# Set up the LOFAR path and run the archivems python script
# This is a separate script for legacy reasons

# parameter is the directory to work in
# second (optical) parameter is column to archive

dir=${2-DATA_DI_CORRECTED}
cd ${1}
python $DDF_DIR/ddf-pipeline/scripts/archivems.py -c ${dir}
