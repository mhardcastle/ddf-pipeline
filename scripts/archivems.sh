#!/bin/bash

# Set up the LOFAR path and run the archivems python script
# We wrap the script up to avoid having to have the LOFAR software
# on path for all the pipeline operations

# parameter is the directory to work in

eval /usr/bin/modulecmd bash load gcc-6.4
source /soft/lofar-071018/init.sh

cd $1
python $DDF_DIR/ddf-pipeline/scripts/archivems.py
