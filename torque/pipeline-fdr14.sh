#!/bin/bash

CFG=${CONFIG:-/home/mjh/pipeline-master/ddf-pipeline/examples/tier1-jul2018.cfg}
echo Using ddf-pipeline config $CFG
source /home/mjh/pipeline-master/init.sh
source $LOFARSOFT
pipeline.py $CFG
