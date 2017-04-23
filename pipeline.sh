#!/bin/bash

CFG=${CONFIG:-/home/mjh/git/ddf-pipeline/tier1.cfg}
echo Using ddf-pipeline config $CFG
source /home/mjh/Wirtinger_amSNR/init.sh
export PATH=/home/mjh/git/ddf-pipeline:$PATH
pipeline.py $CFG
