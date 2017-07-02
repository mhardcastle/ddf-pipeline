#!/bin/bash

CFG=${CONFIG:-/home/mjh/git/public/ddf-pipeline/examples/tier1.cfg}
echo Using ddf-pipeline config $CFG
source /home/mjh/git/public/DDF.sh
pipeline.py $CFG
