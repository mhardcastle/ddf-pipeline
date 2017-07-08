#!/bin/bash

CFG=${CONFIG:-/home/mjh/temp/public/ddf-pipeline/examples/tier1.cfg}
echo Using ddf-pipeline config $CFG
source /home/mjh/temp/public/init.sh
pipeline.py $CFG
