#!/bin/bash

CFG=${CONFIG:-/home/mjh/pipeline-master/ddf-pipeline/examples/quality-example.cfg}
echo Using quality config $CFG
source /home/mjh/pipeline-master/init.sh

quality_pipeline.py $CFG
