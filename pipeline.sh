#!/bin/bash

source /home/mjh/Wirtinger_Pack2/init.sh
export PATH=/home/mjh/git/ddf-pipeline:$PATH
pipeline.py /home/mjh/git/ddf-pipeline/tier1.cfg
