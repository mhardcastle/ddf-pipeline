#!/bin/bash

source /home/mjh/pipeline-py3/init.sh
cd /beegfs/car/mjh/vlow
make_vlow_stokesI.py --field $FIELD
