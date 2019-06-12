#!/bin/bash

source /home/mjh/pipeline-master/init.sh
archivems.sh $WD
cd $WD
set_complete.py
