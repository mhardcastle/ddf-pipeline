#!/bin/bash

export WD=`pwd`
sed -e "s|INSTALLDIR|$WD|" ddf-pipeline/misc/DDF.sh > ddf-pipeline/init.sh
