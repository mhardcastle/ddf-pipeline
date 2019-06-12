#!/bin/bash

CFG=${CONFIG:-/home/mjh/pipeline-master/ddf-pipeline/examples/tier1-jul2018.cfg}
echo Using ddf-pipeline config $CFG
source /home/mjh/pipeline-master/init.sh

mkdir KEEP
mv image_full_ampphase_di_m.NS* KEEP
mv image_full_low.* image_full_low_im.* image_full_low_m.* KEEP
mv full-mask-low.fits external_mask_ext-deep.fits KEEP
mv facet-offset.txt KEEP
mv panstarrs-facet_offsets.npy KEEP
mv DynSpec* KEEP
rm -r *.ddfcache

pipeline.py $CFG
