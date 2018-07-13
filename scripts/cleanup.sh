#!/bin/bash

# run this if you want to do a complete restart
# saves the output of the previous run in the 'old' directory
# paths may need altering to suit user setup

mkdir old
mv image_full_ampphase1m.smooth.int.restored.fits old
mv image_full_ampphase1m.app.restored.fits old
rm image_*
rm crossmatch*
rm cube*
rm last*
rm mask*
rm external_mask*
rm -r logs
rm *.cfg
rm -r *.ddfcache
archive_old_solutions.py $DDFACET_DIR/ddf-pipeline/tier1.cfg
remove_columns.py big-mslist.txt
