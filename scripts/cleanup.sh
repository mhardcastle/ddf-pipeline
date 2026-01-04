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
rm -r *.ddfcache
rm -r SOLSDIR
rm DD*.npz
rm MaskDiffuse.pickle Noise*fits MaskDiffuse.reg MaskDiffuse.fits
rm Predict_*
rm _has*
rm *.last
rm local_*.txt
rm pop.myPickle
# commented this as it was removing the parset
# rm *.cfg
# archive_old_solutions.py $DDFACET_DIR/ddf-pipeline/tier1.cfg
remove_columns.py big-mslist.txt
CleanSHM.py
#rm -rf /tmp/*
