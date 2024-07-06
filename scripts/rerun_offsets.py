#!/usr/bin/env python
# Rerun the offset fitting in the working directory. This will re-make
# the pslocal files using the code from pipeline.py and then re-run
# the offsets using offsets.py.

from __future__ import print_function
from get_cat import get_cat,download_required
from offsets import do_offsets

method='pslocal'
if download_required(method):
    get_cat(method)

do_offsets({'method':method,'mode':'normal','fit':'mcmc','cellsize':1.5,'no_astrometry':True})
