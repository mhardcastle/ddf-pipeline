# ddf-pipeline

Prerequisites: DDFacet and KillMS (duh). astropy for use of TGSS in
mask-making.

emcee and pyrap must be on path for bootstrap, and
mpi4py if you want to use MPI to speed up emcee.

The bootstrap code also expects the script directory to be on user's
PATH.

Most options are configurable with the config file -- see options.py
for possible options and their defaults.
