# extraction and self-calibration

ddf-pipeline contains tools to allow you to extract uv data for a
small region of the LoTSS sky and self-calibrate using either your own
tools or a self-calibration script using wsclean that we provide. The
extraction and self-calibration process is described by van Weeren et
al (2021):
https://ui.adsabs.harvard.edu/abs/2021A%26A...651A.115V/abstract . If
you use these routines please cite that paper as well as the relevant
LoTSS survey papers.

The singularity image is the recommended way of installing the
dependencies, including DDFacet, wsclean and DP3.

## system requirements

With the default configuration the extraction routines will have a
memory footprint of around 100 GB and will make good use of large
numbers of cores (because of the `DDFacet` dependency). It may be
possible to reduce the memory requirements using options that are
passed to `DDFacet`. The self-calibration routines use only a small
amount of RAM and could in principle be run on a desktop PC, although
they benefit from large numbers of cores.

## requirements for download

The extraction pipeline will download from several possible locations. DR2 files are currently located in the SURF Science Data Repository and are downloadable without credentials using an API provided by the SDR. Other files (non-DR2 fields) are located in cloud storage and the extraction pipeline uses `rclone` (https://rclone.org/) to
retrieve archive files from cloud storage. You need a special token (a
'macaroon') to download from the archive, again only available to LoTSS collaborators. If using this, place the macaroon in a directory and
ensure that the environment variable `RCLONE_CONFIG_DIR` points to
this directory.

## running an extraction

To run an extraction call the `extraction.py` script with the name of a field or object you wish to extract, and optionally the field size, right ascension and declination.

Examples:

```extraction.py NGC507```

Look up the object NGC507 to determine a position and do a
default-size extraction of a region of 0.5 degrees square around the
catalogued position.

```extraction.py myfield 0.4 286.1918961 59.8494461```

Extract a region of 0.4 degrees square around the specified RA and DEC, naming the working directory `myfield`.

```extraction.py --chunkhours=1.0 --ncpu=8 NGC507```

As first example, but specify options to be passed to the `sub-sources-outside-region.py` script to control memory footprint.

An extraction takes a few hours on a machine with >100 GB RAM and 32 cores.

## running self-calibration

Extraction will create a directory with the downloaded pipeline output
and concatenated broad-band measurement sets, one per observation
(which may mean more than one per field). These measurements sets will be in per-observation directories with names `*.dysco.sub.shift.avg.weights.ms.archive?`.

Inside the working directory for the target, created by the extraction, do

```
runwsclean.py -b NAME.ds9.reg -i NAME *.dysco.sub.shift.avg.weights.ms.archive?
```

where the `NAME` is the name of your target. This assumes that all
measurement sets extracted are equally good. The self-calibration
script will go through 9 iterations of phase and then amplitude and
phase self-calibration, taking between hours and days depending on the
number of measurement sets, the size and complexity of the field, and
the speed of your machine and underlying file system. You can monitor
its progress by looking at PNG files in the working directory
(`NAME_*.png` and `plotlosoto*/*.png`). If a particular extracted
measurement set is bad you may want to consider excluding it from the
self-calibration. The final output image will be stored in your
working directory as `NAME_9-MFS.image.fits`.

## flux scale

The bootstrap flux scale corrections applied by ddf-pipeline and
described by Shimwell et al (2019) are applied to the input
data. However, you are advised to check the flux scale by comparison
with the released DR2 mosaics (where an additional correction has been
applied) or with a known flux density for your target or other objects
in the field.
