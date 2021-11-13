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
dependencies, including DDFacet, wsclean and DPPP.

## system requirements

With the default configuration the extraction routines will have a
memory footprint of around 100 GB and will make good use of large
numbers of cores (because of the `DDFacet` dependency. It may be
possible to reduce the memory requirements using options that are
passed to `DDFacet`. The self-calibration routines use only a small
amount of RAM and could in principle be run on a desktop PC, although
they benefit from large numbers of cores.

## requirements for download

The extraction pipeline uses `rclone` (https://rclone.org/) to
retrieve archive files from cloud storage. You need a token (a
'macaroon') to download from the archive and you can obtain this from
the LoTSS website https://lofar-surveys.org/releases.html . Place it
in a directory and ensure that the environment variable
`RCLONE_CONFIG_DIR` points to this directory.

## running an extraction

To run an extraction call the `extract.py` script with the name of a field or object you wish to extract, and optionally the field size, right ascension and declination.

Examples:

```extract.py NGC507```

Look up the object NGC507 to determine a position and do a
default-size extraction of a region of 0.5 degrees square around the
catalogued position.

```extract.py myfield 0.5 0.4 286.1918961 59.8494461```

Extract a region of 0.4 degrees square around the specified RA and DEC, naming the working directory `myfield`.

```extract.py --chunkhours=1.0 --ncpu=8 NGC507```

As first example, but specify options to be passed to the `sub-sources-outside-region.py` script to control memory footprint.

An extraction takes a few hours on a machine with >100 GB RAM and 32 cores.

## running self-calibration

Extraction will create a directory with the downloaded pipeline output
and concatenated broad-band measurement sets, one per observation
(which may mean more than one per field).

