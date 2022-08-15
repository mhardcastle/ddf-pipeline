# ddf-pipeline manual

This is the users' manual for ddf-pipeline, a pipeline for the
direction-dependent self-calibration and imaging of LOFAR data.

## introduction

ddf-pipeline carries out several iterations of direction-dependent
self-calibration on your LOFAR imaging data using DDFacet and KillMS
For more on KillMS see Tasse 2014a,b
(http://adsabs.harvard.edu/abs/2014arXiv1410.8706T
http://adsabs.harvard.edu/abs/2014A%26A...566A.127T) and Smirnov &
Tasse 2015 (http://adsabs.harvard.edu/abs/2015MNRAS.449.2668S); for
DDFacet see Tasse et al 2017 (https://arxiv.org/abs/1712.02078). The
objective of ddf-pipeline is to fully reduce LOFAR (continuum, Stokes
I) imaging data to a science-ready state with no human
intervention. ddf-pipeline is currently in use by the LOFAR Surveys
KSP for reduction of the LoTSS survey (Shimwell et al 2017
http://adsabs.harvard.edu/abs/2017A%26A...598A.104S). Version 2 of
ddf-pipeline is the version used for the second data release (DR2) of
LoTSS.

ddf-pipeline was written by Martin Hardcastle, Tim Shimwell, Cyril
Tasse and Wendy Williams and is described by Shimwell et al 2019
(https://www.aanda.org/articles/aa/full_html/2019/02/aa33559-18/aa33559-18.html)
and Tasse et al 2021
(https://ui.adsabs.harvard.edu/abs/2021A%26A...648A...1T/abstract).
Scientific users of ddf-pipeline are requested to cite the relevant
papers and refer to the ddf-pipeline github.

## who is this release for?

This release of ddf-pipeline is aimed at people who want to
reduce wide-field LOFAR data quickly and to a good standard. If you
are interested in a single bright source in a LOFAR field and are not
in a hurry, you may find factor
(https://github.com/lofar-astron/factor) to be more useful to you.

## getting support

Support for ddf-pipeline and the code that backs it up is provided on
a best-efforts basis &mdash; it is not supported by ASTRON and all of
the programmers have other work to do. A certain level of technical
expertise is required to run ddf-pipeline at all.

Please request support by raising an issue on the github, **not** by
direct e-mail to the programmers.

## hardware prerequisites

For reduction of LOFAR data you will need a Linux machine with at
least 192 GB of physical RAM. 256 GB RAM and 32 cores works well; this
will give a full reduction in 4-5 days. Hyperthreading does not give a
significant advantage and can be disabled.

Several TB of fast storage attached (ideally) directly to the node are
necessary for a reduction of a 48-MHz 8-h dataset.

As root you or your sysadmin should do:

```
mount -o remount -o size=256g /dev/shm
echo never > /sys/kernel/mm/transparent_hugepage/defrag
```

Also ensure that your machine has some swap &mdash; 100 GB or so is fine.

## software prerequisites

In addition to the software directly required for ddf-pipeline (see below), you
will need:

* git (for installation)
* pyrap
* pybdsf (possibly from the LOFAR tree but an up-to-date separate build is recommended)
* astropy
* astropy-healpix
* pyregion
* emcee
* reproject (for mosaicing)
* Polygon
* pyephem
* pyFFTW
* SharedArray
* deap
* py-cpuinfo
* cfitsio fpack (only if compression required: fpack must be on PATH)
* DP3 and dysco (for MS compression)

Recent versions of numpy and numexpr are important for DDFacet and KillMS.

If you are doing HBA surveys processing and wish to use the database
integration features you must ensure that `surveys_db.py` from
https://github.com/mhardcastle/lotss-query is on the PYTHONPATH. This
introduces dependencies on sshtunnel and MySQLdb.

Infrastructure for the LoTSS HBA survey is (mostly) in https://github.com/mhardcastle/lotss-hba-survey .

## installation

Given that the whole ddf-pipeline infrastructure has many
dependencies, we include a singularity recipe `ddf-py3.singularity` in
the repository which you may use to build a singularity
image. Alternatively, following the steps in the post-install part of
this script will allow you to install from scratch with root access on
a Debian bullseye system (recent Ubuntu versions will work similarly.)

If you are manually installing, once ddf-pipeline is installed in a
source directory you need to source ddf-pipeline/init.sh to get the
pipeline and other utilities on your `PATH`/`PYTHONPATH`. Make sure
your `LD_LIBRARY_PATH` points to the installation location of
libraries such as Dysco.

## directory structure

ddf-pipeline provides the following directory structure:

* scripts: main scripts, should be on your PATH and PYTHONPATH after `init.sh`
* utils: utility code, should be on your PYTHONPATH after `init.sh`
* examples: contains example configuration files
* docs: contains documentation, including this file
* misc: miscellaneous useful files

## what it does

For a full description of version 2 of ddf-pipeline see Tasse et al 2021 (https://ui.adsabs.harvard.edu/abs/2021A%26A...648A...1T/abstract).

## getting ready to run

ddf-pipeline assumes that your LOFAR data have been through prefactor
(https://github.com/lofar-astron/prefactor). If you are using LOFAR
surveys KSP data or co-observing with the KSP, you may be able to
persuade someone to run this for you on SARA.

To image a full LOFAR field we recommend averaging to 2 channels per
subband and 8s in time. 1 channel/16 seconds will take about a factor
2 longer. International baselines should not be present.

If the data have been processed for you by the KSP, then you should be
able to run `download.py L123456` (where here, and from now on,
L123456 is the LOFAR observation ID of your data). This will create a
directory of the same name and fill it with tar files. Work in this
directory from now on, and unpack the tar files with
`unpack.py`. Otherwise, create a working directory, and move or copy
the prefactor-processed measurement sets there.

You now need to make two measurement set lists. One small MS list is
used for the initial calibration and imaging: the full list, which
contains all your measurement sets is used only when a good sky model
has been developed. If you are using KSP data, then `make_mslists.py`
will build these two MS lists for you.

## preparing the configuration file

`pipeline.py` takes as command line options an arbitrary mixture of
names of configuration files (more than one can be specified) and
direct option settings. If you run the script without any options, the
default settings are shown. Most of these are sensible for LOFAR data.

Config files follow the standard Python ConfigParser layout of section
names in square brackets followed by name-value pairs separated by an
equals sign.  An example configuration file is in examples/example.cfg
. A more detailed example is in examples/tier1-jul2018.cfg &mdash; note that
this includes many settings which are the current defaults.

Command-line options are specified by two dashes, the name of the
section, a dash, and the name of the parameter within that section,
followed by an equals sign and the parameter value. They over-ride
settings in files, so when running the code you can do e.g.

```
pipeline.py examples/tier1.cfg --control-restart=False
```

Below we describe each section of the config file in more detail.

### [control]

Options that control the overall running of the pipeline. Defaults are
mostly sensible here. If you expect to have to restart frequently, set
`clearcache=False` otherwise a lot of time will be spent re-making the
cache. Note that by default `clearcache_end=True` and so DDF cache
will be tidied up at the end of a successful run.

Some parts of the pipeline are only enabled by switching them on
here, for example `bootstrap=True` enables bootstrap (see below)

If you have your MS files on slow/high-latency remote storage, you may
want to set `cache_dir` to any fast local storage you have (order of 1
TB may be used during the pipeline run).

### [bootstrap]

Specifies the images and catalogues to be used for bootstrap. See
below for more information on this.

### [data]

Specifies the MS lists to be used. This must be set. `mslist` refers
to your short MS list and `full_mslist` to the MS list containing all
the data. If no `full_mslist` is set then only the steps involving the
short MS list will be carried out. `colname` is the column of the MS
to image (usually CORRECTED_DATA).

### [image]

The properties of the image to be made. These parameters are generally
sensible for high-declination LOFAR data. Note that there are
different values of the robust parameter and the PSF size for the
intermediate (self-cal) maps and the final maps made with the full
bandwidth. If you want to make a low-resolution image from the final
self-calibrated data, set `low_imsize`, `low_cell` and `low_psf_arcsec`.

### [machine]

Parameters of the machine to be used, particularly the number of
CPUs. Normally these should be calculated automatically by
ddf-pipeline.

### [masking]

This section allows control over the masks made for cleaning. A useful
option is to use the TGSS ADR catalogue, by setting `tgss=TGSS.fits`
where TGSS.fits is the filename of the full FITS catalogue downloaded
from http://tgssadr.strw.leidenuniv.nl/doku.php (see Intema et al 2016
https://arxiv.org/abs/1603.04368 for more on TGSS). This will ensure
that all bright TGSS sources are automatically included in the mask,
which helps the self-calibration to converge faster. If a TGSS path is
specified, the other default options are probably sensible. Enable the
use of extended masks for extended sources with `tgss_extended=True`
&mdash; this is slightly more risky but can be helpful if you have a bright
extended source far off-axis in your field.

### [offsets]

Controls whether offsets from the optical frame are corrected for (see
below).

### [solutions]

These control KillMS (most directly translate into KillMS input
parameters). Generally the default values should be sensible.

### [compression]

Controls the compression of some final images and measurement sets. If
`compress_ms` is `True` then DP3 will be run to compress the output data.

## running the code

Once you have set up a config file, make sure you are in the directory
where you want to work (and probably where your MS are stored) and do

```
pipeline.py my_config_file.cfg
```

The code will then start the self-calibration cycle. Come back in 3-4 days.

## the output

The code produces a large amount of output, most of it no use to
you. Good images to look at (use a FITS viewer like ds9) are files of
the form `image_*_app.restored.fits`. These are the final imaging
products, in apparent flux units, for each of the self-calibration
steps. When the code is finished,
`image_full_ampphase1m.int.restored.fits` is the full-resolution image
in physical units (unless you have specified correction for offsets,
see below). A file `summary.txt` is also created summarizing the
run. Log files for the individual steps are stored by default in the
'logs' directory &mdash; you can use the `analyse_logs.py` script to
do some very basic profiling of the time taken by the various steps.

## advanced topics

### bootstrap

Flux scale bootstrap (see Hardcastle et al 2016 http://adsabs.harvard.edu/abs/2016MNRAS.462.1910H) requires the following in the config file:

```
[control]
bootstrap=True

[bootstrap]
catalogues=['cat1.fits','cat2.fits'...]
radii=[40,10..]
frequencies=[74e6,327e6,...]
```

Catalogues, radii, frequencies are lists in Python list
format.

* catalogues must be a list of existing files in PyBDSM format
(i.e. they must contain RA, Dec, Total_flux, E_Total_flux). Some suitable
files are available at http://www.extragalactic.info/bootstrap/ . 

* radii are the separation distances in arcsec to use for each catalogue. If not specified this will default to 10 arcsec for each, but this is probably not what you want.

* frequencies are the frequencies of each catalogue in Hz.

Optionally you may supply in the `bootstrap` block a list of names which will be used in the
crossmatch tables to identify the catalogues (`names`) and a
list of group identifiers for each catalogue (`groups`). If `groups` is used, matches from one or more of the catalogues in each group are required for a source to be fitted. The default is that a match in each individual catalogue is required.

Bootstrap operates on the `mslist` specified and so needs that to contain enough measurement sets to do the fitting. 5 or 6 is a good compromise.

If bootstrap runs successfully then a `SCALED_DATA` column will be
generated and used thereafter for imaging.

Results can be plotted using the `plot_factors.py` script.

Bootstrap improves the flux scale, but testing during DR2 suggests
that there is still a reasonable amount of flux scale scatter after
the process, probably at least 10 per cent (see Shimwell et al in
prep.) and maybe more on a per field basis. Ideally you should check
the quality of your flux scale by reference to some known object in
the field, or several of them.

### offsets

ddf-pipeline can correct for astrometric offsets using a version of
the method of Smith et al (2011) &mdash;
http://adsabs.harvard.edu/abs/2011MNRAS.416..857S . Empirically we
have found that the best results come from aligning with PanSTARRS,
which works anywhere in the Northern hemisphere. The correction method 'pslocal' uses a version of the PanSTARRS catalogue optimized for this and made available through the Hertfordshire cluster, and is recommended.

Set
```
[offsets]
method=pslocal
fit=mcmc
```
to determine and correct for the per-direction offset from the
PanSTARRS frame. This is particularly useful if you intend to do
optical crossmatching or any form of mosaicing. Output images after
shifting have names of the form
`image_full_ampphase1m_shift.*.facetRestored.fit`

### restart

If the pipeline crashes, then if `[control] restart=True`, the
default, it will attempt to pick up from where it left off on
rerun. Normally this is what you want.

The option `[control] redofrom` can be used to start again from after
a specified step in the pipeline. Available options are `start` and
`dirin`.  `start` will clear all but the input MSs, `dirin` will
retain the initial direction-independent images.  All old files are
moved to an `archive_dir` directory. `[control] restart` must be True
for this to work.

### signal handling

By default, the pipeline interprets SIGUSR1 as a request to stop
running at the next convenient point, i.e. normally when a KillMS or
DDF would otherwise be about to start. Initiate this process with
e.g. `pkill -USR1 -f pipeline.py`.

### quality pipeline

Once the pipeline has run you can use `quality_pipeline.py` to get
some basic quality information, including estimates of positional and
flux scale offsets. Copy `quality_pipeline.cfg` from the examples
directory and amend suitably. By default this looks for the standard
pipeline products in the current working directory. Catalogues used in
the quality pipeline must conform (at least roughly) to the PyBDSM
format. See http://www.extragalactic.info/bootstrap/ for examples.

TGSS and FIRST catalogues must be provided here to allow the code to run.

### mosaicing

Make a mosaic of adjacent observations from the pipeline with the
`mosaic.py` script. At a minimum it needs directories specified with
the `--directories` argument which contain standard pipeline
output. You can also determine the noise from the image with the
`--find_noise` option.

