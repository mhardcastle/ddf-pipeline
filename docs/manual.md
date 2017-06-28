# ddf-pipeline manual

This is the users' manual for ddf-pipeline, a pipeline for the
direction-dependent self-calibration and imaging of LOFAR data.

## what it does

ddf-pipeline carries out several iterations of direction-dependent
self-calibration on your LOFAR imaging data using DDFacet and KillMS
For more on KillMS see Tasse 2014a,b
(http://adsabs.harvard.edu/abs/2014arXiv1410.8706T
http://adsabs.harvard.edu/abs/2014A%26A...566A.127T) and Smirnov &
Tasse 2015 (http://adsabs.harvard.edu/abs/2015MNRAS.449.2668S); for
DDFacet see Tasse et al 2017 submitted. The objective of ddf-pipeline is to
fully reduce LOFAR (continuum, Stokes I) imaging data to a
science-ready state with no human intervention. ddf-pipeline is
currently in use by the LOFAR Surveys KSP for reduction of the LoTSS
survey (Shimwell et al 2017
http://adsabs.harvard.edu/abs/2017A%26A...598A.104S).

ddf-pipeline was written by Martin Hardcastle, Tim Shimwell, Cyril
Tasse and Wendy Williams and will be described by Shimwell et al (in prep).
Scientific users of ddf-pipeline are requested to cite the relevant
papers and refer to the ddf-pipeline github.

## hardware prerequisites

For reduction of LOFAR data you will need a Linux machine with at
least 192 GB of physical RAM. 256 GB RAM and 32 cores works well; this
will give a full reduction in 3-4 days. Hyperthreading does not give a
significant advantage and can be disabled.

Several TB of fast storage attached (ideally) directly to the node are
necessary for a reduction of a 48-MHz 8-h dataset.

As root you or or your sysadmin should do:

```
mount -o remount -o size=256g /dev/shm
echo never > /sys/kernel/mm/transparent_hugepage/defrag
```

Also ensure that your machine has some swap &mdash; 100 GB or so is fine.

## software prerequisites

In addition to the software directly required for ddf-pipeline (see below), you
will need:

* pyrap
* pybdsf (possibly from the LOFAR tree)
* astropy
* pyregion
* emcee
* reproject (for mosaicing)

Recent versions of numpy and numexpr are important for DDFacet and KillMS.

## installation

Currently we recommend that you install ddf-pipeline and its
prerequisites KillMS, DDFacet and SkyModel directly from the Github
repositories, so that you can easily get updates with `git pull`.

1. Make a suitable directory and clone all of the Github repos to it:

```
mkdir DDF
cd DDF
git clone https://github.com/mhardcastle/ddf-pipeline.git
git clone https://github.com/cyriltasse/killMS2.git
git clone https://github.com/cyriltasse/SkyModel.git
git clone https://github.com/cyriltasse/DDFacet.git
```

2. Compile the DDFacet / killMS code by following the instructions in
   the relevant documentation.
   
3. Copy the `DDF.sh` script from this directory to the root directory
   of your working directory and modify it appropriately to refer to
   your installation directory.
   
4. Source the modified `DDF.sh` file. The ddf-pipeline scripts
   directory, DDFacet and KillMS should now all be on your path. We
   assume that you have done this from here on.

## directory structure

ddf-pipeline provides the following directory structure:

* scripts: should be on your PATH and PYTHONPATH
* utils: should be on your PYTHONPATH
* examples: contains example configuration files
* docs: contains documentation, including this file
* misc: miscellaneous useful files

## getting ready to run

ddf-pipeline assumes that your LOFAR data have been through prefactor
(https://github.com/lofar-astron/prefactor) or an equivalent
(e.g. Reinout van Weeren's equivalent scripts). If you are using LOFAR
surveys KSP data, you may be able to persuade someone to run this for
you on SARA.

To image a full LOFAR field we recommend averaging to 2 channels per
subband and 8s in time.

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

(The script `run_pipeline.py` packages up the downloading, unpacking
and making of the measurement sets. It will need modifying for your
local environment before use and so we do not describe it here.)

## preparing the configuration file

`pipeline.py` takes as command line options an arbitrary mixture of
names of configuration files (more than one can be specified) and
direct option settings. If you run the script without any options, the
default settings are shown. Most of these are sensible for LOFAR data.

Config files follow the standard Python ConfigParser layout of section
names in square brackets followed by name-value pairs separated by an
equals sign.  An example configuration file is in examples/example.cfg
. A more detailed example is in examples/tier1.cfg &mdash; note that
this includes many settings which are the current defaults.

Command-line options are specified by two dashes, the name of the
section, a dash, and the name of the parameter within that section,
followed by an equals sign and the parameter value. They over-ride
settings in files, so you can do e.g.

```
pipeline.py examples/tier1.cfg --control-restart=False
```

Below we describe each section of the config file in more detail.

### [control]

Options that control the overall running of the pipeline. Defaults are
mostly sensible here. If you expect to have to restart frequently, set
`clearcache=False` otherwise a lot of time will be spent re-making the
cache.

Some parts of the pipeline are only enabled by switching them on
here, for example `bootstrap=True` enables bootstrap (see below)

### [bootstrap]

Specifies the images and catalogues to be used for bootstrap. See
below for more information on this.

### [data]

Specifies the MS lists to be used. This must be set. `mslist` refers
to your short MS list and `full_mslist` to the MS list containing all
the data. `colname` is the column of the MS to image (usually
CORRECTED_DATA).

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
that all bright TGSS sources are automatically included, which helps
the self-calibration to converge faster. If a TGSS path is specified,
the other default options are probably sensible. Enable the use of masks for extended sources with `tgss_extended=True`.

### [offsets]

Controls whether offsets from the optical frame are corrected for (see
below).

### [solutions]

These control KillMS (most directly translate into KillMS input
parameters). Generally these values should be sensible. An important
choice is whether to set a `uvmin` (in km) for self-calibration. The
default (no uvmin) gives maps with good noise properties but tends to
absorb unmodelled extended structure. We currently use `uvmin=1.5`
which does a better job of preserving large-scale extended flux at the
cost of some additional structure in the large-scale noise. Most other
options should be left at their default settings.

## bootstrap

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

## offsets
Set
```
[offsets]
method=panstarrs
fit=mcmc
```
to determine and correct for the per-direction offset from the
PanSTARRS frame. This is particularly useful if you intend to do
optical crossmatching or any form of mosaicing.

## restart

The option `[control] redofrom` can be used to start again from after
a specified step in the pipeline. Available options are start, dirin,
phase, ampphase, fullow, full. Solutions and, if necessary, bootstrap
results will be removed before the restart; all old files are moved to
an `archive_dir` directory. `[control] restart` must be True for this
to work.

## signal handling

By default, the pipeline interprets SIGUSR1 as a request to stop
running at the next convenient point, i.e. normally when a KillMS or
DDF would otherwise be about to start. Initiate this process with
e.g. `pkill -USR1 -f pipeline.py`.

## quality pipeline

Once the pipeline has run you can use `quality_pipeline.py` to get
some basic quality information, including estimates of positional and
flux scale offsets. Copy `quality_pipeline.cfg` from the examples
directory and amend suitably. By default this looks for the standard
pipeline products in the current working directory. Catalogues used in
the quality pipeline must conform (at least roughly) to the PyBDSM
format. See http://www.extragalactic.info/bootstrap/ for examples.

TGSS and FIRST catalogues must be provided here to allow the code to run.

## mosaicing

Make a mosaic of adjacent observations from the pipeline with the
`mosaic.py` script. At a minimum it needs directories specified with
the `--directories` argument which contain standard pipeline
output. You can also determine the noise from the image with the
`--find_noise` option and remove global offsets with respect to FIRST
(quality pipeline must have been run first) with the `--shift` option.

