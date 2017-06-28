# ddf-pipeline manual

This is the users' manual for ddf-pipeline, a pipeline for the
direction-dependent self-calibration and imaging of LOFAR data.

## prerequisites

In addition to the software directly required for ddf-pipeline, you
will need:

* pyrap
* pybdsf (possibly from the LOFAR tree)
* astropy
* pyregion
* emcee
* reproject (for mosaicing)

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

