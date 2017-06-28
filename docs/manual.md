# ddf-pipeline manual

This is the users' manual for ddf-pipeline, a pipeline for the
direction-dependent self-calibration and imaging of LOFAR data.

## prerequisites

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
   
4. Source the modified `DDF.sh` file.
