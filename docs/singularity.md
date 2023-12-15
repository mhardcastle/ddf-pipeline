# building the singularity image

This is the recommended way of installing ddf-pipeline and all
dependencies. This includes DP3 and wsclean as they are used in the
self-calibration scripts. It also pulls in the LoTSS repositories
`lofar-hba-survey` and `lotss-query` as production use of ddf-pipeline
requires them.

## you will need

* A recent version of singularity
* The recipe in the root directory of this repository

## to build

* Copy the recipe to a build directory -- you do not need any of the rest of the repository.
* Identify a large temporary directory for singularity's working files, say `/my/big/file/system`
* Run `export SINGULARITY_CACHEDIR=/my/big/file/system ; export SINGULARITY_TMPDIR=/my/big/file/system`
* Run `singularity build --fakeroot ddf.sif ddf-py3.singularity`

(You either need to be root or to have 'fakeroot' permission set up. If the former, remove the `--fakeroot` option above.)

## to run

``singularity shell -B /my/data/directory ddf.sif`` will give you an environment inside the singularity, binding in the file systems or directories you plan to work in. Note that to run ddf-pipeline itself you will also need the catalogues used for calibration, masking and bootstrap work. Set the environment variable `DDF_PIPELINE_CATALOGS` to the location of your catalogues.

