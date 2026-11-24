#!/bin/bash
# Shared environment for the ddf-pipeline SLURM launchers (slurm/jean-zay.slurm,
# slurm/adastra.slurm). Source it from a launcher; everything platform-specific
# (SBATCH options, module loads) stays in the launcher.

# Resolve the checkout from this file's own location, so the job does not depend
# on the directory sbatch was called from.
ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_DIR="$(realpath -- "$ENV_DIR/..")"   # the ddf-pipeline checkout
SRC="$(realpath -- "$ENV_DIR/../..")"     # holds the DDFacet, killMS, ... siblings

UV_ENV="$REPO_DIR/.venv"

# Simulated dataset. Override from the command line, e.g.
#   DDF_SIM_DIR=/path/to/Sim sbatch slurm/jean-zay.slurm
# SCRATCH is the Jean Zay name, SCRATCHDIR the Adastra one.
DATA_DIR="${DDF_SIM_DIR:-${SCRATCH:-${SCRATCHDIR:-}}/ddf-sim-dataset/Sim}"
INPUT_DIR="$DATA_DIR"
OUTPUT_DIR="$DATA_DIR"
CFG="${DDF_PIPELINE_CFG:-$REPO_DIR/examples/Simul-Hackaton.cfg}"

source "$UV_ENV/bin/activate"

# Intra-node parallelism comes from --NCPU (DDFacet/kMS workers), not from BLAS/OpenMP.
export OPENBLAS_NUM_THREADS=1
export OPENBLAS_MAX_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMEXPR_MAX_THREADS="${SLURM_CPUS_PER_TASK:-40}"
export PYTHONPATH_FIRST=1
export PYTHONHASHSEED=0

export DDF_LOCAL_DEV=1
export DDF_PIPELINE_CATALOGS="${DDF_PIPELINE_CATALOGS:-$REPO_DIR/CATALOGS}"

export PYTHONPATH=$SRC/DDFacet:$PYTHONPATH
export PYTHONPATH=$SRC/killMS:$PYTHONPATH
export PYTHONPATH=$SRC/killMS/FindCluster:$PYTHONPATH
export PYTHONPATH=$SRC/killMS/DynSpecMS:$PYTHONPATH
export PYTHONPATH=$SRC/lotss-query:$PYTHONPATH
export PYTHONPATH=$REPO_DIR:$PYTHONPATH
export PYTHONPATH=$REPO_DIR/utils:$PYTHONPATH
export PYTHONPATH=$REPO_DIR/scripts:$PYTHONPATH

export LD_LIBRARY_PATH=$SRC/DDFacet/DDFacet/cbuild:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SRC/killMS/LOFARBeam/build

export PATH=$SRC/killMS/killMS:$PATH
export PATH=$SRC/DDFacet/DDFacet:$PATH
export PATH=$SRC/killMS/drawMS:$PATH
export PATH=$SRC/killMS/AnalyseDynSpecMS:$PATH
export PATH=$SRC/DDFacet/SkyModel:$PATH
export PATH=$SRC/killMS/DynSpecMS:$PATH
export PATH=$SRC/killMS/FindCluster:$PATH
export PATH=$SRC/killMS/DIDAP:$PATH
export PATH=$REPO_DIR/scripts:$PATH

# Printed to the job's .out file: first thing to check when a run cannot find its
# inputs.
echo "env.sh: REPO_DIR=$REPO_DIR"
echo "env.sh: SRC=$SRC"
echo "env.sh: DATA_DIR=$DATA_DIR"
echo "env.sh: CFG=$CFG"
echo "env.sh: SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-unset}"
