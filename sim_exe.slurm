#!/bin/bash
#SBATCH --job-name=DDF_MPI_%j
#SBATCH --output=log/DDF_MPI_%j.out
#SBATCH --error=log/DDF_MPI_%j.err
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks-per-node=6
#SBATCH --time=00:20:00
#SBATCH --hint=nomultithread
# SLURM options on the Jean Zay supercomputer
#SBATCH --account=doz@cpu
#SBATCH --partition=cpu_p1
# SLURM options on the Adastra supercomputer
##SBATCH --account=cad15976
##SBATCH --constraint=GENOA

# Modules on the Jean Zay supercomputer
module load gsl/2.7.1
module load gcc-native/14
module load craype-network-ofi
module load craype/2.7.35
module load cray-mpich/9.0.1

# Modules on the Adastra supercomputer
#module load gsl/2.7.1
#module load gcc-native/14
#module load craype-network-ofi
#module load craype/2.7.35
#module load cray-mpich/9.0.1

module list

SRC= realpath ../
UV_ENV=$SRC/ddf-pipeline/.venv
SRC_DDF=$SRC
SRC_KMS=$SRC_DDF
SRC_DDFPIPELINE=$SRC_DDF
DATA_DIR=$SCRATCH/ddf-sim-dataset/Sim/
INPUT_DIR=$DATA_DIR
OUTPUT_DIR=$DATA_DIR

source $UV_ENV/bin/activate

export I_MPI_DEBUG=5

export DDF_LOCAL_DEV=1
export DDF_PIPELINE_CATALOGS=$SRC/ddf-pipeline/CATALOGS
export OPENBLAS_NUM_THREADS=1
export OPENBLAS_MAX_THREADS=1
export PYTHONPATH_FIRST=1
export NUMEXPR_MAX_THREADS=40
export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export PYTHONPATH=$SRC_DDF/DDFacet:$PYTHONPATH
export PYTHONPATH=$SRC_KMS/killMS:$PYTHONPATH
export PYTHONPATH=$SRC_KMS/FindCluster:$PYTHONPATH
export PYTHONPATH=$SRC_KMS/DynSpecMS:$PYTHONPATH
export PYTHONPATH=$SRC_KMS/lotss-query:$PYTHONPATH
export PYTHONPATH=$SRC_DDFPIPELINE/ddf-pipeline:$PYTHONPATH
export PYTHONPATH=$SRC_DDFPIPELINE/ddf-pipeline/utils:$PYTHONPATH
export PYTHONPATH=$SRC_DDFPIPELINE/ddf-pipeline/scripts:$PYTHONPATH

export LD_LIBRARY_PATH=$SRC_DDF/DDFacet/DDFacet/cbuild:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SRC_KMS/LOFARBeam/build
# Use on Adastra supercomputer
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAY_MPICH_DIR/lib-abi-mpich

export PATH=$SRC_KMS/killMS/killMS:$PATH
export PATH=$SRC_DDF/DDFacet/DDFacet:$PATH
export PATH=$SRC_KMS/drawMS:$PATH
export PATH=$SRC_KMS/AnalyseDynSpecMS:$PATH
export PATH=$SRC_DDF/DDFacet/SkyModel:$PATH
export PATH=$SRC_KMS/DynSpecMS:$PATH
export PATH=$SRC_KMS/FindCluster:$PATH
export PATH=$SRC_KMS/DIDAP:$PATH
export PATH=$SRC_DDFPIPELINE/ddf-pipeline/scripts:$PATH

export PS1="\[\e[1;92m\][src@${SRC_KMS}]\[\e[m\] \u@\h:\w\$ "

cd $OUTPUT_DIR
ls

# DDF with simulated dataset
srun --verbose $SRC_DDFPIPELINE/ddf-pipeline/scripts/pipeline.py $SRC_DDFPIPELINE/ddf-pipeline/examples/Simul-Hackaton.cfg
