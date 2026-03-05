#!/bin/bash
#SBATCH --job-name=DDF_MPI
#SBATCH --account= # Setup your account with H100 partition
#SBATCH --constraint=h100
#SBATCH --qos=qos_gpu_h100-t4
#SBATCH --nodes=2
#SBATCH --cpus-per-task=96
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --time=20:00:00
#SBATCH --hint=nomultithread
#SBATCH --output=DDF_MPI_%j.out
#SBATCH --error=DDF_MPI_%j.err

module purge

module load arch/h100
module load anaconda-py3/2024.06
module load gcc/11.4.1
module load openmpi/4.1.6
module load gsl/2.7.1

conda activate ddf_mpi

SRC=$WORK/workspace/dev
SRC_DDF=$SRC/VE_MPI/sources
SRC_KMS=$SRC_DDF
SRC_DDFPIPELINE=$SRC_DDF
DATASET=$WORK/workspace/dataset/LoTSS
CFG=$WORK/workspace/dataset/config/lotss_mpi.cfg
OUTPUT_DIR=$SCRATCH/LoTSS

export DDF_LOCAL_DEV=1
export DDF_PIPELINE_CATALOGS=/lustre/fswork/projects/rech/hfk/uvh54iy/workspace/dev/CATALOGS

export OPENBLAS_NUM_THREADS=1
export OPENBLAS_MAX_THREADS=1
export PYTHONPATH_FIRST=1
export NUMEXPR_MAX_THREADS=96
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
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

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

srun pipeline.py $CFG