#!/bin/bash

# Use on Jean Zay supercomputer
#module load gcc/11.4.1
#module load openmpi/4.1.6
#module load cmake/3.21.3
#module load gsl/2.7.1

# Use on Adastra supercomputer
#module load gsl/2.7.1
#module load gcc-native/14
#module load craype-network-ofi
#module load craype/2.7.35
#module load cray-mpich/9.0.1

# install UV
export PATH=$HOME/.local/bin:$PATH
if ! command -v uv; then
    curl -LsSf https://astral.sh/uv/install.sh | sh
fi

# Use on supercomputers (Jean Zay, Adastra, ...)
#export UV_CACHE_DIR=$PWD/uv_cache
#export CC=gcc
#export CXX=g++

# download DDFacet/killMS and ddf-pipeline
# Use https method on supercomputers
export METHOD_CLONE="https://github.com/"
export METHOD_CLONE="git@github.com:"

git clone ${METHOD_CLONE}cyriltasse/DDFacet -b HackathonRennes  ../DDFacet || true
git clone ${METHOD_CLONE}cyriltasse/killMS -b HackathonRennes ../killMS || true
# git clone ${METHOD_CLONE}dguibert/ddf-pipeline -b Hackaton_mpipool_test_NancepMPI_Herts ../ddf-pipeline || true

uv venv -p 3.10
source .venv/bin/activate
uv sync --extra mpi-support --refresh-package ddfacet --refresh-package killms --active --verbose
uv pip install https://github.com/dguibert/LOFARBeam/releases/download/v0.1-10-gc49afaf/lofarbeam-0.0.0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
