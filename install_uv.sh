#!/bin/bash

# install UV
export PATH=$HOME/.local/bin:$PATH
if ! command -v uv; then
    curl -LsSf https://astral.sh/uv/install.sh | sh
fi


# Use on HPC centers (Jean Zay, Adastra, ...)
#export UV_CACHE_DIR=$PWD/uv_cache
# Use on Adastra HPC center
#export CC=gcc

# download DDFacet/killMS and ddf-pipeline
# Use https method on some HPC centers
export METHOD_CLONE="https://github.com/"
export METHOD_CLONE="git@github.com:"


git clone ${METHOD_CLONE}cyriltasse/DDFacet -b MassiveMerge_PR_MergeSSD3_NancepMPI_APP  ../DDFacet || true
git clone ${METHOD_CLONE}cyriltasse/killMS -b APP_Predict_Compress_PolSmooth_HybridSM_OpFit_MultiField_MPI_MultiChain ../killMS || true
# git clone ${METHOD_CLONE}dguibert/ddf-pipeline -b Hackaton_mpipool_test_NancepMPI_Herts ../ddf-pipeline || true

uv venv -p 3.10
source .venv/bin/activate
uv sync --extra mpi-support --refresh-package ddfacet --refresh-package killms --active --verbose


