#!/bin/bash

SCRIPT_DIR="$(dirname "$(realpath -- "${BASH_SOURCE[0]}")")"

export DDF_PIPELINE_CATALOGS=$USER/CATALOGS
export SOURCE_DIR="$(realpath -- "$SCRIPT_DIR/..")"

source $SCRIPT_DIR/.venv/bin/activate

export OPENBLAS_NUM_THREADS=1
export OPENBLAS_MAX_THREADS=1
export PYTHONPATH_FIRST=1
export NUMEXPR_MAX_THREADS=96
export PYTHONHASHSEED=0

export PYTHONPATH=$SOURCE_DIR/DDFacet:$PYTHONPATH
export PYTHONPATH=$SOURCE_DIR/killMS:$PYTHONPATH
export PYTHONPATH=$SOURCE_DIR/ddf-pipeline:$PYTHONPATH
export PYTHONPATH=$SOURCE_DIR/ddf-pipeline/utils:$PYTHONPATH
export PYTHONPATH=$SOURCE_DIR/ddf-pipeline/scripts:$PYTHONPATH
export PYTHONPATH=$SOURCE_DIR/lotss-query:$PYTHONPATH


export LD_LIBRARY_PATH=$SCRIPT_DIR/.venv/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$SOURCE_DIR/DDFacet/DDFacet/cbuild:$LD_LIBRARY_PATH


export PATH=$SOURCE_DIR/killMS/killMS:$PATH
export PATH=$SOURCE_DIR/DDFacet/DDFacet:$PATH
export PATH=$SOURCE_DIR/DDFacet/SkyModel:$PATH
export PATH=$SOURCE_DIR/DynSpecMS:$PATH
export PATH=$SOURCE_DIR/FindCluster:$PATH
export PATH=$SOURCE_DIR/ddf-pipeline/scripts:$PATH

export PS1="\[\e[1;92m\][src@${SOURCE_DIR}]\[\e[m\] \u@\h:\w\$ "

echo "DDF exec: "$(which DDF.py)
echo "kMS exec: "$(which kMS.py)
