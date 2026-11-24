#!/bin/bash
set -eo pipefail

# configuration
export PATH=$HOME/.local/bin:$PATH
export METHOD_CLONE="https://github.com/"
export UV_LINK_MODE=copy
export UV_CACHE_DIR="${SCRATCH:-$HOME/scratch}/uv_cache"
mkdir -p "$UV_CACHE_DIR"

# uv installation
if ! command -v uv > /dev/null 2>&1; then
    curl -LsSf https://astral.sh/uv/install.sh | sh
fi

# modules
set +u
module load python/3.10.4
module load gcc/11.4.1
module load openmpi/4.1.6
module load cmake/3.21.3
module load gsl/2.7.1
set -u

# sources — pinned to the commits validated on Jean Zay
git clone ${METHOD_CLONE}cyriltasse/DDFacet -b HackathonRennes ../DDFacet || true
git -C ../DDFacet fetch --quiet --all || true
git -C ../DDFacet checkout 675ce535422c70bb1499d00a2437657a8bfcb9ff

git clone ${METHOD_CLONE}cyriltasse/killMS -b HackathonRennes ../killMS || true
git -C ../killMS fetch --quiet --all || true
git -C ../killMS checkout c8f178848e6effaede8bdfb97dfe7c3303dedaeb

# python environment
uv venv -p 3.10
source .venv/bin/activate

# packages installation
uv sync --extra mpi-support --refresh-package ddfacet --refresh-package killms --verbose
uv pip install https://github.com/dguibert/LOFARBeam/releases/download/v0.1-10-gc49afaf/lofarbeam-0.0.0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
