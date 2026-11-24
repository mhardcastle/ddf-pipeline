#!/bin/bash
#SBATCH --job-name=DDF_MPI_%j
#SBATCH --output=DDF_MPI_%j.out
#SBATCH --error=DDF_MPI_%j.err
#SBATCH --account=cad15976
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --threads-per-core=1
#SBATCH --exclusive
#SBATCH --time=00:20:00
#
# Adastra GENOA scalar node: 2 x EPYC 9654 (96 cores each) = 192 cores.
# Layout: one MPI rank on the node, 192 threads inside it (NCPU_DDF =
# NCPU_killms = 192, taken from the process affinity). The pipeline's MPI layer
# distributes work per node, so extra ranks on one node would replay the same job.
# Cheaper shared-mode run: drop --exclusive and add --cpus-per-task=48 on the
# sbatch command line; the job then lands in the genoa-shared partition.

module load craype/2.7.35
module load craype-network-ofi
module load craype-x86-genoa
module load gcc-native/14
module load cray-mpich/9.0.1
module load gsl/2.7.1

module list

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "$SCRIPT_DIR/env.sh"

# The venv's mpi4py is built against an ABI-compatible MPICH: expose Cray MPICH's
# ABI library to the loader.
if [ -n "${CRAY_MPICH_DIR:-}" ]; then
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAY_MPICH_DIR/lib-abi-mpich
fi

cd "$DATA_DIR"
srun "$REPO_DIR/scripts/pipeline.py" "$CFG"
