#!/bin/bash
#SBATCH --job-name=DDF_BIG_%j
#SBATCH --output=DDF_BIG_%j.out
#SBATCH --error=DDF_BIG_%j.err
#SBATCH --account=doz@cpu
#SBATCH --partition=cpu_p1
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --hint=nomultithread
#
# Multi-node launcher (real dataset, and the scaling matrix). One MPI rank per
# node, all cores as threads. Submit from the repo root. Override on the command
# line, e.g. for the scaling study:
#   DDF_SIM_DIR=<sim run dir> DDF_PIPELINE_CFG=<sim cfg> \
#     sbatch --nodes=3 --time=00:40:00 slurm/jean-zay-big.slurm
# For the real run:
#   DDF_SIM_DIR=<big MS dir> DDF_PIPELINE_CFG=$PWD/examples/lotss_mpi-no_low.cfg \
#     sbatch slurm/jean-zay-big.slurm

module load gcc/11.4.1
module load openmpi/4.1.6
module load cmake/3.21.3
module load gsl/2.7.1

module list

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
[ -f "$REPO_DIR/env.sh" ] && REPO_DIR="$(realpath -- "$REPO_DIR/..")"
if [ ! -f "$REPO_DIR/slurm/env.sh" ]; then
    echo "error: cannot find slurm/env.sh under $REPO_DIR; submit from the repo root" >&2
    exit 1
fi
source "$REPO_DIR/slurm/env.sh"

cd "$DATA_DIR"
srun --verbose "$REPO_DIR/scripts/pipeline.py" "$CFG"
rc=$?
if [ "$rc" -ne 0 ]; then
    echo "pipeline exited $rc; skipping smoke check" >&2
    exit "$rc"
fi
if [ "${DDF_SMOKE_CHECK:-0}" = "1" ]; then
    python "$REPO_DIR/scripts/smoke_check.py" "$DATA_DIR"
fi
