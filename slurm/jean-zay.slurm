#!/bin/bash
#SBATCH --job-name=DDF_MPI_%j
#SBATCH --output=DDF_MPI_%j.out
#SBATCH --error=DDF_MPI_%j.err
#SBATCH --account=doz@cpu
#SBATCH --partition=cpu_p1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --time=00:40:00
#SBATCH --hint=nomultithread
#
# Jean Zay cpu_p1 node: 2 x Xeon Gold 6248 (20 cores each) = 40 cores, no SMT.
# Layout: one MPI rank on the node, 40 threads inside it (NCPU_DDF = NCPU_killms
# = 40, taken from the process affinity). The pipeline's MPI layer distributes
# work per node, so extra ranks on one node would replay the same job.

module load gcc/11.4.1
module load openmpi/4.1.6
module load cmake/3.21.3
module load gsl/2.7.1

module list

# SLURM copies the job script to a spool directory, so $BASH_SOURCE does not
# point at the checkout. Resolve from the submit directory instead; submit from
# the repo root:  cd ddf-pipeline && sbatch slurm/jean-zay.slurm
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
