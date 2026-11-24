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
#SBATCH --time=00:20:00
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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "$SCRIPT_DIR/env.sh"

cd "$DATA_DIR"
srun --verbose "$REPO_DIR/scripts/pipeline.py" "$CFG"
