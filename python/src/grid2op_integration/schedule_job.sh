#!/bin/bash

# ======= SLURM HEADER =======
#SBATCH --job-name=grid2op   # Job name
#SBATCH --output=logs/%x_%j.out    # Stdout (%j=jobid)
#SBATCH --error=logs/%x_%j.err     # Stderr
#SBATCH --partition=slurm       # Queue/partition (e.g., gpu, compute)
#SBATCH --account=pymoods       # allocation/account
#SBATCH --time=7-00:00:00        # Time limit hrs:min:sec
## SBATCH --gres=gpu:1          # Number of nodes

module load python/miniconda25.1.1
module load gcc/11.2.0
module load openmpi/4.1.4
module load cmake/3.29.0

source /qfs/projects/gridpack_wind/grid2op_interface/pyenv_gridpack/bin/activate

if [[ "${2}" == "train" ]]; then
    mpirun -np ${SLURM_NTASKS:-1} python -u fidvr_study.py --total-steps 100000 --mode ${2} --log-sim-during-train --max-steps 10 --grid_xml ${1}
else
    mpirun -np ${SLURM_NTASKS:-1} python -u fidvr_study.py \
    --mode ${2} --max-steps 10 --grid_xml "${1}"
fi