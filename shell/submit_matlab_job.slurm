#!/bin/bash

# default settings are for tiger

#stop if error
set -e

# Pick a cells_file
matlab_command=$1

#SBATCH -o slurm_outfiles/%A.out
#SBATCH -e slurm_outfiles/%A.err
#SBATCH -t 1439 # one minute shy of a day by default
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=1     # number of tasks per node
#SBATCH --cpus-per-task=7        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --gres=gpu:1             # number of gpus per node


module load matlab/R2017a

matlab -nodesktop -nodisplay -nosplash -r "${matlab_command}"
