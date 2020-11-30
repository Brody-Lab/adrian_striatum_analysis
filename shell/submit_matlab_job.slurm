#!/bin/bash
#SBATCH --array=0
#SBATCH -o slurm_outfiles/%A_%a.out
#SBATCH -e slurm_outfiles/%A_%a.err;
#SBATCH -t 300                   # 5 hrs by default
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=1      # number of tasks per node
#SBATCH --cpus-per-task=7        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --partition=gpu          # send to gpu partition
#SBATCH --mem-per-cpu=16G         # important for long sessions or small time binning
#SBATCH --mail-type=FAIL,TIME_LIMIT,ARRAY_TASKS         # email only upon error or time limit
#SBATCH --mail-user=abondy@princeton.edu

pwd; hostname; date
env | grep SLURM | sort

# these default settings are for tiger

#abort immediately upon error
set -e

matlab_command=$1

module purge
module load matlab/R2017a

if [ $SLURM_ARRAY_TASK_COUNT -eq 1 ]; then
  echo "Matlab command is: ${matlab_command}"
  matlab -nodesktop -nodisplay -nosplash -r "${matlab_command}"
else
  echo "Matlab command is: id=$SLURM_ARRAY_TASK_ID;${matlab_command}"
  matlab -nodesktop -nodisplay -nosplash -r "id=$SLURM_ARRAY_TASK_ID;${matlab_command}"
fi

date
