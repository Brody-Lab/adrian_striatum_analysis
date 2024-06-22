#!/bin/bash
#SBATCH --array=0
#SBATCH -o slurm_outfiles/%A_%a.out
#SBATCH -e slurm_outfiles/%A_%a.err;
#SBATCH -t 180                    # 3 hours by default
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=1      # number of tasks per node
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=7G         # important for long sessions or small time binning
#SBATCH --mail-type=FAIL,TIME_LIMIT,ARRAY_TASKS,REQUEUE         # email 
#SBATCH --mail-user=abondy@princeton.edu

pwd; hostname; date
env | sort

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command failed with exit code $?."' EXIT

matlab_command=$1

if [ $SLURM_ARRAY_TASK_COUNT -eq 1 ]; then
  echo "Matlab command is: ${matlab_command}"
  matlab -nodisplay -batch "${matlab_command}" || exit 202
else
  echo "Matlab command is: id=$SLURM_ARRAY_TASK_ID;${matlab_command};exit"
  matlab -nodisplay -batch "id=$SLURM_ARRAY_TASK_ID;${matlab_command};exit" || exit 202
fi

echo "Finished executing matlab command successfully!"

date