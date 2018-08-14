#!/bin/bash
#
#SBATCH -J LB-Impute
#SBATCH --ntasks-per-node=5
#SBATCH --nodes 1
#SBATCH --verbose	      # check man page for more verbosity
#SBATCH --share               # allow for shared resources with other jobs
#SBATCH --time=0              # request that no time limit be imposed
#SBATCH --array=1-36

# must call with dataset name in sbatch
./lb-impute_impute.run $1 $SLURM_ARRAY_TASK_ID
