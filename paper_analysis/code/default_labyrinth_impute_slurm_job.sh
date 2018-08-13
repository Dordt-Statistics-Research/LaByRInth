#!/bin/bash
#
#SBATCH -J LB-Impute
#SBATCH --ntasks-per-node=25
#SBATCH --nodes 1
#SBATCH --verbose	      # check man page for more verbosity
#SBATCH --share               # allow for shared resources with other jobs
#SBATCH --time=0              # request that no time limit be imposed

./default_labyrinth_impute.run $1 > $1.out 2>&1
