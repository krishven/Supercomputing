#!/bin/bash
#SBATCH -J mpipm  # Job Name
#SBATCH -o Test.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 5	        # Total number of  tasks requested
#SBATCH -N 5
#SBATCH -p development  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 01:28:00     # Run time (hh:mm:ss) - 1.5 hours
export CILK_NWORKERS=32
ibrun tacc_affinity ./dsort > distsort_op
