#!/bin/bash

#SBATCH --time=168:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=1536M   # memory per CPU core
#SBATCH --mail-user=kimball.benjamin1@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

python supercomputer_transmutation_effects.py luad ZFHX4 10000