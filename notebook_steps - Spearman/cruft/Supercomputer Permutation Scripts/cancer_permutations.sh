#!/bin/bash

#SBATCH --time=120:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=6000M   # memory per CPU core
#SBATCH -J "permutation_10k"   # job name
#SBATCH --mail-user=kimball.benjamin1@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE


python supercomputer_permutation_tests.py ccrcc 10000
python supercomputer_permutation_tests.py en 10000
python supercomputer_permutation_tests.py luad 10000
python supercomputer_permutation_tests.py hnscc 10000
python supercomputer_permutation_tests.py lscc 10000
