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

python supercomputer_transmutation_effects.py luad TP53 10000
python supercomputer_transmutation_effects.py luad EGFR 10000
python supercomputer_transmutation_effects.py luad MUC16 10000
python supercomputer_transmutation_effects.py luad RYR2 10000
python supercomputer_transmutation_effects.py luad TTN 10000
python supercomputer_transmutation_effects.py luad LRP1B 10000
python supercomputer_transmutation_effects.py luad KRAS 10000
python supercomputer_transmutation_effects.py luad CSMD3 10000
python supercomputer_transmutation_effects.py luad USH2A 10000
python supercomputer_transmutation_effects.py luad ZFHX4 10000
python supercomputer_transmutation_effects.py lscc TP53 10000
python supercomputer_transmutation_effects.py lscc TTN 10000
python supercomputer_transmutation_effects.py lscc CSMD3 10000
python supercomputer_transmutation_effects.py lscc USH2A 10000
python supercomputer_transmutation_effects.py lscc ZFHX4 10000
python supercomputer_transmutation_effects.py lscc RYR2 10000
python supercomputer_transmutation_effects.py lscc LRP1B 10000
python supercomputer_transmutation_effects.py lscc SYNE1 10000
python supercomputer_transmutation_effects.py lscc MUC16 10000
python supercomputer_transmutation_effects.py lscc FAM135B 10000
python supercomputer_transmutation_effects.py hnscc TP53 10000
python supercomputer_transmutation_effects.py hnscc TTN 10000
python supercomputer_transmutation_effects.py hnscc CDKN2A 10000
python supercomputer_transmutation_effects.py hnscc FAT1 10000
python supercomputer_transmutation_effects.py hnscc NOTCH1 10000
python supercomputer_transmutation_effects.py hnscc CSMD3 10000
python supercomputer_transmutation_effects.py hnscc DNAH5 10000
python supercomputer_transmutation_effects.py hnscc MUC16 10000
python supercomputer_transmutation_effects.py hnscc RYR2 10000
python supercomputer_transmutation_effects.py hnscc CTNNA2 10000
python supercomputer_transmutation_effects.py en PTEN 10000
python supercomputer_transmutation_effects.py en PIK3CA 10000
python supercomputer_transmutation_effects.py en ARID1A 10000
python supercomputer_transmutation_effects.py en PIK3R1 10000
python supercomputer_transmutation_effects.py en KRAS 10000
python supercomputer_transmutation_effects.py en CTNNB1 10000
python supercomputer_transmutation_effects.py en CTCF 10000
python supercomputer_transmutation_effects.py en KMT2B 10000
python supercomputer_transmutation_effects.py en TP53 10000
python supercomputer_transmutation_effects.py en ZFHX3 10000
python supercomputer_transmutation_effects.py ccrcc VHL 10000
python supercomputer_transmutation_effects.py ccrcc PBRM1 10000
python supercomputer_transmutation_effects.py ccrcc KDM5C 10000
python supercomputer_transmutation_effects.py ccrcc BAP1 10000
python supercomputer_transmutation_effects.py ccrcc SETD2 10000
python supercomputer_transmutation_effects.py ccrcc TTN 10000
python supercomputer_transmutation_effects.py ccrcc CSMD3 10000
python supercomputer_transmutation_effects.py ccrcc USH2A 10000
python supercomputer_transmutation_effects.py ccrcc MUC16 10000
python supercomputer_transmutation_effects.py ccrcc BIRC6 10000
