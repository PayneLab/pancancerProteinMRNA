make_list = [["luad","EGFR"],
["luad","MUC16"],
["luad",'RYR2'],
["luad",'TTN'],
["luad",'LRP1B'],
['luad','KRAS'],
['luad','CSMD3'],
['luad','USH2A'],
['luad','ZFHX4'],
['lscc','TTN'],
['lscc','CSMD3'],
['lscc','USH2A'],
['lscc','ZFHX4'],
['lscc','RYR2'],
['lscc','LRP1B'],
['lscc','SYNE1'],
['lscc','MUC16'],
['lscc','FAM135B'],
['hnscc','TTN'],
['hnscc','CDKN2A'],
['hnscc','FAT1'],
['hnscc','NOTCH1'],
['hnscc','CSMD3'],
['hnscc','DNAH5'],
['hnscc','MUC16'],
['hnscc','RYR2'],
['hnscc','CTNNA2'],
['en','PIK3CA'],
['en','ARID1A'],
['en','PIK3R1'],
['en','KRAS'],
['en','CTNNB1'],
['en','CTCF'],
['en','KMT2B'],
['en','TP53'],
['en','ZFHX3'],
['ccrcc','PBRM1'],
['ccrcc','KDM5C'],
['ccrcc','BAP1'],
['ccrcc','SETD2']]

full_header = """#!/bin/bash

#SBATCH --time=168:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=1536M   # memory per CPU core
#SBATCH --mail-user=kimball.benjamin1@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE"""

for item in make_list:
    file_create_name = "trans_" + item[0] + "_" + item[1] + ".sh"
    command_string = "\n\npython supercomputer_transmutation_effects.py " + item[0] + " " + item[1] + " 10000"
    with open(file_create_name,'w') as writeFile:
        writeFile.write(full_header)
        writeFile.write(command_string)
