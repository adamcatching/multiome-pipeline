#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:1
#SBATCH --array=0-6

module purge
module load apptainer
module load snakemake/7.7.0

# Decision tree for analysis type
# 1. Single cell
# 2. Multiome
# 3. Bulk
# 4. 

# Results for 2
# QC metric values
# Data directory
# Results directory
# Location


# Pull the containers
apptainer pull envs/snapATAC2.sif oras://quay.io/adamcatchingdti/snapatac2
apptainer pull envs/single_cell_cpu.sif oras://quay.io/adamcatchingdti/single_cell_cpu

snakemake --cores all --profile profile/snakemake_profile -n
