#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:1
#SBATCH --array=0-6

module purge
module load snakemake/7.7.0

snakemake --cores all --profile profile/snakemake_profile --use-conda -f annotate 