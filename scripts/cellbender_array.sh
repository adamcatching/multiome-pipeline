#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:1
#SBATCH --array=0-246

"""
This takes all of the transferred Cellranger output files and runs cellbender
"""

# Load modules
module load cellbender
module load CUDA/12.1

# Define the directories to run cell bender within
output_file_base=$(echo /data/CARD_singlecell/SN_atlas/cellbender/data/samples/*/)

# Convert the locations of the directories to an array
out_dirs=(`echo ${output_file_base}`);

# Iterate through the array of sample directories
cd ${out_dirs[$SLURM_ARRAY_TASK_ID]}; cellbender remove-background --input raw_feature_bc_matrix.h5 --output cellbender_gex_counts.h5 --fpr 0 --cuda; cd ..