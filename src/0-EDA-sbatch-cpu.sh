#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=300g
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:50

# IMPORTS
module purge
source myconda

# RUN
python 0-EDA.py