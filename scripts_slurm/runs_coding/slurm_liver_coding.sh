#!/bin/bash
#SBATCH --job-name=codingLiver
#SBATCH --mem=128g
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --output=../../output_biowulf/terminal/oe_liver_%j.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_coding.py data_liver.h5ad liver Curated_annotation feature_type
