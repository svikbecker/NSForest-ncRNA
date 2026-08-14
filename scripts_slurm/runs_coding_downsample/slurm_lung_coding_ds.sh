#!/bin/bash
#SBATCH --job-name=codingLung_ds
#SBATCH --mem=256g
#SBATCH --cpus-per-task=16
#SBATCH --time=3-00:00:00
#SBATCH --output=../../output_biowulf/terminal/oe_lung_coding_ds.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_coding_downsample.py data_lung_ds.h5ad lung ann_finest_level feature_type
