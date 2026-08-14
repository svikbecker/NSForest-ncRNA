#!/bin/bash
#SBATCH --job-name=codingBreast_ds
#SBATCH --mem=256g
#SBATCH --cpus-per-task=16
#SBATCH --time=3-00:00:00
#SBATCH --output=../../output_biowulf/terminal/oe_breast_coding_ds.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_coding_downsample.py data_breast_ds_filter.h5ad breast level2 feature_type
