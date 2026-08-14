#!/bin/bash
#SBATCH --job-name=codingKidney_ds
#SBATCH --mem=256g
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --output=../../output_biowulf/terminal/oe_kidney_coding_ds.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_coding_downsample.py data_kidney_ds_filter.h5ad kidney subclass.full feature_type