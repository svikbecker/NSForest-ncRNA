#!/bin/bash
#SBATCH --job-name=codingKidney
#SBATCH --mem=128g
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --output=../../output_biowulf/terminal/oe_kidney_%j.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_coding.py data_kidney_ds_filter.h5ad kidney subclass.full feature_type