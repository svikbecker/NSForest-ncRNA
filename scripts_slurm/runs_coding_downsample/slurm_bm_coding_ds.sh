#!/bin/bash
#SBATCH --job-name=codingBm_ds
#SBATCH --mem=256g
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --output=../../output_biowulf/terminal/oe_bm_coding_ds.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_coding_downsample.py data_bm.h5ad bm Curated_annotation feature_type
