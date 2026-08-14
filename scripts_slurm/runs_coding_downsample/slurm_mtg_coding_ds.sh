#!/bin/bash
#SBATCH --job-name=codingMtg_ds
#SBATCH --mem=256g
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --output=../../output_biowulf/terminal/oe_mtg_coding_ds.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_coding_downsample.py data_mtg_ds.h5ad mtg cluster_id biotype
