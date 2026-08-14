#!/bin/bash
#SBATCH --job-name=codingMtg
#SBATCH --mem=64g
#SBATCH --cpus-per-task=8
#SBATCH --output=../../output_biowulf/terminal/oe_mtg_%j.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_coding.py data_mtg_ds.h5ad mtg cluster_id biotype
