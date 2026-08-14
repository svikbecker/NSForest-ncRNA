#!/bin/bash
#SBATCH --job-name=codingRetina
#SBATCH --mem=128g
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --output=../../output_biowulf/terminal/oe_retina_%j.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_coding.py data_retina_all_ds.h5ad retina author_cell_type feature_type
