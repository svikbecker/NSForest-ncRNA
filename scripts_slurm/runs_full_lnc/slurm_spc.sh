#!/bin/bash
#SBATCH --job-name=lncSpc
#SBATCH --mem=64g
#SBATCH --cpus-per-task=8
#SBATCH --output=../../output_biowulf/terminal/oe_spc_%j.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_full_lnc.py data_spc_ds.h5ad spc ThirdAnnotation feature_type