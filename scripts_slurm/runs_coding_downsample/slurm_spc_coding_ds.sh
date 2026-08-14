#!/bin/bash
#SBATCH --job-name=codingSpc_ds
#SBATCH --mem=128g
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --output=../../output_biowulf/terminal/oe_spc_coding_ds.out

module purge
source myconda
conda activate nsforest

python ../../scripts_py/nsforest_coding_downsample.py data_spc_ds.h5ad spc ThirdAnnotation feature_type