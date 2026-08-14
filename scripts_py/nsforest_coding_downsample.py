"""
python nsforest_coding_downsample.py {data file name} {tissue} {cluster header} {subset col}
"""

# libraries
import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import random
import scanpy as sc
import anndata as ad
import nsforest as ns
import io
import contextlib
import warnings

## HEADER FUNCTION
def print_header(title):
    width = 20
    print("\n" + "=" * width)
    print(f"{title.upper():^{width}}")
    print("=" * width)

# dictionary with num to downsample to (num of positive lncRNA for cluster_header annotation),
# as obtained in full-lncRNA run.
ds_dict = {"bm" : 79,
           "breast" : 64,
           "kidney" : 203,
           "liver" : 41,
           "lung" : 83,
           "m1" : 262,
           "mtg" : 386,
           "retina" : 703,
           "spc" : 703}

### CONFIGURATION -- Set the paths to the data folder and output folder
print_header("Configuring Environment")

# change working directory to the script's directory
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
os.chdir(script_dir)

data_folder = "../data_clean/" # path to folder containing the input data file (.h5ad format)
filename = sys.argv[1]
file = data_folder + filename

output_folder = f"../output_biowulf/{sys.argv[2]}/"
preprocessed_folder = "../data_preprocessed/"
store_preprocessed_subset = False # in-script toggle for saving the preprocessed subset 
                                  # (not required for further analysis in this project)

seed = 0 # random seed for reproducibility
np.random.seed(seed) # set np seed
random.seed(seed) # set random seed

cluster_header = sys.argv[3]   # column name in adata.obs that contains the cluster labels
subset_col = sys.argv[4]       # column name in adata.var that contains the gene feature type
subset_gene = "protein_coding" # feature type to subset the data by (EX: "lncRNA" or "protein_coding")

print("Tissue:", sys.argv[2])
print("Cluster Header:", sys.argv[3])

### IMPORT DATA -- Load the data and explore the dataset
print_header("Importing/Loading Data")
adata = sc.read_h5ad(file) # load the data into an AnnData object
print(adata)

### NORMALIZE/SCALE DATA -- Normalize and scale (log1p)
### cannot do this with retina data, already done + no raw counts
if adata.raw is not None and hasattr(adata.raw, "X"):
    print_header("Normalize/Scale Data")

    if sys.argv[2] in ["mtg", "m1"]:
        adata.raw._var.index = adata.raw.var["gene"]
        adata.X = adata.raw[:, adata.var_names].X.copy()
    else:
        adata.X = adata.raw.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    sc.pp.log1p(adata)
    print("Data normalized and log1p transformed.")
else:
    print("No Normalizing/Scaling, data passed has no raw.X")

### CREATE GENE TYPE SUBSET -- Subset the data to only include genes of a specific feature type (EX: "lncRNA" or "protein_coding")
# adata_subset = adata[:, adata.var[subset_col] == subset_gene].copy()
# OR: compute subset_gene mask on the adata object
print_header(f"Create {subset_gene} Subset")
subset_mask = adata.var[subset_col].isin([subset_gene])
adata_subset = adata[:, subset_mask]
print(adata_subset)

# Preprocessing: Cluster Medians
print_header("Preprocess: Cluster Medians")
adata_subset = ns.pp.prep_medians(adata_subset, cluster_header)

# Preprocessing: Binary Scores
print_header("Preprocess: Binary Score")
adata_subset = ns.pp.prep_binary_scores(adata_subset, cluster_header)

## SAVE PREPROCESSED DATA AS NEW .h5ad
print_header("Preprocess: Saving Preprocessed Data")
    
if store_preprocessed_subset:
    filepp = filename.replace(".h5ad", f"_subset_{subset_gene}_preprocessed.h5ad")
    print(f"Saving new anndata object as...\n{preprocessed_folder + filepp}")
    adata_subset.write_h5ad(preprocessed_folder + filepp)
else: 
    print(f"Not saving the preprocessed full protein-coding subset.")

# DIVIDE CODING GENES INTO FOLDS
print_header("Divide Coding Genes into Folds")

# get number of folds
# num_folds = min(int(adata_subset.varm[f'binary_scores_{cluster_header}'].shape[0] / ds_dict[sys.argv[2]]), 15) # limit
# print((f"There are {adata_subset.varm[f'binary_scores_{cluster_header}'].shape[0]} positive coding genes "
#        f"and {ds_dict[sys.argv[2]]} positive lncRNA genes for {cluster_header}. When downsampling to {ds_dict[sys.argv[2]]}, "
#        f"we can create a maximum of {int(adata_subset.varm[f'binary_scores_{cluster_header}'].shape[0] / ds_dict[sys.argv[2]])} folds."))

num_folds = 10

print(f"We will divide the coding genes into {num_folds} folds.")

# build fold dictionary
fold_dict = {}
for i in range(num_folds):
    fold_dict[f"fold_{i+1}"] = random.sample(list(adata_subset.varm[f'binary_scores_{cluster_header}'].index), ds_dict[sys.argv[2]])

# RUN NSFOREST -- Run NSForest on each fold of coding genes
print_header("Run NSForest on Each Fold")

results_list = [] # initialize empty list to store results
columns = ["clusterName", "f_score", "precision", "recall", "onTarget"] # columns to extract
silent_output = io.StringIO() # captures output from nsforest runs

for fold_name, coding_genes in fold_dict.items():
    print(f"Running NSForest on {fold_name}...")
    # subset adata to only include the genes in the current fold
    adata_fold = adata_subset[:, coding_genes].copy()
    
    # run NSForest on the current fold
    with warnings.catch_warnings(action="ignore"):
        with contextlib.redirect_stdout(silent_output):
            results_fold = ns.nsforesting.NSForest(adata_fold, cluster_header, save_supplementary = False, save = False)

    print(f"Completed NS-Forest run on {fold_name}. Extracting results...\n")
    # extract clusterName, f_score, precision, recall, and onTarget from results_fold and add to results_df
    results_fold_df = pd.DataFrame(results_fold)
    results_fold_df["fold"] = fold_name
    results_list.append(results_fold_df)

# CONCATENATE/AGG RESULTS -- Obtain average performance measures across folds
print_header("Concatenating/Aggregating Results")

results_df = pd.concat(results_list)

results_agg = results_df.groupby(['clusterName']).agg(
    avg_precision=pd.NamedAgg(column = "precision", aggfunc="mean"),
    avg_recall=pd.NamedAgg(column = "recall", aggfunc="mean"),
    avg_f_score=pd.NamedAgg(column = "f_score", aggfunc="mean"),
    avg_onTarget=pd.NamedAgg(column = "onTarget", aggfunc="mean"),
).reset_index(inplace=False)

print("First 10 rows of aggregated Results:")
print(results_agg.head(10))

## SAVE FULL AND AGGREGATED RESULTS
print_header("Saving Full and Aggregated Results")
    
# full results
filefull = f"dscoding_results_full_{sys.argv[2]}.csv"
print(f"Saving full results as...\n{output_folder + filefull}")
results_df.to_csv(output_folder + filefull, index=False)

# full results
fileagg = f"dscoding_results_agg_{sys.argv[2]}.csv"
print(f"Saving full results as...\n{output_folder + fileagg}")
results_agg.to_csv(output_folder + fileagg, index=False)


