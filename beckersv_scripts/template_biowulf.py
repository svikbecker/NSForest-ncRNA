# libraries
import sys
import os
from pathlib import Path
import numpy as np
import random
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import plotly.io as pio
pio.renderers.default = "notebook"
import nsforest as ns
from nsforest import utils

## HEADER FUNCTION!
def print_header(title):
    width = 20
    print("\n" + "=" * width)
    print(f"{title.upper():^{width}}")
    print("=" * width)

### CONFIGURATION -- Set the paths to the code folder, data folder, and output folder
print_header("Configuring Environment")

code_folder = "[]/" # path to the NSForest-ncRNA folder
sys.path.insert(0, os.path.abspath(code_folder))

data_folder = "[]/" # path to folder containing the input data file (.h5ad format)
filename = "[].h5ad"
file = data_folder + filename

output_folder = "[]/"
preprocessed_folder = "[]/"

seed = 0 # random seed for reproducibility
np.random.seed(seed) # set np seed
random.seed(seed) # set random seed

cluster_header = "author_cell_type" # column name in adata.obs that contains the cluster labels
subset_col = "feature_type"         # column name in adata.var that contains the gene feature type
subset_gene = "lncRNA"              # feature type to subset the data by (EX: "lncRNA" or "protein_coding")

### IMPORT DATA -- Load the data and explore the dataset
print_header("Importing/Loading Data")
adata = sc.read_h5ad(file) # load the data into an AnnData object

### NORMALIZE/SCALE DATA -- Normalize and scale (log1p)
print_header("Normalize/Scale Data")

adata.layers["counts"] = adata.raw.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4, layer='counts', inplace=True)
sc.pp.log1p(adata, layer='counts')
adata.X = adata.layers["counts"]
    
### CREATE GENE TYPE SUBSET -- Subset the data to only include genes of a specific feature type (EX: "lncRNA" or "protein_coding")# keep only lncRNA genes
# adata_subset = adata[:, adata.var[subset_col] == subset_gene].copy()
# OR: compute lncrna mask on the adata object
print_header("Create lncRNA Subset")
lncrna_mask = adata.var[subset_col].isin([subset_gene])
adata_subset = adata[:, lncrna_mask]

### PREPROCESS DATA -- Dendrogram, cluster median, and binary score generation.

## DENDROGRAM
print_header("Preprocess: Dendrograms")
# full adata
if not adata.obsm or "X_pca" not in adata.obsm:
    sc.pp.pca(adata, random_state=seed)

ns.pp.dendrogram(adata, cluster_header, save = "svg", pl_kwargs = {'show': False}, output_folder = output_folder, outputfilename_suffix = f"{cluster_header}_full")
Path(code_folder + f"dendrogram_{cluster_header}_full.svg").rename(
    output_folder + f"dendrogram_{cluster_header}_full.svg"
)

# subset adata
if not adata_subset.obsm or "X_pca" not in adata_subset.obsm:
    sc.pp.pca(adata_subset, random_state=seed)

ns.pp.dendrogram(adata_subset, cluster_header, pl_kwargs = {'show': False}, save = "svg", output_folder = output_folder, outputfilename_suffix =  f"{cluster_header}_subset")
Path(code_folder + f"dendrogram_{cluster_header}_subset.svg").rename(
    output_folder + f"dendrogram_{cluster_header}_subset.svg"
)

## CLUSTER MEDIAN
print_header("Preprocess: Cluster Medians")
adata = ns.pp.prep_medians(adata, cluster_header)
adata_subset = ns.pp.prep_medians(adata_subset, cluster_header)

## BINARY SCORE
print_header("Preprocess: Binary Score")
adata = ns.pp.prep_binary_scores(adata, cluster_header)
adata_subset = ns.pp.prep_binary_scores(adata_subset, cluster_header)

## SAVE PREPROCESSED DATA AS NEW .h5ad
print_header("Preprocess: Saving Preprocessed Data")
# full adata
filepp = filename.replace(".h5ad", "_preprocessed.h5ad")
print(f"Saving new anndata object as...\n{preprocessed_folder + filepp}")
adata.write_h5ad(preprocessed_folder + filepp)
    
# subset adata
filepp = file.replace(".h5ad", f"_subset_{subset_gene}_preprocessed.h5ad")
print(f"Saving new anndata object as...\n{filepp}")
adata_subset.write_h5ad(preprocessed_folder + filepp)

### VISUALIZE PREPROCESSING
print_header("Preprocess: Visualizing Preprocessed Data")
# cluster medians (unscaled)
ns.pp.plot_varm(adata, f"medians_{cluster_header}", show = False, nonzero = True, save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_medians_{cluster_header}.svg").rename(
    output_folder + f"histogram_medians_{cluster_header}_full.svg"
)

ns.pp.plot_varm(adata_subset, f"medians_{cluster_header}", show = False, nonzero = True, save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_medians_{cluster_header}.svg").rename(
    output_folder + f"histogram_medians_{cluster_header}_subset.svg"
)

# cluster medians (log scale)
ns.pp.plot_varm(adata, f"medians_{cluster_header}", show = False, scale = "log", save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_medians_{cluster_header}.svg").rename(
    output_folder + f"histogram_medians_{cluster_header}_full_log.svg"
)

ns.pp.plot_varm(adata_subset, f"medians_{cluster_header}", show = False, scale = "log", save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_medians_{cluster_header}.svg").rename(
    output_folder + f"histogram_medians_{cluster_header}_subset_log.svg"
)

# binary scores (unscaled)
ns.pp.plot_varm(adata, f"binary_scores_{cluster_header}", show = False, nonzero = True, save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_binary_scores_{cluster_header}.svg").rename(
    output_folder + f"histogram_binary_scores_{cluster_header}_full.svg"
)

ns.pp.plot_varm(adata_subset, f"binary_scores_{cluster_header}", show = False, nonzero = True, save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_binary_scores_{cluster_header}.svg").rename(
    output_folder + f"histogram_binary_scores_{cluster_header}_subset.svg"
)

# binary scores (log scale)
ns.pp.plot_varm(adata, f"binary_scores_{cluster_header}", show = False, scale = "log", save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_binary_scores_{cluster_header}.svg").rename(
    output_folder + f"histogram_binary_scores_{cluster_header}_full_log.svg"
)

ns.pp.plot_varm(adata_subset, f"binary_scores_{cluster_header}", show = False, scale = "log", save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_binary_scores_{cluster_header}.svg").rename(
    output_folder + f"histogram_binary_scores_{cluster_header}_subset_log.svg"
)

### RUN NSFOREST -- Run NSForest on the full dataset and the subset dataset
# full adata
print_header("NS-Forest: Full Data Run")
outputfilename_prefix = cluster_header
results = ns.nsforesting.NSForest(adata, cluster_header, save_supplementary = True, save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

# subset adata
print_header("NS-Forest: lncRNA Subset Run")
outputfilename_prefix_subset = cluster_header + "_subset_" + subset_gene
results_subset = ns.nsforesting.NSForest(adata_subset, cluster_header, save_supplementary = True, save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

### VISUALIZE NSFOREST RESULTS -- Visualize the results of NSForest for the full dataset and the subset dataset
print_header("Preparing for Visualization")
to_plot = results.copy()
to_plot_subset = results_subset.copy()

# dendrogram full adata
dendrogram = [] # custom dendrogram order
dendrogram = list(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])
to_plot["clusterName"] = to_plot["clusterName"].astype("category")
to_plot["clusterName"] = to_plot["clusterName"].cat.set_categories(dendrogram)
to_plot = to_plot.sort_values("clusterName")
to_plot = to_plot.rename(columns = {"NSForest_markers": "markers"})

# dendrogram subset adata
dendrogram_subset = [] # custom dendrogram order
dendrogram_subset = list(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])
to_plot_subset["clusterName"] = to_plot_subset["clusterName"].astype("category")
to_plot_subset["clusterName"] = to_plot_subset["clusterName"].cat.set_categories(dendrogram_subset)
to_plot_subset = to_plot_subset.sort_values("clusterName")
to_plot_subset = to_plot_subset.rename(columns = {"NSForest_markers": "markers"})

# marker dictionary full adata
markers_dict = dict(zip(to_plot["clusterName"], to_plot["markers"]))
markers_dict

# marker dictionary subset adata
markers_dict_subset = dict(zip(to_plot_subset["clusterName"], to_plot_subset["markers"]))
markers_dict_subset

print_header("Visualization: Selected Marker")

# dotplot
ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = "svg", pl_kwargs = {'show': False}, output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
ns.pl.dotplot(adata_subset, markers_dict_subset, cluster_header, dendrogram = dendrogram_subset, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix_subset)

# stacked violin plot
ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = "svg", pl_kwargs = {'show': False}, output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
ns.pl.stackedviolin(adata_subset, markers_dict_subset, cluster_header, dendrogram = dendrogram_subset, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix_subset)

# heatmap
ns.pl.matrixplot(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = "svg", pl_kwargs = {'show': False}, output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
ns.pl.matrixplot(adata_subset, markers_dict_subset, cluster_header, dendrogram = dendrogram_subset, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix_subset)

print_header("Visualization: Classification Metrics")

# classification metrics
ns.pl.boxplot(results, ["f_score", "precision", "recall", "onTarget"], save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
ns.pl.boxplot(results_subset, ["f_score", "precision", "recall", "onTarget"], save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

# metrics vs cluster size
# f_score
ns.pl.scatter_w_clusterSize(results, "f_score", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
ns.pl.scatter_w_clusterSize(results_subset, "f_score", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

# precision
ns.pl.scatter_w_clusterSize(results, "precision", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
ns.pl.scatter_w_clusterSize(results_subset, "precision", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

# recall
ns.pl.scatter_w_clusterSize(results, "recall", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
ns.pl.scatter_w_clusterSize(results_subset, "recall", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

# onTarget
ns.pl.scatter_w_clusterSize(results, "onTarget", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)
ns.pl.scatter_w_clusterSize(results_subset, "onTarget", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)