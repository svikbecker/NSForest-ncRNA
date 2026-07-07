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
import celltypist as ct


### CONFIGURATION -- Set the paths to the code folder, data folder, and output folder
code_folder = "/Users/vbecker/NSForest-ncRNA" # path to the NSForest-ncRNA folder
sys.path.insert(0, os.path.abspath(code_folder))

data_folder = "beckersv_data/" # path to folder containing the input data file (.h5ad format)
file = data_folder + "data_spc.h5ad"

output_folder = "beckersv_output/spc/"

to_downsample = False  # True if you want to downsample the dataset to a specific number of cells, 
                       # False otherwise
to_downsample_n = None # number of cells to downsample each cluster to, if to_downsample is True

seed = 0 # random seed for reproducibility
np.random.seed(seed) # set np seed
random.seed(seed) # set random seed

cluster_header = "cell_type" # column name in adata.obs that contains the cluster labels
subset_col = "feature_type"         # column name in adata.var that contains the gene feature type
subset_gene = "lncRNA"              # feature type to subset the data by (EX: "lncRNA" or "protein_coding")


### IMPORT DATA -- Load the data and explore the dataset
adata_raw = sc.read_h5ad(file) # load the data into an AnnData object

### DOWNSAMPLE (OPTIONAL) -- Downsample the dataset to a specific number of cells per cluster
if to_downsample:
    adata = ct.samples.downsample_adata(adata_raw, mode = "each", n_cells = to_downsample_n, by = cluster_header,
                                        random_state = seed, return_index = False)
else:
    adata = adata_raw
    
### CREATE GENE TYPE SUBSET -- Subset the data to only include genes of a specific feature type (EX: "lncRNA" or "protein_coding")# keep only lncRNA genes
# adata_subset = adata[:, adata.var[subset_col] == subset_gene].copy()
# OR: compute lncrna mask on the adata_raw object
lncrna_mask = adata_raw.var[subset_col].isin([subset_gene])
adata_subset = adata[:, lncrna_mask]

### PREPROCESS DATA -- Dendrogram, cluster median, and binary score generation.

## DENDROGRAM
# full adata
if not adata.obsm or "X_pca" not in adata.obsm:
    sc.pp.pca(adata, random_state=seed)

ns.pp.dendrogram(adata, cluster_header, save = "svg", output_folder = output_folder, outputfilename_suffix = f"{cluster_header}_full")
Path(f"dendrogram_{cluster_header}.svg").rename(
    output_folder + f"dendrogram_{cluster_header}_full.svg"
)

# subset adata
if not adata_subset.obsm or "X_pca" not in adata_subset.obsm:
    sc.pp.pca(adata_subset, random_state=seed)

ns.pp.dendrogram(adata_subset, cluster_header, save = "svg", output_folder = output_folder, outputfilename_suffix =  f"{cluster_header}_subset")
Path(f"dendrogram_{cluster_header}.svg").rename(
    output_folder + f"dendrogram_{cluster_header}_subset.svg"
)

## CLUSTER MEDIAN
adata = ns.pp.prep_medians(adata, cluster_header)
adata_subset = ns.pp.prep_medians(adata_subset, cluster_header)

## BINARY SCORE
adata = ns.pp.prep_binary_scores(adata, cluster_header)
adata_subset = ns.pp.prep_binary_scores(adata_subset, cluster_header)

### VISUALIZE PREPROCESSING
# cluster medians (unscaled)
ns.pp.plot_varm(adata, f"medians_{cluster_header}", show = False, nonzero = True, save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_medians_{cluster_header}.png").rename(
    output_folder + f"histogram_medians_{cluster_header}_full.png"
)

ns.pp.plot_varm(adata_subset, f"medians_{cluster_header}", show = False, nonzero = True, save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_medians_{cluster_header}.png").rename(
    output_folder + f"histogram_medians_{cluster_header}_subset.png"
)

# cluster medians (log scale)
ns.pp.plot_varm(adata, f"medians_{cluster_header}", show = False, scale = "log", save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_medians_{cluster_header}.png").rename(
    output_folder + f"histogram_medians_{cluster_header}_full_log.png"
)

ns.pp.plot_varm(adata_subset, f"medians_{cluster_header}", show = False, scale = "log", save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_medians_{cluster_header}.png").rename(
    output_folder + f"histogram_medians_{cluster_header}_subset_log.png"
)

# binary scores (unscaled)
ns.pp.plot_varm(adata, f"binary_scores_{cluster_header}", show = False, nonzero = True, save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_binary_scores_{cluster_header}.png").rename(
    output_folder + f"histogram_binary_scores_{cluster_header}_full.png"
)

ns.pp.plot_varm(adata_subset, f"binary_scores_{cluster_header}", show = False, nonzero = True, save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_binary_scores_{cluster_header}.png").rename(
    output_folder + f"histogram_binary_scores_{cluster_header}_subset.png"
)

# binary scores (log scale)
ns.pp.plot_varm(adata, f"binary_scores_{cluster_header}", show = False, scale = "log", save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_binary_scores_{cluster_header}.png").rename(
    output_folder + f"histogram_binary_scores_{cluster_header}_full_log.png"
)

ns.pp.plot_varm(adata_subset, f"binary_scores_{cluster_header}", show = False, scale = "log", save = "svg", output_folder = output_folder)
Path(output_folder + f"histogram_binary_scores_{cluster_header}.png").rename(
    output_folder + f"histogram_binary_scores_{cluster_header}_subset_log.png"
)

### RUN NSFOREST -- Run NSForest on the full dataset and the subset dataset# full adata
outputfilename_prefix = cluster_header
results = ns.nsforesting.NSForest(adata, cluster_header, save_supplementary = True, save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

# subset adata
outputfilename_prefix_subset = cluster_header + "_subset_" + subset_gene
results_subset = ns.nsforesting.NSForest(adata_subset, cluster_header, save_supplementary = True, save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

### VISUALIZE NSFOREST RESULTS -- Visualize the results of NSForest for the full dataset and the subset dataset
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

# dotplot
ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
ns.pl.dotplot(adata_subset, markers_dict_subset, cluster_header, dendrogram = dendrogram_subset, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix_subset)

# stacked violin plot
ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
ns.pl.stackedviolin(adata_subset, markers_dict_subset, cluster_header, dendrogram = dendrogram_subset, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix_subset)

# heatmap
ns.pl.matrixplot(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
ns.pl.matrixplot(adata_subset, markers_dict_subset, cluster_header, dendrogram = dendrogram_subset, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix_subset)

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