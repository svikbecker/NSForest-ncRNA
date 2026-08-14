"""
python nsforest_coding.py {data file name} {structure} {cluster header} {subset col}
"""

# libraries
import sys
import os
from pathlib import Path
import numpy as np
import random
import scanpy as sc
import anndata as ad
import plotly.io as pio
pio.renderers.default = "notebook"
import nsforest as ns

## HEADER FUNCTION!
def print_header(title):
    width = 20
    print("\n" + "=" * width)
    print(f"{title.upper():^{width}}")
    print("=" * width)

### CONFIGURATION -- Set the paths to the data folder and output folder
print_header("Configuring Environment")

# change working directory to the script's directory
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
os.chdir(script_dir)

data_folder = "../data_clean/" # path to folder containing the input data file (.h5ad format)
filename = sys.argv[1]
file = data_folder + filename

output_folder = f"../output/{sys.argv[2]}/"
preprocessed_folder = "../data_preprocessed/"
store_preprocessed_subset = False # in-script toggle for saving the preprocessed subset 
                                  # (not required for further analysis in this project)

seed = 0 # random seed for reproducibility
np.random.seed(seed) # set np seed
random.seed(seed) # set random seed

cluster_header = sys.argv[3]   # column name in adata.obs that contains the cluster labels
subset_col = sys.argv[4]       # column name in adata.var that contains the gene feature type
subset_gene = "protein_coding" # feature type to subset the data by (EX: "lncRNA" or "protein_coding")

### IMPORT DATA -- Load the data and explore the dataset
print_header("Importing/Loading Data")
adata = sc.read_h5ad(file) # load the data into an AnnData object

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
else:
    print_header("No Normalizing/Scaling, data passed has no raw.X")
    
### CREATE GENE TYPE SUBSET -- Subset the data to only include genes of a specific feature type (EX: "lncRNA" or "protein_coding")
# adata_subset = adata[:, adata.var[subset_col] == subset_gene].copy()
# OR: compute subset_gene mask on the adata object
print_header(f"Create {subset_gene} Subset")
subset_mask = adata.var[subset_col].isin([subset_gene])
adata_subset = adata[:, subset_mask]

### PREPROCESS DATA -- Dendrogram, cluster median, and binary score generation.

## DENDROGRAM
print_header("Preprocess: Dendrograms")
# subset adata
if not adata_subset.obsm or "X_pca" not in adata_subset.obsm:
    sc.pp.pca(adata_subset, random_state=seed)

ns.pp.dendrogram(adata_subset, cluster_header, pl_kwargs = {'show': False}, save = "svg", output_folder = output_folder, outputfilename_suffix =  f"{cluster_header}_subset")
Path(f"dendrogram_{cluster_header}_subset.svg").rename(
    output_folder + f"dendrogram_{cluster_header}_subset.svg"
)

## CLUSTER MEDIAN
print_header("Preprocess: Cluster Medians")
adata_subset = ns.pp.prep_medians(adata_subset, cluster_header)

## BINARY SCORE
print_header("Preprocess: Binary Score")
adata_subset = ns.pp.prep_binary_scores(adata_subset, cluster_header)

## SAVE PREPROCESSED DATA AS NEW .h5ad
print_header("Preprocess: Saving Preprocessed Data")
    
# subset adata
if store_preprocessed_subset:
    filepp = filename.replace(".h5ad", f"_subset_{subset_gene}_preprocessed.h5ad")
    print(f"Saving new anndata object as...\n{preprocessed_folder + filepp}")
    adata_subset.write_h5ad(preprocessed_folder + filepp)
else: 
    print(f"Not saving the preprocessed full protein-coding subset.")

### RUN NSFOREST -- Run NSForest on the full dataset and the subset dataset

# subset adata
print_header(f"NS-Forest: {subset_gene} Subset Run")
outputfilename_prefix_subset = cluster_header + "_subset_" + subset_gene
results_subset = ns.nsforesting.NSForest(adata_subset, cluster_header, save_supplementary = True, save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

### VISUALIZE NSFOREST RESULTS -- Visualize the results of NSForest for the full dataset and the subset dataset
print_header("Preparing for Visualization")
to_plot_subset = results_subset.copy()

# dendrogram subset adata
dendrogram_subset = [] # custom dendrogram order
dendrogram_subset = list(adata_subset.uns["dendrogram_" + cluster_header]["categories_ordered"])
to_plot_subset["clusterName"] = to_plot_subset["clusterName"].astype("category")
to_plot_subset["clusterName"] = to_plot_subset["clusterName"].cat.set_categories(dendrogram_subset)
to_plot_subset = to_plot_subset.sort_values("clusterName")
to_plot_subset = to_plot_subset.rename(columns = {"NSForest_markers": "markers"})

# marker dictionary subset adata
markers_dict_subset = dict(zip(to_plot_subset["clusterName"], to_plot_subset["markers"]))
markers_dict_subset

print_header("Visualization: Selected Marker")

# dotplot
ns.pl.dotplot(adata_subset, markers_dict_subset, cluster_header, dendrogram = dendrogram_subset, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix_subset)

# stacked violin plot
ns.pl.stackedviolin(adata_subset, markers_dict_subset, cluster_header, dendrogram = dendrogram_subset, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix_subset)

# heatmap
ns.pl.matrixplot(adata_subset, markers_dict_subset, cluster_header, dendrogram = dendrogram_subset, save = "svg", output_folder = output_folder, outputfilename_suffix = outputfilename_prefix_subset)

print_header("Visualization: Classification Metrics")

# classification metrics
ns.pl.boxplot(results_subset, ["f_score", "precision", "recall", "onTarget"], save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

# metrics vs cluster size
# f_score
ns.pl.scatter_w_clusterSize(results_subset, "f_score", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

# precision
ns.pl.scatter_w_clusterSize(results_subset, "precision", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

# recall
ns.pl.scatter_w_clusterSize(results_subset, "recall", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)

# onTarget
ns.pl.scatter_w_clusterSize(results_subset, "onTarget", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix_subset)