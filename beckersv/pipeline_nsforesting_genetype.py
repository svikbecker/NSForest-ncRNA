"""
Author: Sasha "Vik" Becker, Graduate Student, George Washington University
PI: Dr. Richard Scheuermann, National Library of Medicine, NIH
Date: Summer 2026 

This script contains the code for running a comparative 
NS-Forest based analysis of a full dataset and a specificed gene type subset.

The steps to run this script are as follows:
1. Initialize a NSForestComparison object with the appropriate file paths.
2. Set up the desired data using the set_data method.
3. Set the gene type(s) subset using the set_genetype_subset method.
4. Run the NS-Forest analysis on both the full dataset and the gene type(s) subset using
   the run_comparison method.

"""

# SETUP
import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import plotly.io as pio
pio.renderers.default = "notebook"
import nsforest as ns
from nsforest import utils

class NSForestComparison:
    def __init__(self, code_folder: str, data_folder: str, output_folder: str):
        self.code_folder = code_folder
        self.data_folder = data_folder
        self.output_folder = output_folder
        self.dataset_file = None # to be set in set_data
        self.adata = None # to be set in set_data
        self.annotation_file = None # to be set in set_data
        self.cluster_header = None # to be set in set_data
        self.genetype_subset = None # to be set in set_genetype_subset
        return
        
    def set_data(self, dataset_file: str, annotation_file: str = "", cluster_header: str = "",
                 print_info: bool = True):
        
        """
        Set up the data for analysis by loading the specified dataset file.

        Parameters
        ----------
        dataset_file : str
            The relative file path of the dataset to be loaded (e.g., 'dataset.h5ad').
            Supports all file formats compatible with Scanpy's read function 
            (e.g., .h5ad, .loom, .csv, etc.).
            
        annotation_file : str
            The relative file path of the gene type annotation file (e.g., 'gencode_annotation.csv').
            
        cluster_header : str
            The name of the column in the dataset that contains the cluster labels.
            
        print_info : bool, optional
            If True, prints information about the loaded dataset. Default is True.
        
        Returns
        -------
        None

        """
        
        # load the dataset
        self.dataset_file = dataset_file
        self.adata = sc.read(os.path.join(self.data_folder, self.dataset_file))
        
        if print_info:
            print(f"Dataset '{self.dataset_file}' loaded successfully. Here are some details:")
            print(f"Number of cells: {self.adata.n_obs}")
            print(f"Number of genes: {self.adata.n_vars}")
            print(f"Available columns in adata.obs: {self.adata.obs.columns.tolist()}")
            print(f"Available columns in adata.var: {self.adata.var.columns.tolist()}")
            
        # if annotation_file isn't empty, load the gene type information from gencode annotation
        if annotation_file not in ["", None]:
            self.annotation_file = annotation_file
            gencode = pd.read_csv(os.path.join(self.data_folder, self.annotation_file))
            
            # set gene_name as the index for joining
            gencode = gencode.set_index("gene_name")

            # add annotations to adata.var
            self.adata.var = self.adata.var.join(gencode, how="left")

            # if print_info, check result
            if print_info:
                print("Gene type annotation added to adata.var. Here are the first few rows:")
                print(self.adata.var.head())
        
        return
        