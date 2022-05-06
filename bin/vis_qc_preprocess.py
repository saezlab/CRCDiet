from genericpath import sameopenfile
import os
import utils
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy.external as sce
from anndata._core.anndata import AnnData
from utils import OUT_DATA_PATH, PLOT_PATH, DATA_PATH, read_raw_sc_sample
import os
import argparse
import plotting
from tabulate import tabulate
import warnings
from utils import printmd
sc.settings.verbosity = 0



warnings.simplefilter(action='ignore')


meta = utils.get_meta_data("visium")


def get_threshold_dict():
    """This functions keeps the threshold used to filter the data"""

    df_threshold = {"mt_thr": 10, # mitochondrial gene threshold
                "rp_thr": 5, # ribosomal gene threshold
                "gene_thr": 200,
                "cell_thr": 3,
                "gene_qnt": 0.99}

    return df_threshold


def filter_cells_genes(adata, sample_id):
    """Perform basic filtering on cells and genes

    This function takes sample id as input and performs cell/gene filtering and remove mitochondrial genes and saves the AnnData for each sample <sample_id>_<condition>_filtered.h5ad

    Args:
        sample_id (str): the name of the folder where the sample files stored


    Returns:
        anndata: return the AnnData object

    """
    row = meta[meta["sample_id"]==sample_id]
    condition = str(row["condition"].values[0])

    # get threshold used to filter data
    df_threshold = get_threshold_dict()

    pre_filter_shape = np.shape(adata.X)

    print("Calculating QC metrics...")
    # calculate qc metrics
    adata.var["mt"] = adata.var_names.str.contains("^mt-")
    adata.var["rp"] = adata.var_names.str.contains("^Rp[sl]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rp"], inplace=True)

    sc.pl.spatial(adata, img_key="hires", title=condition, color="total_counts", color_map="bwr", size=1.0, alpha_img=0.5, wspace = 0.1, hspace = 1.0, show=True)
    sc.pl.spatial(adata, img_key="hires", title=condition, color="n_genes_by_counts", color_map="bwr", size=1.0, alpha_img=0.5, wspace = 0.1, hspace = 1.0, show=True)
    sc.pl.spatial(adata, img_key="hires", title=condition, color="pct_counts_rp", color_map="bwr", size=1.0, alpha_img=0.5, wspace = 0.1, hspace = 1.0, show=True)
    sc.pl.spatial(adata, img_key="hires", title=condition, color="pct_counts_mt", color_map="bwr", size=1.0, alpha_img=0.5, wspace = 0.1, hspace = 1.0, show=True)

    # , title=condition, color=adata.obs.pct_counts_mt, color_map="bwr", size=1.0, alpha_img=0.5, wspace = 0.1, hspace = 1.0, show=True)
    """
    adata = adata[adata.obs.pct_counts_mt < df_threshold["mt_thr"], :]
    adata = adata[adata.obs.pct_counts_rp > df_threshold["rp_thr"], :]
    adata = adata[adata.obs.doublet_score < df_threshold["doublet_thr"], : ]
    adata = adata[adata.obs.n_genes_by_counts < gene_quant_thr, : ]
    """
    # TODO: Visualize features
    # ST.FeaturePlot(se, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 2, grid.ncol = 1, cols = c("darkblue", "cyan", "yellow", "red", "darkred"), show.sb = F)

    # Plot before the preprocessing
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
    sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
    sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])
    fig.savefig(os.path.join(plot_path, "preprocess", f"basic_stats_before_filtering_{sample_id}.png") , dpi=300)

    # remove mitochondrial genes
    adata = adata[:,~adata.var["mt"]]
    
    # remove ribosomal genes
    adata = adata[:,~adata.var["rp"]]

    # number og genes at each change it to 300
    # each spot has at least 500 
    sc.pp.filter_cells(adata, min_genes=300)
    sc.pp.filter_cells(adata, min_counts=500)
    sc.pp.filter_genes(adata, min_cells=5)

    # Plot after the preprocessing
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
    sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
    sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])
    fig.savefig(os.path.join(plot_path, "preprocess", f"basic_stats_after_filtering_{sample_id}.png") , dpi=300)

    post_filter_shape = np.shape(adata.X)
    adata.obs["batch"] = sample_id
    print(f"{sample_id}:\nAnnData shape before filtering {pre_filter_shape}")
    print(f"AnnData shape before filtering {post_filter_shape}")
    adata.write(os.path.join(out_data_path, f"{sample_id}_{condition}_filtered.h5ad"))

    return adata



def create_filtered_adata_files():

    for _, row in meta.iterrows():
    
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata = utils.read_raw_visium_sample(sample_id)
        printmd(f"<h4 style='color:black' align='center'>=============== Processing {condition} ===============")
        print(f"")
        
        filter_cells_genes(adata, sample_id)
        # break


create_filtered_adata_files()