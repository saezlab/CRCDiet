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


sc.settings.verbosity = 1

warnings.simplefilter(action='ignore')


S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
PLOT_PATH =  os.path.join(S_PATH, "../plots", "qc_preprocess")

Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = PLOT_PATH

meta = utils.get_meta_data("sc")


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

    # adata = utils.read_raw_sample(condition)
    
    pre_filter_shape = np.shape(adata.X)

    doublet_thr = 0.2
    mt_thr = 10 # mitochondrial gene threshold
    rp_thr = 5 # ribosomal gene threshold
    gene_qnt = 0.99

    print("Calculating QC metrics...")
    # calculate qc metrics
    adata.var["mt"] = adata.var_names.str.contains("^mt-")
    adata.var["rp"] = adata.var_names.str.contains("^Rp[sl]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rp"], inplace=True)

    # calculate doublet scores
    sce.pp.scrublet(adata, verbose=False)
    gene_thr = np.quantile(adata.obs.n_genes_by_counts, gene_qnt)
    
    #  “total_counts”. Sum of counts for a gene.
    #  “n_genes_by_counts”. The number of genes with at least 1 count in a cell. Calculated for all cells.
    plt.figure();
    fig, axs = plt.subplots(2, 4, figsize=(30, 10))# , figsize=(100, 20))
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, ax=axs[0][0])
    
    # sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=False, ax=axs[0][1])
    plotting.plot_mt_vs_counts(adata, axs[0][1], mt_thr=mt_thr)
    plotting.plot_ngenes_vs_counts(adata, axs[0][2], gene_thr=gene_thr)

    # sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', legend_fontsize="xx-large", show=False, ax=axs[0][2])
    
    plotting.plot_doublet_scores(adata, axs[0][3], doublet_thr=doublet_thr, fontsize=11)
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[1][0])
    sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1][1])
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1][2])
    sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[1][3])
    plt.show();
    # plt.show()
    # fig.savefig(os.path.join(PLOT_PATH, f"basic_stats_before_filtering_{condition}.png"), dpi=300)
    # plt.clf()
    
    
    # adata.obs.doublet_score < doublet_thr
    # number og genes at each change it to 200 
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    # filter based on total counts
    gene_thr = np.quantile(adata.obs.n_genes_by_counts, gene_qnt)
    adata = adata[adata.obs.pct_counts_mt < mt_thr, :]
    adata = adata[adata.obs.pct_counts_rp > rp_thr, :]
    adata = adata[adata.obs.doublet_score < doublet_thr, : ]
    adata = adata[adata.obs.n_genes_by_counts < gene_thr, : ]
    post_filter_shape = np.shape(adata.X)
    
    
    # Assume they are coming from the same batch
    adata.obs["batch"] = 0
    adata.obs["condition"] = condition
    
    
    print(tabulate([[condition, "Before filtering", pre_filter_shape[0], pre_filter_shape[1]], [condition, "After filtering", post_filter_shape[0], post_filter_shape[1]]], headers=["Sample ID", 'Stage', "# of cells", "# of genes"], tablefmt='fancy_grid'))
    # print( {pre_filter_shape}")
    # print(f"AnnData shape after filtering {post_filter_shape}")

    print("Recalculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rp"], inplace=True)
    plt.figure();
    fig, axs = plt.subplots(2, 4, figsize=(30, 10))# , figsize=(100, 20))
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, ax=axs[0][0])
    
    # sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=False, ax=axs[0][1])
    plotting.plot_mt_vs_counts(adata, axs[0][1], mt_thr=mt_thr)
    plotting.plot_ngenes_vs_counts(adata, axs[0][2], gene_thr=gene_thr)

    # sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', legend_fontsize="xx-large", show=False, ax=axs[0][2])
    
    plotting.plot_doublet_scores(adata, axs[0][3], doublet_thr=doublet_thr, fontsize=11)
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[1][0])
    sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1][1])
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1][2])
    sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[1][3])
    plt.show()
    plt.clf();
    # fig.savefig(os.path.join(PLOT_PATH, f"basic_stats_after_filtering_{sample_id}.png"), dpi=300)

    """del adata.obs["predicted_doublet"]
    adata.write(os.path.join(OUT_DATA_PATH, f"{sample_id}_{condition}_filtered.h5ad"))
    # print(os.path.join(out_data_path, f"{sample_id}_{condition}_filtered.h5ad"))"""
    

    return adata

def create_filtered_adata_files():

    for _, row in meta.iterrows():
    
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata = utils.read_raw_sc_sample(sample_id)
        printmd(f"<h3 style='color:grey'>=============== Processing {sample_id} ===============</h1>")
        print(f"")
        
        filter_cells_genes(adata, sample_id)
        # break



def get_processed_sample_from_adata_file(sample_id):
    """Given samples id get filtered adata object

    This function takes sample id as input and returns the filtered AnnData object

    Args:
        sample_id (str): the name of the folder where the sample files stored
    
    Returns:
        Filtered AnnData object of the sample

    """
    # row = meta[meta["sample_id"]==sample_id]
    adata = sc.read(os.path.join(OUT_DATA_PATH, f"{sample_id}_filtered.h5ad"))

    return adata




    


