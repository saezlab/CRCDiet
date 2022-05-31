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

'''
Run QC and filtering per sample
'''

############################### BOOOORIING STUFF BELOW ###############################
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=150)
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths("sc_qc_preprocess")
############################### BOOOORIING STUFF ABOVE ###############################

meta = utils.get_meta_data("sc")

def get_threshold_dict():
    """This functions keeps the threshold used to filter the data
    
    Returns:
        dict: return the dictionary of used thresholds
    """

    threshold_dict = {"mt_thr": 20, # mitochondrial gene threshold
                "rp_thr": 5, # ribosomal gene threshold
                "doublet_thr": 0.2, #doublet threshold
                "gene_thr": 200,
                "cell_thr": 3,
                "gene_qnt": 0.99}

    return threshold_dict



def filter_cells_genes(adata, sample_id):
    """Perform basic filtering on cells and genes

    This function takes sample id as input and performs cell/gene filtering and remove mitochondrial genes\
         and saves the AnnData for each sample <sample_id>_<condition>_filtered.h5ad

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

    # calculate doublet scores
    sce.pp.scrublet(adata, verbose=False)
    # calculate threshold to remove extreme outliers
    gene_quant_thr = np.quantile(adata.obs.n_genes_by_counts, df_threshold["gene_qnt"])
    
    #  “total_counts”. Sum of counts for a gene.
    #  “n_genes_by_counts”. The number of genes with at least 1 count in a cell. Calculated for all cells.
    # plt.figure();
    fig, axs = plt.subplots(2, 4, figsize=(30, 10));
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, ax=axs[0][0])
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[0][1])
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0][2])
    sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[0][3])

    plotting.plot_mt_vs_counts(adata, axs[1][0], mt_thr=df_threshold["mt_thr"])
    plotting.plot_rp_vs_counts(adata, axs[1][1], rp_thr=df_threshold["rp_thr"])
    plotting.plot_doublet_scores(adata, axs[1][2], doublet_thr=df_threshold["doublet_thr"], fontsize=11)
    plotting.plot_ngenes_vs_counts(adata, axs[1][3], gene_thr=gene_quant_thr)
    fig.savefig(os.path.join(PLOT_PATH, f"basic_stats_before_filtering_{sample_id}.png"));
    plt.show();
    # plt.clf();
    

    # number og genes at each change it to 200 
    sc.pp.filter_cells(adata, min_genes = df_threshold["gene_thr"])
    sc.pp.filter_genes(adata, min_cells=df_threshold["cell_thr"])

    # filter based on total counts    
    adata = adata[adata.obs.pct_counts_mt < df_threshold["mt_thr"], :]
    adata = adata[adata.obs.pct_counts_rp > df_threshold["rp_thr"], :]
    adata = adata[adata.obs.doublet_score < df_threshold["doublet_thr"], : ]
    adata = adata[adata.obs.n_genes_by_counts < gene_quant_thr, : ]

    # Filter MALAT1 and Gm42418
    adata = adata[:, ~adata.var_names.str.startswith('Malat1')]
    adata = adata[:, ~adata.var_names.str.startswith('Gm42418')]
    
    # remove mitochondrial genes
    adata = adata[:,~adata.var["mt"]]
    
    # remove ribosomal genes
    adata = adata[:,~adata.var["rp"]]

    post_filter_shape = np.shape(adata.X)

    adata.obs["condition"] = condition    
    
    print(tabulate([[condition, "Before filtering", pre_filter_shape[0], pre_filter_shape[1]],\
                    [condition, "After filtering", post_filter_shape[0], post_filter_shape[1]]],\
                    headers=["Sample ID", 'Stage', "# of cells", "# of genes"], tablefmt='fancy_grid'))
    # print( {pre_filter_shape}")
    # print(f"AnnData shape after filtering {post_filter_shape}")
    sc.set_figure_params(figsize=(8, 8)) 
    print("Recalculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rp"], inplace=True)
    print("Plotting highest expressed genes after QC and filtering...")
    sc.pl.highest_expr_genes(adata, n_top=20, show=True, save=f"basic_stats_after_filtering_{sample_id}.png")

    adata.layers["raw"] = adata.X.copy()
    adata.layers["sqrt_norm"] = np.sqrt(
    sc.pp.normalize_total(adata, inplace=False)["X"])
    
    del adata.obs["predicted_doublet"]
    print("Saving filtered AnnData file...")
    adata.write(os.path.join(OUT_DATA_PATH, f"{sample_id}_filtered.h5ad"))
    
    return adata

def create_filtered_adata_files():

    for _, row in meta.iterrows():
    
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata = utils.read_raw_sc_sample(sample_id)
        printmd(f"<h4 style='color:black' align='center'>=============== Processing {condition} ===============")
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




    


