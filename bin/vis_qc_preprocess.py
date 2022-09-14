from genericpath import sameopenfile
import os
import utils
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from utils import OUT_DATA_PATH, PLOT_PATH, DATA_PATH
import os
import plotting
from tabulate import tabulate
import warnings
from utils import printmd
from matplotlib import rcParams
import matplotlib as mpl


############################### BOOOORIING STUFF BELOW ###############################
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths("visium_qc_preprocess")
############################### BOOOORIING STUFF ABOVE ###############################


sample_type = "visium"
meta = utils.get_meta_data("visium")


def get_threshold_dict():
    """This functions keeps the threshold used to filter the data"""

    df_threshold = {"mt_thr": 10, # mitochondrial gene threshold
                # "rp_thr": 3, # ribosomal gene threshold
                "gene_thr": 200,
                "cell_thr": 5}

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
    mpl.rcParams["image.cmap"]= plt.cm.Spectral
    sc.pl.spatial(adata, img_key="hires", color = ["total_counts", "n_genes_by_counts",'pct_counts_mt', 'pct_counts_rp'],  size=1.25, alpha_img=0.5, wspace = 0.3, show=True, save=f"vis_{sample_id}.pdf")
    plt.show();

    fig, axs = plt.subplots(1, 5, figsize=(30, 10));
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, ax=axs[0])
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1])
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[2])
    plotting.plot_mt_vs_counts(adata, axs[3], mt_thr=df_threshold["mt_thr"])
    plotting.plot_rp_vs_counts(adata, axs[4]) # , rp_thr=df_threshold["rp_thr"])
    fig.savefig(os.path.join(PLOT_PATH, f"vis_basic_stats_before_filtering_{sample_id}.pdf"), dpi=300);
    plt.show();


    # number og genes at each change it to 300
    # each spot has at least 500 
    sc.pp.filter_cells(adata, min_genes=300)
    sc.pp.filter_cells(adata, min_counts=500)
    sc.pp.filter_genes(adata, min_cells=5)

    adata = adata[adata.obs.pct_counts_mt < df_threshold["mt_thr"], :]
    print("After mit filter ", np.shape(adata.X))
    # adata = adata[adata.obs.pct_counts_rp > df_threshold["rp_thr"], :]
    # print("After ribosomal filter ", np.shape(adata.X))

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
    
    sc.set_figure_params(figsize=(8, 8)) 
    print("Recalculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rp"], inplace=True)
    print("Plotting highest expressed genes after QC and filtering...")
    sc.pl.highest_expr_genes(adata, n_top=20, show=True, save=f"basic_stats_after_filtering_{sample_id}.pdf")
    
    adata.layers["raw"] = adata.X.copy()
    adata.layers["sqrt_norm"] = np.sqrt(
    sc.pp.normalize_total(adata, inplace=False)["X"]
)

    print("Saving filtered AnnData file...")
    adata.write(os.path.join(OUT_DATA_PATH, f"{sample_id}_filtered.h5ad"))

    return adata



def create_filtered_adata_files():

    for _, row in meta.iterrows():
    
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata = utils.read_raw_visium_sample(sample_id)
        printmd(f"<h4 style='color:black' align='center'>=============== Processing {condition} ===============")
        filter_cells_genes(adata, sample_id)
        # break

