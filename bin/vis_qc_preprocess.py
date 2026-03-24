import os
import utils
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from utils import OUT_DATA_PATH, PLOT_PATH, DATA_PATH
import plotting
from tabulate import tabulate
import warnings
import argparse
from utils import printmd
import matplotlib as mpl


############################### BOOOORIING STUFF BELOW ###############################
# Warning settings
warnings.simplefilter(action='ignore')
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
sc.settings.verbosity = 0

parser = argparse.ArgumentParser(prog='visium qc preprocess', description='Run QC and preprocess...')
# Get necesary paths and create folders if necessary
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)

args = vars(parser.parse_args())
sample_type = args['sample_type'] # visium, visium_helminth
analysis_name = args['analysis_name'] # sc_integrate
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name) # visium_qc_preprocess, visium_helminth_qc_preprocess

############################### BOOOORIING STUFF ABOVE ###############################

# Read meta data
meta = utils.get_meta_data(sample_type)
print(meta)

def get_threshold_dict():
    """This functions keeps the threshold used to filter the data"""

    df_threshold = {"mt_thr": 10, # mitochondrial gene threshold
                "gene_thr": 300,
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
    adata.var_names = [vn.upper() for vn in adata.var_names]
    df_empty_spots = pd.read_csv(f"../data/visium_outputs/{sample_id}/EmptySpots_{sample_id}.csv")
    adata = adata[~adata.obs_names.isin(df_empty_spots["Barcode"]),:]
    print(pre_filter_shape, adata.shape)
    print("Calculating QC metrics...")
    # calculate qc metrics
    adata.var["mt"] = adata.var_names.str.lower().str.contains("^mt-")
    adata.var["rp"] = adata.var_names.str.lower().str.contains("^rp[sl]")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rp"], inplace=True)
    # mpl.rcParams["image.cmap"]= plt.cm.Spectral
    #sc.pl.spatial(adata, img_key="hires", cmap = "Spectral", color = ["total_counts", "n_genes_by_counts",'pct_counts_mt', 'pct_counts_rp'],  size=1.25, alpha_img=0.5, wspace = 1.0, hspace = 1.0, show=True, save=f"vis_{sample_id}.pdf")
    #plt.show();

    fig, axs = plt.subplots(1, 5, figsize=(30, 10));
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, ax=axs[0])
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1])
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[2])
    plotting.plot_mt_vs_counts(adata, axs[3], mt_thr=df_threshold["mt_thr"])
    plotting.plot_rp_vs_counts(adata, axs[4]) # , rp_thr=df_threshold["rp_thr"])
    fig.savefig(os.path.join(PLOT_PATH, f"vis_basic_stats_before_filtering_{sample_id}.pdf"), dpi=300);
    plt.show();


    sc.pp.filter_cells(adata, min_genes=df_threshold["gene_thr"])
    sc.pp.filter_cells(adata, min_counts=500)
    sc.pp.filter_genes(adata, min_cells=df_threshold["cell_thr"])

    adata = adata[adata.obs.pct_counts_mt < df_threshold["mt_thr"], :]
    print("After mit filter ", np.shape(adata.X))
    # adata = adata[adata.obs.pct_counts_rp > df_threshold["rp_thr"], :]
    # print("After ribosomal filter ", np.shape(adata.X))

    # Filter MALAT1 and Gm42418
    adata = adata[:, ~adata.var_names.str.lower().str.startswith('malat1')]
    adata = adata[:, ~adata.var_names.str.lower().str.startswith('gm42418')]
    adata = adata[:, ~adata.var_names.str.lower().str.startswith('mtrp')]
    
    # remove mitochondrial genes
    adata = adata[:,~adata.var["mt"]]
    # remove ribosomal genes
    adata = adata[:,~adata.var["rp"]]

    post_filter_shape = np.shape(adata.X)
    adata.obs["condition"] = condition   
    adata.obs["sample"] = sample_id      
    
    print(tabulate([[condition, "Before filtering", pre_filter_shape[0], pre_filter_shape[1]],\
                    [condition, "After filtering", post_filter_shape[0], post_filter_shape[1]]],\
                    headers=["Sample ID", 'Stage', "# of cells", "# of genes"], tablefmt='fancy_grid'))
    
    sc.set_figure_params(figsize=(8, 8)) 
    print("Recalculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rp"], inplace=True)
    print("Plotting highest expressed genes after QC and filtering...")
    sc.pl.highest_expr_genes(adata, n_top=20, show=True, save=f"basic_stats_after_filtering_{sample_id}.pdf")
    
    adata.layers["counts"] = adata.X.copy()

    adata.layers["log1p_transformed"] = sc.pp.log1p(
    sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)["X"])

    print("Saving filtered AnnData file...")
    adata.write(os.path.join(OUT_DATA_PATH, f"{sample_id}_filtered.h5ad"))

    return adata



def create_filtered_adata_files(raw=True):

    for _, row in meta.iterrows():
    
        sample_id = row["sample_id"]
        condition = row["condition"]
        if raw:
            adata = utils.read_raw_visium_sample(sample_id)
        else:
            adata = sc.read_h5ad(f"/net/data.isilon/ag-saez/bq_arifaioglu/home/Projects/CRCDiet/data/out_data/{sample_id}.h5ad")
        printmd(f"<h4 style='color:black' align='center'>=============== Processing {condition} ===============")
        filter_cells_genes(adata, sample_id)

create_filtered_adata_files(raw=True)

# python vis_qc_preprocess.py -an vis_microbiota_qc_preprocess -st visium_microbiota
# python vis_qc_preprocess.py -an visium_qc_preprocess -st visium