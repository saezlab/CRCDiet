from genericpath import sameopenfile
from operator import index
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns 
import decoupler as dc  
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import os
import sys
from sklearn.metrics import silhouette_score, pairwise_distances
import sys
import warnings
import utils
from utils import printmd
import matplotlib as mpl

sc.settings.verbosity = 0


warnings.simplefilter(action='ignore')

# Read integrated object
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Run annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']

###############################

sample_type = "sc"

S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
PLOT_PATH =  os.path.join(S_PATH, "../plots", "sc_annotate")

Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = PLOT_PATH
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80)

adata = sc.read_h5ad(input_path)

adata_concat = utils.get_filtered_concat_data(sample_type)
# https://github.com/scverse/scanpy/issues/2239
# if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
# KeyError: 'base'
adata.uns['log1p']["base"] = None

l_param, _ = adata.uns["leiden_best_silh_param"]
l_param = f"{l_param:.2f}"
sample_type = "sc"

l_param_list = [0.20]
for l_param in l_param_list:
    l_param = f"{l_param:.2f}"
    printmd(f"## Clusters with resolution param: {l_param} <a class='anchor' id='seventh-bullet-1'></a>")
    
    adata_concat.obs[f"leiden_{l_param}"] = adata.obs[f"leiden_{l_param}"]
    adata_concat.obsm["X_umap"] = adata.obsm["X_umap"]
    mpl.rcParams["figure.figsize"] = (10,10)
    sc.pl.umap(adata_concat, color=f"leiden_{l_param}", palette=sc.pl.palettes.default_20, size=4 , show=True);
    # plt.show();
    # change below anndata objects to "anndata" to run on only HVGs
    mpl.rcParams["figure.figsize"] = (5,5)
    sc.tl.rank_genes_groups(adata_concat, groupby=f"leiden_{l_param}", method='wilcoxon', key_added = f"wilcoxon_{l_param}")
    mpl.rcParams['axes.titlesize'] = 20
    sc.pl.rank_genes_groups(adata_concat, n_genes=25, sharey=False, key=f"wilcoxon_{l_param}", show=True, groupby=f"leiden_{l_param}", save=f'{sample_type}_one_vs_rest_{l_param}')#
    # mpl.rcParams['axes.titlesize'] = 60
    # sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=f"wilcoxon_{l_param}", show=True, groupby=f"leiden_{l_param}", save=f'{sample_type}_deg_clusters_dotplot_{l_param}')

    wc = sc.get.rank_genes_groups_df(adata_concat, group=None, key=f"wilcoxon_{l_param}", pval_cutoff=0.01, log2fc_min=0)[["group", "names", "scores","logfoldchanges"]]
    # print(l_param)
    # print(wc.to_csv(os.path.join(output_path, f'{sample_type}_deg_leiden_res_{l_param}.csv'), index=False))
    