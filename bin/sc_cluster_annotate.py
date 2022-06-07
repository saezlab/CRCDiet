from genericpath import sameopenfile
from operator import index
from imageio import save
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

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Run annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
parser.add_argument('-of', '--output_file', help='Output file name', required=False)
args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # "sc_cluster_annotate"
output_file = args['output_file']

# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
###############################

sample_type = "sc"
adata = sc.read_h5ad(input_path)

adata_concat = utils.get_filtered_concat_data(sample_type)
"""for item in adata_concat.var_names:
    if item.startswith("Pd"):
        print(item)"""
adata_concat = adata_concat[adata.obs_names,:]
# https://github.com/scverse/scanpy/issues/2239
# if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
# KeyError: 'base'
adata.uns['log1p']["base"] = None

l_param, _ = adata.uns["leiden_best_silh_param"]
l_param = f"{l_param:.2f}"


l_param_list = [0.30] # for major cell types
l_param_list = [0.20] # for SC data
step = 0.10
for l_param in l_param_list:
#for l_param in np.arange(0.1, 1.01, step):
    l_param = f"{l_param:.2f}"
    printmd(f"## Clusters with resolution param: {l_param} <a class='anchor' id='seventh-bullet-1'></a>")
    
    adata_concat.obs[f"leiden_{l_param}"] = adata.obs[f"leiden_{l_param}"]
    adata_concat.obsm["X_umap"] = adata.obsm["X_umap"]
    mpl.rcParams['figure.dpi']= 300
    mpl.rcParams["figure.figsize"] = (10,10)
    mpl.rcParams["legend.fontsize"]  = 'xx-small'
    mpl.rcParams["legend.loc"]  = "upper right"
    mpl.rcParams['axes.facecolor'] = "white"

    sc.pl.umap(adata_concat, color=f"leiden_{l_param}", palette=sc.pl.palettes.default_20, size=4 ,legend_loc='on data', show=True, save=f'{sample_type}_leiden_{l_param}_ondata')
    sc.pl.umap(adata_concat, color=f"leiden_{l_param}", palette=sc.pl.palettes.default_20, size=4 , show=False, save=f'{sample_type}_leiden_{l_param}_umap')

    # "Pdgrfra",
    markers_dot_plot = markers_dot_plot = ["Epcam", "Agr2", "Fabp2", "Krt14", "Pdgfra", "Myh11", "Ano1", "Lyve1", "Esam", "Ptprc", "Itgax", "Cd3g", "Mzb1", "Jchain", "Il17rb", "Cpa3", "S100a9", "Mki67"]
    sc.pl.dotplot(adata, markers_dot_plot, groupby=f'leiden_{l_param}', swap_axes=True, dendrogram=True,  show=True, save=f'{sample_type}_clusters_marker_{l_param}_dotplot_dendogram')
    sc.pl.dotplot(adata, markers_dot_plot, groupby=f'leiden_{l_param}', swap_axes=True, dendrogram=False,  show=True, save=f'{sample_type}_clusters_marker_{l_param}_dotplot')
    # change below anndata objects to "anndata" to run on only HVGs
    mpl.rcParams['figure.dpi']= 300
    mpl.rcParams["figure.figsize"] = (5,5)
    sc.tl.rank_genes_groups(adata_concat, groupby=f"leiden_{l_param}", method='wilcoxon', key_added = f"wilcoxon_{l_param}")
    mpl.rcParams['axes.titlesize'] = 20
    sc.pl.rank_genes_groups(adata_concat, n_genes=25, sharey=False, key=f"wilcoxon_{l_param}", show=True, groupby=f"leiden_{l_param}", save=f'{sample_type}_one_vs_rest_{l_param}')#
    # mpl.rcParams['axes.titlesize'] = 60
    # sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=f"wilcoxon_{l_param}", show=True, groupby=f"leiden_{l_param}", save=f'{sample_type}_deg_clusters_dotplot_{l_param}')

    wc = sc.get.rank_genes_groups_df(adata_concat, group=None, key=f"wilcoxon_{l_param}", pval_cutoff=0.01, log2fc_min=0)[["group", "names", "scores","logfoldchanges"]]
    # print(l_param)
    # print(wc.to_csv(os.path.join(output_path, f'{sample_type}_deg_leiden_res_{l_param}.csv'), index=False))
