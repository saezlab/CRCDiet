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
from utils import printmd
import anndata
import utils
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

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # "sc_cluster_annotate"


# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
###############################

sample_type = "visium"
meta = utils.get_meta_data(sample_type)
adata = sc.read_h5ad(input_path)

adata.raw = anndata.AnnData(adata.layers['counts'], obs=adata.obs, var=adata.var)
# print(adata.layers['counts'])
# https://github.com/scverse/scanpy/issues/2239
# if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
# KeyError: 'base'
# adata.uns['log1p']["base"] = None

# Reasanable clusters after visual inspection
# 0.10, 0.40, 0.60

# This is not used, because I wanted to provide a few clustering results based on different res parameter
l_param, _ = adata.uns["leiden_best_silh_param"]

l_param_list = [0.10, 0.40, 0.60]
l_param_list = [0.40]
for l_param in l_param_list:
    l_param = f"{l_param:.2f}"
    printmd(f"## Clusters with resolution param: {l_param} <a class='anchor' id='seventh-bullet-1'></a>")

    rows, cols = (1, 6)
    fig, ax = plt.subplots(rows, cols, figsize=(20,20))

    for ind, row in meta.iterrows():
        
        fig_row, fig_col = int(ind/cols), ind%cols
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata_raw = utils.read_raw_visium_sample(sample_id)
        adata_temp = adata[adata.obs["condition"]==condition,:]    
        adata_temp.obs.index = pd.Index("-".join(cl.split("-")[:-1]) for cl in adata_temp.obs.index.values)
        adata_raw = adata_raw[adata_temp.obs.index,:]
        adata_raw = adata_raw[:, list(adata.var_names)]
        adata_raw.obs[f"leiden_{l_param}"] = adata_temp.obs[f"leiden_{l_param}"]
        mpl.rcParams['axes.titlesize'] = 12
        sc.pl.spatial(adata_raw, img_key="hires", color =f"leiden_{l_param}",  title=condition, size=1.25, alpha_img=0.3, ax = ax[ind], show=False)
    fig.tight_layout()
    plt.savefig(f"{PLOT_PATH}/spatial_cluster.pdf")
    plt.show();
    
    

    # TODO: Use whole transcriptome instead of HVGs
    # Comment out the section below for running DEGs on HVGs
    # TODO: Refactor this script, it is ugly and inefficient
    adata_concat = []
    for ind, row in meta.iterrows():
        
        fig_row, fig_col = int(ind/cols), ind%cols
        sample_id = row["sample_id"]
        condition = row["condition"]
        # print(sample_id)
        tmp = sc.read_h5ad(os.path.join(OUT_DATA_PATH,f"{sample_id}_filtered.h5ad"))
        # Fetch sample metadata
        m = meta[meta['sample_id'] == sample_id]
        # Add metadata to adata
        for col in m.columns:
            tmp.obs[col] = m[col].values[0]
        # Append
        adata_concat.append(tmp)
        del tmp
    
    # Merge objects and delete list
    adata_concat = adata_concat[0].concatenate(adata_concat[1:], join='outer')

    adata_concat.obs[f"leiden_{l_param}"] = adata.obs[f"leiden_{l_param}"]
    # Comment out the section above for running DEGs on HVGs
    

    
    # change below anndata objects to "anndata" to run on only HVGs
    sc.tl.rank_genes_groups(adata_concat, groupby=f"leiden_{l_param}", method='wilcoxon', key_added = f"wilcoxon_{l_param}")
    mpl.rcParams['axes.titlesize'] = 20
    sc.pl.rank_genes_groups(adata_concat, n_genes=25, sharey=False, key=f"wilcoxon_{l_param}", show=True, groupby=f"leiden_{l_param}", save=f'{sample_type}_one_vs_rest_{l_param}.pdf')#
    # mpl.rcParams['axes.titlesize'] = 60
    # sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=f"wilcoxon_{l_param}", show=True, groupby=f"leiden_{l_param}", save=f'{sample_type}_deg_clusters_dotplot_{l_param}')

    wc = sc.get.rank_genes_groups_df(adata_concat, group=None, key=f"wilcoxon_{l_param}", pval_cutoff=0.01, log2fc_min=0)[["group", "names", "scores","logfoldchanges"]]
    # print(l_param)
    # print(wc.to_csv(os.path.join(output_path, f'{sample_type}_deg_leiden_res_{l_param}.csv'), index=False))

"""marker_genes = None
own_markers = True
marker_db="PanglaoDB"

if own_markers:
    # adata = adata.raw.to_adata()
    # TODO: make this parametric
    marker_genes = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
    marker_genes["genesymbol"] = marker_genes["genesymbol"].str.capitalize()
    marker_genes['weight'] = 1
    marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]

else:
    
    # print(marker_genes)
    # marker_genes = dc.get_resource('PanglaoDB')
    # marker_db="PanglaoDB"
    marker_genes = pd.read_csv(os.path.join(DATA_PATH,"PanglaoDB_markers_27_Mar_2020.tsv"), sep="\t")
    print(marker_genes)
    # Filter by canonical_marker and human
    
    marker_genes = marker_genes[marker_genes["species_official"].str.contains("mouse") & marker_genes["canonical_marker"]==1.0]
    
    print(marker_genes["genesymbol"])
    # Remove duplicated entries
    marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]

dc.run_ora(mat=adata, net=marker_genes, source='cell_type', target='genesymbol', min_n=3, verbose=True)
acts = dc.get_acts(adata, obsm_key='ora_estimate')

# print("============= ORA p-Vals ============= ")
# print(adata.obsm["ora_pvals"])
# print("============= ora_estimate ============= ")
# print(adata.obsm["ora_estimate"])
for c_t in list(adata.obsm["ora_estimate"].columns):
    adata.obs[c_t] = adata.obsm["ora_estimate"][c_t]

for c_t in list(adata.obsm["ora_estimate"].columns):
    rows, cols = (1, 6)
    fig, ax = plt.subplots(rows, cols, figsize=(20,2))

    for ind, row in meta.iterrows():
        
        # fig_row, fig_col = int(ind/cols), ind%cols
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata_raw = utils.read_raw_visium_sample(sample_id)
        adata_temp = adata[adata.obs["condition"]==condition,:]    
        adata_temp.obs.index = pd.Index("-".join(cl.split("-")[:-1]) for cl in adata_temp.obs.index.values)
        adata_raw = adata_raw[adata_temp.obs.index,:]
        adata_raw = adata_raw[:, list(adata.var_names)]
        adata_raw.obs[c_t] = adata_temp.obs[c_t]
        adata_raw.obs[f"leiden_{l_param}"] = adata_temp.obs[f"leiden_{l_param}"]
        mpl.rcParams["image.cmap"]= plt.cm.magma
        mpl.rcParams['axes.titlesize'] = 14
        sc.pl.spatial(adata_raw, img_key="hires", color =c_t, title=f"{c_t}_{condition}", size=1.25, alpha_img=0.5, ax = ax[ind], show=False) # ax[fig_row][fig_col]
        cbar = ax[ind].collections[0].colorbar
        cbar.set_ticks([])
        # plt.tight_layout()
    # plt.tight_layout()
    plt.show();

mpl.rcParams['axes.titlesize'] = 20
cell_types = list(set(marker_genes["cell_type"]))
sc.pl.umap(adata, color=list(adata.obsm["ora_pvals"].columns))

dict_mean_enr = dict()

mean_enr = dc.summarize_acts(acts, groupby=f'leiden_{l_param}')
"""
# python vis_cluster_annotate.py -i ../data/out_data/visium_integrated.h5ad -o ../data/out_data
