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
from copy import deepcopy

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
sc.settings.verbosity = 0
warnings.simplefilter(action='ignore')
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Parse arguments
parser = argparse.ArgumentParser(prog='Progeny', description='Pathway activity estimation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # "sc_pse_func_analys"  


# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 

sample_type = "visium"
meta = utils.get_meta_data(sample_type)
adata_integ_clust = sc.read_h5ad(input_path)


# Retrieve PROGENy model weights
progeny = dc.get_progeny(organism='mouse', top=500)
progeny["target"] = progeny["target"].str.upper()
adata_integ_clust.X = adata_integ_clust.layers['log1p_transformed']
adata_integ_clust.var.index = pd.Index(gen.upper() for gen in adata_integ_clust.var.index.values)
dc.run_mlm(mat=adata_integ_clust, net=progeny, source='source', target='target', weight='weight', verbose=True, use_raw=False)
adata_integ_clust.obsm['dorothea_mlm_estimate'] = adata_integ_clust.obsm['mlm_estimate'].copy()
adata_integ_clust.obsm['dorothea_mlm_pvals'] = adata_integ_clust.obsm['mlm_pvals'].copy()

### This block is for the comparative pathway activity estimation 
"""adata_integ_clust = adata_integ_clust[adata_integ_clust.obs["leiden_0.10"]=="6",:]
adata_integ_clust = adata_integ_clust[adata_integ_clust.obs["condition"].isin(["CD-AOM-DSS", "HFD-AOM-DSS", "LFD-AOM-DSS"]),:]
acts_integrated = dc.get_acts(adata_integ_clust, obsm_key='dorothea_mlm_estimate')
mean_acts = dc.summarize_acts(acts_integrated, groupby="condition", min_std=0)
sns.clustermap(mean_acts, xticklabels=mean_acts.columns, vmin=-2, vmax=2, colors_ratio=0.50, cmap='coolwarm')
plt.savefig(f"{PLOT_PATH}/{sample_type}_aom_dss_only_tumor_cluster_pathway_activity_est_cmap.pdf")
plt.show();"""
### This block is for the comparative pathway activity estimation 

adata_dict = dict()


for _, row in meta.iterrows():

    sample_id = row["sample_id"]
    condition = row["condition"]

    # printmd(f"<h4 style='color:black' align='center'>=============== {condition} ===============")

    adata_filtered = sc.read_h5ad(os.path.join(OUT_DATA_PATH, f"{sample_id}_filtered.h5ad"))
    adata_filtered.X = adata_filtered.layers['log1p_transformed']
    adata_filtered.var.index = pd.Index(gen.upper() for gen in adata_filtered.var.index.values)
    adata_filtered = adata_filtered[:,adata_filtered.var_names.intersection(adata_integ_clust.var_names)]

    adata_temp = adata_integ_clust[adata_integ_clust.obs["condition"]==condition,:].copy()
    cols = adata_integ_clust[adata_integ_clust.obs["condition"]==condition,:].obsm['mlm_estimate'].columns
    mlm_estimate_df =  pd.DataFrame(adata_integ_clust[adata_integ_clust.obs["condition"]==condition,:].obsm['mlm_estimate'].values, index=adata_filtered.obs_names, columns=cols)
    adata_filtered.obsm['progeny_mlm_estimate'] = mlm_estimate_df #adata_temp.obsm['mlm_estimate'].copy()
    adata_dict[sample_id] = adata_filtered
    
    """mpl.rcParams["image.cmap"]= plt.cm.Spectral
    fig, ax = plt.subplots(2, int(len(adata_filtered.obsm['progeny_mlm_estimate'].columns)/2), figsize=(60,15))

    for ind, col in enumerate(adata_filtered.obsm['progeny_mlm_estimate'].columns):
        fig_row, fig_col = int(ind/7), ind%7
        mpl.rcParams["image.cmap"]= plt.cm.magma
        mpl.rcParams['axes.titlesize'] = 15
        sc.pl.spatial(acts, color=col, size=1.25, alpha_img=0.5, wspace = 2.5, hspace=1.5, ax = ax[fig_row][fig_col], show=False)
        cbar = ax[fig_row][fig_col].collections[0].colorbar
        cbar.set_ticks([])
        cbar = None
    
    plt.savefig(f'{PLOT_PATH}/{sample_id}_pathway.pdf');"""
    # plt.show();
    # plt.clf();
    

    # sc.pl.umap(acts, color=adata_filtered.obsm['progeny_mlm_estimate'].columns, vcenter=0, cmap='coolwarm', show=True) # , save=f'{sample_type}_pathway_activity_est')

    # mean_acts = dc.summarize_acts(acts, groupby='condition', min_std=0)
    # adata_samp_cc = adata_cc_merged[adata_cc_merged.obs["condition"]==condition,:]
    # sns.clustermap(mean_acts, xticklabels=mean_acts.columns, vmin=-2, vmax=2, cmap='coolwarm')
    # plt.savefig(f"{PLOT_PATH}/{sample_type}_cell_type_pathway_activity_est_cmap.pdf")
    # plt.show();
    # , save=f'{sample_type}_pathway_activity_est_clustermap'
print(len(meta))
for ind, col in enumerate(adata_filtered.obsm['progeny_mlm_estimate'].columns):

    sc.pl.umap(acts_integrated, color=[col, 'leiden_0.10'], show=False, cmap='coolwarm', vcenter=0, save=f"{col}.pdf")

    fig, ax = plt.subplots(1, len(meta), figsize=(40,5))
    print(ind, col)
    for ind, row in meta.iterrows():
        sample_id = row["sample_id"]
        condition = row["condition"]
        mpl.rcParams["image.cmap"]= plt.cm.magma_r
        acts = dc.get_acts(adata_dict[sample_id], obsm_key='progeny_mlm_estimate')
        sc.pl.spatial(acts, color=col, size=1.25, title=condition, alpha_img=0.5, wspace = 2.5, hspace=1.5, ax=ax[ind], show=False)
        cbar = ax[ind].collections[0].colorbar
        cbar.set_ticks([])
        cbar = None
        
    plt.savefig(f'{PLOT_PATH}/{col}_pathway_activity.pdf');
    # plt.show();
    # plt.clf();

# python vis_downstream_analysis.py -i ../data/out_data/visium_integrated_clustered.h5ad  -o ../data/out_data/ -an visium_pathway_act_est