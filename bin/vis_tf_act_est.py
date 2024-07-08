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
parser = argparse.ArgumentParser(prog='Progeny', description='PathTFway activity estimation')
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
net = dc.get_collectri(organism='mouse', split_complexes=False)
net["target"] = net["target"].str.upper()
net["source"] = net["source"].str.upper()

# print(adata_integ_clust.var_names)
# print(net)

adata_integ_clust.X = adata_integ_clust.layers['log1p_transformed']
# adata_integ_clust.var.index = pd.Index(gen.upper() for gen in adata_integ_clust.var.index.values)
dc.run_mlm(mat=adata_integ_clust, net=net, source='source', target='target', weight='weight', verbose=True, use_raw=False)
adata_integ_clust.obsm['dorothea_mlm_estimate'] = adata_integ_clust.obsm['mlm_estimate'].copy()
adata_integ_clust.obsm['dorothea_mlm_pvals'] = adata_integ_clust.obsm['mlm_pvals'].copy()
acts_integrated = dc.get_acts(adata_integ_clust, obsm_key='dorothea_mlm_estimate')
# print(adata_integ_clust.obsm['dorothea_mlm_estimate'])


adata_dict = dict()

for _, row in meta.iterrows():

    sample_id = row["sample_id"]
    condition = row["condition"]

    # printmd(f"<h4 style='color:black' align='center'>=============== {condition} ===============")

    adata_filtered = sc.read_h5ad(os.path.join(OUT_DATA_PATH, f"{sample_id}_filtered.h5ad"))
    adata_filtered.X = adata_filtered.layers['log1p_transformed']
    # adata_filtered.var.index = pd.Index(gen.upper() for gen in adata_filtered.var.index.values)
    adata_filtered = adata_filtered[:,adata_filtered.var_names.intersection(adata_integ_clust.var_names)]
    print(sample_id, adata_filtered)
    adata_temp = adata_integ_clust[adata_integ_clust.obs["condition"]==condition,:].copy()
    cols = adata_integ_clust[adata_integ_clust.obs["condition"]==condition,:].obsm['mlm_estimate'].columns
    mlm_estimate_df =  pd.DataFrame(adata_integ_clust[adata_integ_clust.obs["condition"]==condition,:].obsm['mlm_estimate'].values, index=adata_filtered.obs_names, columns=cols)
    adata_filtered.obsm['dorothea_mlm_estimate'] = mlm_estimate_df #adata_temp.obsm['mlm_estimate'].copy()
    adata_dict[sample_id] = adata_filtered

print(mlm_estimate_df)

plot = False
if plot:
    print(len(meta))
    for ind, col in enumerate(adata_filtered.obsm['dorothea_mlm_estimate'].columns):

        sc.pl.umap(acts_integrated, color=[col, 'leiden_0.10'], show=False, cmap=mpl.cm.magma_r, vcenter=0, save=f"_{col}.pdf")
        
        fig, ax = plt.subplots(1, len(meta), figsize=(40,5))
        print(ind, col)
        for ind, row in meta.iterrows():
            sample_id = row["sample_id"]
            condition = row["condition"]
            # mpl.rcParams["image.cmap"]= plt.cm.magma_r
            acts = dc.get_acts(adata_dict[sample_id], obsm_key='dorothea_mlm_estimate')
            sc.pl.spatial(acts, color=col, size=1.25, title=condition, alpha_img=0.5, wspace = 2.5, color_map=mpl.cm.magma_r, hspace=1.5, ax=ax[ind], show=False)
            cbar = ax[ind].collections[0].colorbar
            cbar.set_ticks([])
            cbar = None
            
        plt.savefig(f'{PLOT_PATH}/{col}_tf_activity.pdf');
        # plt.show();
        # plt.clf();


df = dc.rank_sources_groups(acts_integrated, groupby='leiden_0.10', reference='rest', method='t-test_overestim_var')
acts_integrated
n_markers = 5
source_markers = df.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()

dc.plot_network(
    net=net,
    n_sources=["TCF7L2", "SOX11", "ID3", "ZBTB17", "XBP1"],
    n_targets=15,
    node_size=15,
    s_cmap='white',
    t_cmap='white',
    c_pos_w='darkgreen',
    c_neg_w='darkred',
    figsize=(15, 15),
    save=f'{PLOT_PATH}/cluster_8_network.pdf'
)


sc.pl.matrixplot(acts_integrated, source_markers, 'leiden_0.10', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r', save=f"{sample_type}_tf.pdf")
# python vis_tf_act_est.py -i ../data/out_data/visium_integrated_clustered.h5ad  -o ../data/out_data/ -an visium_tf_act_est