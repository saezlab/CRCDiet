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
parser = argparse.ArgumentParser(prog='Dorothea', description='TF activity estimation')
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

sample_type = "sc"
meta = utils.get_meta_data(sample_type)

adata_merged = sc.read_h5ad(input_path)
adata_integ_clust = sc.read_h5ad(os.path.join(output_path, f'{sample_type}_integrated_cluster_scannot.h5ad'))
adata_merged.X = adata_merged.layers['log1p_transformed']

# Retrieve PROGENy model weights
net = dc.get_dorothea(organism='mouse', levels=['A','B','C'])
net["target"] = net["target"].str.upper()
adata_integ_clust.X = adata_integ_clust.layers['log1p_transformed']
adata_integ_clust.var.index = pd.Index(gen.upper() for gen in adata_integ_clust.var.index.values)
dc.run_mlm(mat=adata_integ_clust, net=net, source='source', target='target', weight='weight', verbose=True, use_raw=False)
adata_integ_clust.obsm['dorothea_mlm_estimate'] = adata_integ_clust.obsm['mlm_estimate'].copy()
adata_integ_clust.obsm['dorothea_mlm_pvals'] = adata_integ_clust.obsm['mlm_pvals'].copy()
acts_integrated = dc.get_acts(adata_integ_clust, obsm_key='dorothea_mlm_estimate')
# print(adata_integ_clust.obsm['dorothea_mlm_estimate'])

for ind, col in enumerate(adata_integ_clust.obsm['dorothea_mlm_estimate'].columns):

    sc.pl.umap(acts_integrated, color=col, show=False, cmap='coolwarm', vcenter=0, save=f"{col}.pdf")

mean_acts = dc.summarize_acts(acts_integrated, groupby='cell_type_0.20', min_std=0.75)

sns.clustermap(mean_acts, xticklabels=mean_acts.columns, vmin=-5, vmax=5, cmap='coolwarm')
plt.savefig(f"{PLOT_PATH}/{sample_type}_cell_type_pathway_activity_est_cmap.pdf")
#plt.show();

mean_acts = dc.summarize_acts(acts_integrated, groupby='condition', min_std=0.75)

sns.clustermap(mean_acts, xticklabels=mean_acts.columns, vmin=-5, vmax=5, cmap='coolwarm')
plt.savefig(f"{PLOT_PATH}/{sample_type}_condition_pathway_activity_est_cmap.pdf")
#plt.show();
# , save=f'{sample_type}_pathway_activity_est_clustermap'


# Write to file
# adata_integ_clust.write(os.path.join(output_path, f'{sample_type}_integrated_clustered.h5ad'))

# python sc_tf_act_est.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data/ -an sc_tf_act_est




