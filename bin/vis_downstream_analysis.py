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

adata_merged = sc.read_h5ad(input_path)
adata_integ_clust = sc.read_h5ad(os.path.join(output_path, f'{sample_type}_integrated_clustered.h5ad'))
adata_merged.X = adata_merged.layers['log1p_transformed']

adata_merged.var.index = pd.Index(gen.upper() for gen in adata_merged.var.index.values)
adata_merged.obsm["X_umap"] = adata_integ_clust.obsm["X_umap"]

# print(adata_integ_clust)

# print(meta)

adata_merged = adata_merged[adata_merged.obs["condition"].isin(["CD-no-AOM-DSS", "HFD-no-AOM-DSS", "LFD-no-AOM-DSS"]),:]
adata_integ_clust = adata_integ_clust[adata_integ_clust.obs["condition"].isin(["CD-no-AOM-DSS", "HFD-no-AOM-DSS", "LFD-no-AOM-DSS"]),:]
# print(adata_merged)





# Retrieve PROGENy model weights
progeny = dc.get_progeny(organism='mouse', top=500)

progeny["target"] = progeny["target"].str.upper()

dc.run_mlm(mat=adata_merged, net=progeny, source='source', target='target', weight='weight', verbose=True, use_raw=False)

adata_integ_clust.obsm['progeny_mlm_estimate'] = adata_merged.obsm['mlm_estimate'].copy()
adata_integ_clust.obsm['progeny_mlm_pvals'] = adata_merged.obsm['mlm_pvals'].copy()

acts = dc.get_acts(adata_integ_clust, obsm_key='progeny_mlm_estimate')

sc.pl.umap(acts, color=adata_integ_clust.obsm['progeny_mlm_estimate'].columns, vcenter=0, cmap='coolwarm', show=True) # , save=f'{sample_type}_pathway_activity_est')

mean_acts = dc.summarize_acts(acts, groupby='condition', min_std=0)
# adata_samp_cc = adata_cc_merged[adata_cc_merged.obs["condition"]==condition,:]
sns.clustermap(mean_acts, xticklabels=mean_acts.columns, vmin=-2, vmax=2, cmap='coolwarm')
# plt.savefig(f"{PLOT_PATH}/{sample_type}_cell_type_pathway_activity_est_cmap.pdf")
plt.show();
# , save=f'{sample_type}_pathway_activity_est_clustermap'


