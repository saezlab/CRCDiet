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
parser.add_argument('-gb', '--group_by', help='Group by cell type, condition etc. for comparison', required=False)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)
parser.add_argument('-samp', '--samples', help='Samples to subset', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # "sc_pse_func_analys"  
group_by = args['group_by']
sample_type = args['sample_type']
samples = args['samples']

# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 

meta = utils.get_meta_data(sample_type)

# adata_merged = sc.read_h5ad(input_path)

adata_merged = utils.get_filtered_concat_data(sample_type)
if "log1p_transformed" not in adata_merged.layers:
    adata_merged.layers["log1p_transformed"] = sc.pp.normalize_total(adata_merged, inplace=False, target_sum=1e6)["X"]
    sc.pp.log1p(adata_merged, layer="log1p_transformed")

# print(adata_merged.obs_names)
adata_integ_clust = sc.read_h5ad(os.path.join(output_path, f'{sample_type}_integrated_clustered.h5ad'))
adata_merged.X = adata_merged.layers['log1p_transformed']
print(adata_merged.var_names)

# run only on HVGs
adata_merged = adata_merged[adata_integ_clust.obs_names,adata_integ_clust.var_names]

# print(adata_integ_clust)

adata_merged.var.index = pd.Index(gen.upper() for gen in adata_merged.var.index.values)
adata_merged.obsm["X_umap"] = adata_integ_clust.obsm["X_umap"]

# Retrieve PROGENy model weights
progeny = dc.get_progeny(organism='mouse', top=500)

progeny["target"] = progeny["target"].str.upper()

dc.run_mlm(mat=adata_merged, net=progeny, source='source', target='target', weight='weight', verbose=True, use_raw=False)

# adata_integ_clust.obsm['progeny_mlm_estimate'] = adata_merged.obsm['mlm_estimate'].copy()
# adata_integ_clust.obsm['progeny_mlm_pvals'] = adata_merged.obsm['mlm_pvals'].copy()

adata_merged.obsm['progeny_mlm_estimate'] = adata_merged.obsm['mlm_estimate'].copy()
adata_merged.obsm['progeny_mlm_pvals'] = adata_merged.obsm['mlm_pvals'].copy()


adata_merged.obsm['progeny_mlm_estimate'].to_csv(f"{OUT_DATA_PATH}/{analysis_name}_mlm_estimate.csv")

if samples:

    lst_samples = samples.split(",")
    adata_merged = adata_merged[adata_merged.obs["condition"].isin(lst_samples)]
    

acts = dc.get_acts(adata_merged, obsm_key='mlm_estimate')

sc.pl.umap(acts, color=adata_merged.obsm['mlm_estimate'].columns, vcenter=0, cmap='coolwarm', save=f'{sample_type}_pathway_activity_est')


mean_acts = dc.summarize_acts(acts, groupby=group_by, min_std=0)

sns.clustermap(mean_acts, xticklabels=mean_acts.columns, vmin=-2, vmax=2, cmap='coolwarm')
plt.savefig(f"{PLOT_PATH}/{sample_type}_{group_by}_pathway_activity_est_cmap.pdf")
plt.show();


# Write to file
# adata_integ_clust.write(os.path.join(output_path, f'{sample_type}_integrated_clustered.h5ad'))

# python sc_pathway_act_est.py -i ../data/out_data/sc_b_cells_integrated.h5ad -o ../data/out_data/ -an sc_bcells_pathway_act_est -gb condition -samp "CD-AOM-DSS-Immune,HFD-AOM-DSS-Immune,LFD-AOM-DSS-Immune"
# python sc_pathway_act_est.py -i ../data/out_data/sc_epicells_integrated_clustered.h5ad -o ../data/out_data/ -an sc_epi_cells_aom_noaom_pathway_act_est -gb condition -st sc_epicells
# python sc_pathway_act_est.py -i ../data/out_data/sc_epicells_only_integrated_clustered.h5ad -o ../data/out_data/ -an sc_epi_cells_aom_noaom_pathway_act_est -gb condition -st sc_epicells



