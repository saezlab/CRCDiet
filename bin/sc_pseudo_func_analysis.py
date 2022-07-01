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
parser = argparse.ArgumentParser(prog='cluster', description='Cell type annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
# parser.add_argument('-sc', '--sample_column', help='Sample column', required=True)
parser.add_argument('-con', '--condition', help='Condition sample', required=True)
parser.add_argument('-gc', '--group_column', help='Group column', required=True)
parser.add_argument('-ref', '--reference', help='Reference sample', required=True)
parser.add_argument('-lp', '--leiden_param', type=float, help='Final Leiden parameter', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # "sc_pse_func_analys"
condition = args['condition']
group_col = args['group_column']
reference = args['reference']
l_param = args['leiden_param']

# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 

sample_type = "sc"
meta = utils.get_meta_data(sample_type)

adata_integ_clust = sc.read_h5ad(input_path)
# l_param, _ = adata_integ_clust.uns["leiden_best_silh_param"]
if l_param:
    l_param = f"{l_param:.2f}"

# Retrieve PROGENy model weights
progeny = dc.get_progeny(organism='mouse', top=500)

progeny["target"] = progeny["target"].str.capitalize()
# print(progeny)

adata = utils.get_filtered_concat_data(sample_type)

# Store raw counts in layers
adata.layers['counts'] = adata.X

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['normalized'] = adata.X

# print(set(adata.obs["condition"].values))
for ind in range(6):

    adata.var[f'mt-{ind}'] = adata.var[f'mt-{ind}'].astype(str)
    adata.var[f'rp-{ind}'] = adata.var[f'mt-{ind}'].astype(str)

# adata.obs[f"leiden_{l_param}"] = adata_integ_clust.obs[f"leiden_{l_param}"]
adata.obsm["X_umap"] = adata_integ_clust.obsm["X_umap"]
adata.obs[group_col] = adata_integ_clust.obs[group_col] # f"cell_type_{l_param}"

# adata.obs['selection'] = pd.Categorical(adata.obs[f"leiden_{l_param}"]=="1")
# adata.obs['selection'] = adata.obs['selection'].astype(str)

# "{'CD-AOM-DSS-Epi_plus_DN', 'CD-AOM-DSS-Immune', 'HFD-AOM-DSS-Immune', 'LFD-AOM-DSS-Immune', 'LFD-AOM-DSS-Epi_plus_DN', 'HFD-AOM-DSS-Epi_plus_DN'}"

padata = dc.get_pseudobulk(adata, sample_col='condition', groups_col=group_col, layer='counts', min_prop=0.2, min_smpls=3)

logFCs, pvals = dc.get_contrast(adata,
                                group_col=group_col, # f"cell_type_{l_param}"
                                condition_col='condition',
                                # condition='HFD-AOM-DSS-Epi_plus_DN',
                                # reference='LFD-AOM-DSS-Epi_plus_DN',
                                condition= condition, # 'HFD-AOM-DSS-Immune',
                                reference= reference,  #'LFD-AOM-DSS-Immune',
                                method='t-test'
                               )

deg = dc.format_contrast_results(logFCs, pvals)
# print(deg)
"""
How the activities changes in "condition" with respect to "reference"
"""

# Infer pathway activities with mlm
pathway_acts, pathway_pvals = dc.run_mlm(mat=deepcopy(logFCs).astype('float64'), net=progeny, source='source', target='target', weight='weight')

sns.clustermap(pathway_acts, center=0, cmap='coolwarm')
plt.savefig(f"{PLOT_PATH}/comparative_pathway_act_est_{condition}_wrt_{reference}_{group_col}.pdf")
plt.show();

# dorothea = dc.get_dorothea()
"""
It looks like JAK-STAT is active across cell types in COVID-19 compared to healthy. 
To further explore how the target genes behave, we can plot them in a volcano plot:
"""

"""mpl.rcParams["font.size"] = 7
logFCs = logFCs[(logFCs>-10) & (logFCs<10)]
for c_t in logFCs.index:

    dc.plot_volcano(deepcopy(logFCs).astype('float64'), deepcopy(pvals).astype('float64'), 'TGFb', c_t, progeny, top=10, source='source', target='target', weight='weight', sign_thr=0.05, lFCs_thr=0.5)
    plt.show()"""


"""sc.pl.umap(adata[adata[: , 'Gata3'].X > 7.5, :], size=2.0 )
sc.pl.umap(adata[adata[: , 'Foxp3'].X > 4.0, :], size=2.0 )
sc.pl.umap(adata, size=1.5 )"""

## Dorothea Part

"""# Retrieve DoRothEA gene regulatory network
dorothea = dc.get_dorothea()
dorothea["target"] = dorothea["target"].str.capitalize()

# Infer pathway activities with mlm
tf_acts, tf_pvals = dc.run_mlm(mat=deepcopy(logFCs).astype('float64'), net=dorothea, source='source', target='target', weight='weight')
tf_acts

# Get top 5 most active/inactive sources
top = 5
top_idxs = set()
for row in tf_acts.values:
    sort_idxs = np.argsort(-np.abs(row))
    top_idxs.update(sort_idxs[:top])
top_idxs = np.array(list(top_idxs))
top_names = tf_acts.columns[top_idxs]
names = tf_acts.index.values
top = pd.DataFrame(tf_acts.values[:,top_idxs], columns=top_names, index=names)

# Plot
sns.clustermap(top, center=0, cmap='coolwarm')
plt.show()
"""


# Comparative Pseudobulk Pathway Analysis - Major Cell Types <a class="anchor" id="seventh-bullet"></a>

## Immune cells
### <i>HFD-AOM-DSS-Immune </i> compared <i>LFD-AOM-DSS-Immune </i><a class="anchor" id="seventh-bullet"></a>
# %run sc_pseudo_func_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an 'pseudobulk_functional_annot' -gc 'major_cell_types' -con 'HFD-AOM-DSS-Immune' -ref 'LFD-AOM-DSS-Immune'
### <i>HFD-AOM-DSS-Immune </i> compared <i>CD-AOM-DSS-Immune </i><a class="anchor" id="seventh-bullet"></a>
# %run sc_pseudo_func_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an 'pseudobulk_functional_annot' -gc 'major_cell_types' -con 'HFD-AOM-DSS-Immune' -ref 'CD-AOM-DSS-Immune'
### <i>LFD-AOM-DSS-Immune </i> compared <i>CD-AOM-DSS-Immune </i><a class="anchor" id="seventh-bullet"></a>
# %run sc_pseudo_func_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an 'pseudobulk_functional_annot' -gc 'major_cell_types' -con 'LFD-AOM-DSS-Immune' -ref 'CD-AOM-DSS-Immune'

## Epithelial cells
### <i>HFD-AOM-DSS-Epi_plus_DN </i> compared <i>LFD-AOM-DSS-Epi_plus_DN </i><a class="anchor" id="seventh-bullet"></a>
# %run sc_pseudo_func_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an 'pseudobulk_functional_annot' -gc 'major_cell_types' -con 'HFD-AOM-DSS-Epi_plus_DN' -ref 'LFD-AOM-DSS-Epi_plus_DN'
### <i>HFD-AOM-DSS-Epi_plus_DN </i> compared <i>CD-AOM-DSS-Epi_plus_DN </i><a class="anchor" id="seventh-bullet"></a>
# %run sc_pseudo_func_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an 'pseudobulk_functional_annot' -gc 'major_cell_types' -con 'HFD-AOM-DSS-Epi_plus_DN' -ref 'CD-AOM-DSS-Epi_plus_DN'
### <i>LFD-AOM-DSS-Epi_plus_DN </i> compared <i>CD-AOM-DSS-Epi_plus_DN </i><a class="anchor" id="seventh-bullet"></a>
# %run sc_pseudo_func_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an 'pseudobulk_functional_annot' -gc 'major_cell_types' -con 'LFD-AOM-DSS-Epi_plus_DN' -ref 'CD-AOM-DSS-Epi_plus_DN'



# python sc_pseudo_func_analysis.py -i ../data/out_data/_integrated_cluster_scannot.h5ad -o ../data/out_data

# jupyter nbconvert sc_00_pipeline.ipynb --to html_embed --template toc2 --output ../docs/crcdiet_pipelines/crcdiet_sc_pipeline.html 