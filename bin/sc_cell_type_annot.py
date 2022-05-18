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

sc.settings.verbosity = 0

warnings.simplefilter(action='ignore')

# Read integrated object
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Cell type annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']

###############################

sample_type = "sc"
meta = utils.get_meta_data("sc")

S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
PLOT_PATH =  os.path.join(S_PATH, "../plots", "sc_cell_type_annot")

Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = PLOT_PATH
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80)

adata_integ_clust = sc.read_h5ad(input_path)
l_param, _ = adata_integ_clust.uns["leiden_best_silh_param"]
l_param = f"{l_param:.2f}"

adata = utils.get_filtered_concat_data(sample_type)
adata.obs[f"leiden_{l_param}"] = adata_integ_clust.obs[f"leiden_{l_param}"]

marker_genes = None
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

dc.run_ora(mat=adata, net=marker_genes, use_raw=False, source='cell_type', target='genesymbol', min_n=3, verbose=True)
acts = dc.get_acts(adata, obsm_key='ora_estimate')

# print("============= ORA p-Vals ============= ")
# print(adata.obsm["ora_pvals"])
# print("============= ora_estimate ============= ")
# print(adata.obsm["ora_estimate"])
for c_t in list(adata.obsm["ora_estimate"].columns):
    adata.obs[c_t] = adata.obsm["ora_estimate"][c_t]


adata.obsm["X_umap"] = adata_integ_clust.obsm["X_umap"]

mpl.rcParams['axes.titlesize'] = 20
cell_types = list(set(marker_genes["cell_type"]))
sc.pl.umap(adata, color=list(adata.obsm["ora_pvals"].columns))

dict_mean_enr = dict()

mean_enr = dc.summarize_acts(acts, groupby=f'leiden_{l_param}')

annotation_dict = dc.assign_groups(mean_enr)

adata.obs[f'cell_type_{l_param}'] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{l_param}']]

# sc.pl.umap(adata, color=f'cell_type_{l_param}', show=True, legend_loc="on data", legend_fontsize="xx-small")

# python sc_cell_type_annot.py -i ../data/out_data/sc_integrated_clustered.h5ad -o ../data/out_data

