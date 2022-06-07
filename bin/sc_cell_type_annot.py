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
parser = argparse.ArgumentParser(prog='cluster', description='Cell type annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths("sc_cell_type_annot")
############################### BOOOORIING STUFF ABOVE ############################### 

sample_type = "sc"
meta = utils.get_meta_data("sc")

adata_integ_clust = sc.read_h5ad(input_path)
l_param, _ = adata_integ_clust.uns["leiden_best_silh_param"]
l_param = f"{l_param:.2f}"

adata = utils.get_filtered_concat_data(sample_type)

for ind in range(6):
    adata.var[f'mt-{ind}'] = adata.var[f'mt-{ind}'].astype(str)
    adata.var[f'rp-{ind}'] = adata.var[f'mt-{ind}'].astype(str)

adata.obs[f"leiden_{l_param}"] = adata_integ_clust.obs[f"leiden_{l_param}"]

marker_genes = None
own_cell_type_list = None
own_markers = True
both = True
marker_db="PanglaoDB"

if own_markers:
    # adata = adata.raw.to_adata()
    # TODO: make this parametric
    marker_genes = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
    marker_genes["genesymbol"] = marker_genes["genesymbol"].str.capitalize()
    own_cell_type_list = list(set(marker_genes["cell_type"]))
    
    if both:
        marker_genes_pk = dc.get_resource('PanglaoDB')
        # marker_db="PanglaoDB"
        # marker_genes = pd.read_csv(os.path.join(DATA_PATH,"PanglaoDB_markers_27_Mar_2020.tsv"), sep="\t")
        # print(marker_genes)
        # Filter by canonical_marker and human
        # Filter by canonical_marker and human
        marker_genes_pk = marker_genes_pk[(marker_genes_pk['mouse']=='True')&(marker_genes_pk['canonical_marker']=='True')]

        # Remove duplicated entries
        # marker_genes = marker_genes[marker_genes["species_official"].str.contains("mouse") & marker_genes["canonical_marker"]==1.0]
        marker_genes_pk = marker_genes_pk[['cell_type', 'genesymbol']]
        # print(marker_genes_pk)
        marker_genes_pk = marker_genes_pk[marker_genes_pk.cell_type != "Müller cells"]
        # print(marker_genes_pk)
        # Remove duplicated entries
        marker_genes_pk = marker_genes_pk[~marker_genes_pk.duplicated(['cell_type', 'genesymbol'])]
        marker_genes = pd.concat([marker_genes, marker_genes_pk])
        # print(marker_genes)
    
    
    marker_genes['weight'] = 1
    marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]


else:
    
    # print(marker_genes)
    marker_genes = dc.get_resource('PanglaoDB')
    # marker_db="PanglaoDB"
    # marker_genes = pd.read_csv(os.path.join(DATA_PATH,"PanglaoDB_markers_27_Mar_2020.tsv"), sep="\t")
    # print(marker_genes)
    # Filter by canonical_marker and human
    # Filter by canonical_marker and human
    marker_genes = marker_genes[(marker_genes['mouse']=='True')&(marker_genes['canonical_marker']=='True')]

    # Remove duplicated entries
    # marker_genes = marker_genes[marker_genes["species_official"].str.contains("mouse") & marker_genes["canonical_marker"]==1.0]
    
    # print(marker_genes["genesymbol"])
    # Remove duplicated entries
    marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]


marker_genes["genesymbol"] = marker_genes["genesymbol"].str.capitalize()
marker_genes["cell_type"] = marker_genes["cell_type"].str.capitalize()
marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]

dc.run_ora(mat=adata, net=marker_genes, use_raw=False, source='cell_type', target='genesymbol', min_n=3, verbose=True)
acts = dc.get_acts(adata, obsm_key='ora_estimate')
# print(acts)

# print("============= ORA p-Vals ============= ")
# print(adata.obsm["ora_pvals"])
# print("============= ora_estimate ============= ")
# print(adata.obsm["ora_estimate"])
for c_t in list(adata.obsm["ora_estimate"].columns):
    adata.obs[c_t] = adata.obsm["ora_estimate"][c_t]

# print(list(adata.obsm["ora_estimate"].columns))
adata.obsm["X_umap"] = adata_integ_clust.obsm["X_umap"]

mpl.rcParams['axes.titlesize'] = 10
cell_types = list(set(marker_genes["cell_type"]))
# sc.pl.umap(adata, color=set(adata.obsm["ora_estimate"].columns))#&set(own_cell_type_list))
sc.pl.umap(adata, color=set(adata.obsm["ora_estimate"].columns) &set(own_cell_type_list))

print("Regeneration plots per sample...")
for ind, row in meta.iterrows():
    sample_id = row["sample_id"]
    condition = row["condition"]
    adata_tmp = adata[adata.obs["condition"]==condition,:]
    #adata = sc.read_h5ad(os.path.join(input_path,f"{sample_id}_filtered.h5ad"))
    #adata_samp_cc = adata_cc_merged[adata_cc_merged.obs["condition"]==condition,:]
    sc.pl.umap(adata_tmp, color="Regeneration", title=f"{condition} - Regeneration", save=f"Regeneration - {condition}")



dict_mean_enr = dict()

mean_enr = dc.summarize_acts(acts, groupby=f'leiden_{l_param}')
# print(mean_enr)

annotation_dict = dc.assign_groups(mean_enr)
# print(annotation_dict)
# Manual annotation

annotation_dict['4'] = "Myeloid cells"
annotation_dict['20'] = "Enteroendocrine"
annotation_dict['10'] = "ILC2"
annotation_dict['12'] = "Goblet cells"
annotation_dict['1'] = "Stroma"
annotation_dict['5'] = "Stroma"
annotation_dict['8'] = "Stroma"
# annotation_dict['9'] = "Regeneration"
annotation_dict['9'] = "Myofibroblasts"
annotation_dict['15'] = "Stroma"
annotation_dict['17'] = "Keratynocytes"
annotation_dict['23'] = "Prolif."
annotation_dict['3'] = "Prolif. + Mature enterocytes"

# main cell type annot dict
main_ct_annot_dict = dict()
immnune_clusters = [0, 2, 4, 6, 7, 10, 11, 13, 18, 19, 22, 24]
epithelial_clusters = [17, 23, 3, 21, 20, 12]
stroma_clusters =  [16, 15, 9, 1, 5, 8, 14]
for ind in range(25):
    if ind in immnune_clusters:
        main_ct_annot_dict[str(ind)] = "Immune cells"
    elif ind in epithelial_clusters:
        main_ct_annot_dict[str(ind)] = "Epithelial cells"
    else:
        main_ct_annot_dict[str(ind)] = "Stroma cells"

adata.obs[f'major_cell_types'] = [main_ct_annot_dict[clust] for clust in adata.obs[f'leiden_{l_param}']]

"""{'0': 'B cells', '1': 'Fibroblasts', '10': 'T cells', '11': 'Dendritic cells', '12': 'Goblet cells', 
'13': 'Neutrophils', '14': 'Endothelial cells', '15': 'Fibroblasts', '16': 'Endothelial cells', '17': 'Fibroblasts', 
'18': 'Mast cells', '19': 'Plasma cells', '2': 'T cells', '20': 'Fibroblasts', '21': 'Tuft cells', 
'22': 'Plasma cells', '23': 'Fibroblasts', '24': 'T cells', '3': 'Goblet cells', 
'4': 'Macrophages', '5': 'Fibroblasts', '6': 'Plasma cells', '7': 'T cells', '8': 'Fibroblasts', '9': 'Fibroblasts'}"""
adata.obs[f'cell_type_{l_param}'] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{l_param}']]


mpl.rcParams['figure.dpi']= 300
mpl.rcParams["figure.figsize"] = (10,10)
mpl.rcParams['axes.facecolor'] = "white"

# the number of genes expressed in the count matrix
sc.pl.umap(adata, color=f'cell_type_{l_param}', title= ["Cell type annotation"], show=True, s=10, legend_loc="on data", legend_fontsize="xx-small",  save=f'{sample_type}_cell_type_annot_{l_param}')
sc.pl.umap(adata, color='major_cell_types', title= ["Major cell types"], show=True, s=10, legend_loc="on data", legend_fontsize="xx-small",  save=f'{sample_type}_major_cell_types')
adata_integ_clust.obs[f'cell_type_{l_param}'] = adata.obs[f'cell_type_{l_param}']
adata_integ_clust.obs['major_cell_types'] = adata.obs['major_cell_types']
adata_integ_clust.write(os.path.join(output_path, f'{sample_type}_integrated_cluster_scannot.h5ad'))
# python sc_cell_type_annot.py -i ../data/out_data/sc_integrated_clustered.h5ad -o ../data/out_data

