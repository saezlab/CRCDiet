from copyreg import pickle
from genericpath import sameopenfile
from operator import index
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns 
# import decoupler as dc  
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import os
from sklearn.metrics import silhouette_score, pairwise_distances
import pickle
import warnings
import utils
import matplotlib as mpl

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=50, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Run Clustering and annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
parser.add_argument('-of', '--output_file', help='Output file name', required=False)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # sc_cluster
output_file = args['output_file']
sample_type = args['sample_type']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 

plt.rcParams['figure.dpi']= 300
plt.rcParams['figure.figsize']= (45, 30)


res_param = 0.4
new_res_param = 0.1
adata = sc.read_h5ad(input_path)
print(adata)
adata.uns['log1p']["base"] = None

to_be_subcluster = [6, 7, 13, 28, 30, 36, 38]
for cluster_number in to_be_subcluster:
    sc.tl.leiden(adata, restrict_to=(f'leiden_{res_param:.2f}', [str(cluster_number)]),  resolution=new_res_param, key_added=f'leiden_{res_param:.2f}')

new_res_param = 0.2
to_be_subcluster = [17]
for cluster_number in to_be_subcluster:
    sc.tl.leiden(adata, restrict_to=(f'leiden_{res_param:.2f}', [str(cluster_number)]),  resolution=new_res_param, key_added=f'leiden_{res_param:.2f}')


annotation_dict = {"0" : "Epithelial cells (colon)", "1" : "Epithelial cells (SI enterocytes)", "2" : "Epithelial cells (SI enterocytes)", "3" : "IECs (goblet cells)", "4" : "Epithelial cells (colon)", "5" : "Stroma/Fibroblasts", "6,0" : "Cells with high mitochondiral activity", "6,1" : "??", "6,2" : "Prolif, + Mature enterocytes ", "6,3" : "??", "6,4" : "??", "7,0" : "Immune/Endotheial cells", "7,1" : "T cells", "7,2" : "T cells/NK cells", "7,3" : "Monocyte/Macrophage", "7,4" : "T cells ?", "7,5" : "??", "7,6" : "T cells/NK cells", "8" : "Epithelial cells (SI & colon)", "9" : "IECs (goblet cells)", "10" : "IECs (goblet cells)", "11" : "Epithelial cells?", "12" : "B cells", "13" : "", "13,0" : "IECs (goblet cells)", "13,1" : "IECs (goblet cells)", "13,2" : "IECs (goblet cells)", "13,3" : "IECs (goblet cells)", "14" : "Neuronal cells (enteric neurons?)", "15" : "T cells", "16" : "Neuronal cells (enteric neurons?)", "17,0" : "Macrophages/ APCs?", "17,1" : "Macrophages/ APCs? B cells?", "17,2" : "Macrophages/ APCs?", "17,3" : "Macrophages/ APCs?", "17,4" : "Immune cells", "17,5" : "B cells?", "17,6" : "T cells/B cells", "17,7" : "B cells", "17,8" : "Immune cells capable of producing defensins", "17,9" : "Macophages/DCs?", "18" : "Epithelial cells (SI & colon)", "19" : "Endothelial cells", "20" : "Neuronal cells (enteric neurons?)", "21" : "Immune cells", "22" : "Plasma cells", "23" : "Immune cells", "24" : "Epithelial cells (SI enterocytes)", "25" : "Paneth cells", "26" : "B cells", "27" : "Neuronal cells (enteric neurons?)", "28,0" : "Endothelial cells", "28,1" : "Endothelial cells ", "28,2" : "Endothelial cells ", "28,3" : "Endothelial cells ", "28,4" : "??", "29" : "T cells", "30,0" : "Tuft cells", "30,1" : "IECs/Tuft cells", "30,2" : "Tuft cells", "30,3" : "IECs/Tuft cells", "30,4" : "Tuft cells", "30,5" : "IECs/Tuft cells", "31" : "Goblet cells", "32" : "Immune cells?", "33" : "Enteroendocrine cells", "34" : "Neuronal cells (enteric neurons?)", "35" : "Smooth muscle cells", "36" : "", "36,0" : "Mast cells (Immune cells : granulocutes)", "36,1" : "Neutrophils (Immune cells : granulocutes)", "36,2" : "Immune cells : granulocutes", "36,3" : "Mast cells ", "36,4" : "Neutrophils?", "36,5" : "Neutrophils?", "36,6" : "Immune cells", "37" : "Epithelial cells (SI & colon)", "38,0" : "Epithelial cells (colon)", "38,1" : "B cells? contamination?", "38,2" : "??", "38,3" : "??", "39" : "Paneth/Goblet cell", "40" : "Goblet cells", "41" : "Paneth cells", "42" : "Fibroblasts", "43" : "Fibroblasts", "44" : "Plasma cells", "45": "??", "46" : "Keratinocytes (distal colonic transitional epithelium?)", "47" : "Germinal center associated immune cells", "48" : "B cells (Plasma cells)", "49" : "??", "50" : "Paneth/Goblet cells"}
adata.obs["cell_type"] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{res_param:.2f}']]
print(set(adata.obs[f"leiden_{res_param:.2f}"]))
# sc.pl.umap(adata, color=f'leiden_{res_param:.2f}', palette=sc.pl.palettes.default_20, size=8 , show=False, legend_loc='on data', save = f'{sample_type}_subcluster_final')
sc.pl.umap(adata, color="cell_type", palette=sc.pl.palettes.default_20, size=8 , show=False, legend_loc='on data', save = f'{sample_type}_cell_type_annotation')

"""
sc.tl.rank_genes_groups(adata, groupby=f"leiden_{res_param:.2f}", method='wilcoxon', key_added = f"wilcoxon_{res_param:.2f}")



mpl.rcParams['axes.titlesize'] = 20
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key=f"wilcoxon_{res_param:.2f}", show=False, groupby=f"leiden_{res_param:.2f}", save=f'{sample_type}_one_vs_rest_subcluster_{res_param:.2f}')#
mpl.rcParams['axes.titlesize'] = 60
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=f"wilcoxon_{res_param:.2f}", show=False, groupby=f"leiden_{res_param:.2f}", save=f'{sample_type}_deg_clusters_dotplot_subcluster_{res_param:.2f}')
"""
print(f"Saving the object... {sample_type}_integrated_subclustered.h5ad...")
# Write to file
# adata.write(os.path.join(output_path, f'{sample_type}_integrated_subclustered.h5ad'))

# python sc_subcluster.py -i sc-o ../data/out_data -st atlas  -an atlas_cluster