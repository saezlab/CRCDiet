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


sc.pl.umap(adata, color=f'leiden_{res_param:.2f}', palette=sc.pl.palettes.default_20, size=8 , show=False, legend_loc='on data', save = f'{sample_type}_subcluster_final')
print("RANK")
sc.tl.rank_genes_groups(adata, groupby=f"leiden_{res_param:.2f}", method='wilcoxon', key_added = f"wilcoxon_{res_param:.2f}")



mpl.rcParams['axes.titlesize'] = 20
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key=f"wilcoxon_{res_param:.2f}", show=False, groupby=f"leiden_{res_param:.2f}", save=f'{sample_type}_one_vs_rest_subcluster_{res_param:.2f}')#
mpl.rcParams['axes.titlesize'] = 60
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=f"wilcoxon_{res_param:.2f}", show=False, groupby=f"leiden_{res_param:.2f}", save=f'{sample_type}_deg_clusters_dotplot_subcluster_{res_param:.2f}')

print(f"Saving the object... {sample_type}_integrated_subclustered.h5ad...")
# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_integrated_subclustered.h5ad'))

# python sc_subcluster.py -i ../data/out_data/atlas_integrated_clustered.h5ad -o ../data/out_data -st atlas  -an atlas_cluster