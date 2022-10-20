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

compute_silh = False
res_param = 0.4
meta = utils.get_meta_data(sample_type)

if output_file:
    sample_type = f"{sample_type}_{output_file}"

adata = sc.read_h5ad(input_path)

dist_mat = None
adata.uns['log1p']["base"] = None

if compute_silh:
    if not os.path.exists(f'{OUT_DATA_PATH}/{sample_type}_dist_mat.pickle'):
        print("Calculating distance matrix... ")
        dist_mat = pairwise_distances(adata.obsm["X_pca"], metric='euclidean') #  , n_jobs=-1)

        with open(f'{OUT_DATA_PATH}/{sample_type}_dist_mat.pickle', 'wb') as handle:
            pickle.dump(dist_mat, handle, protocol=pickle.HIGHEST_PROTOCOL)

    else:
        print("Loading precomputed distance matrix... ")
        with open(f'{OUT_DATA_PATH}/{sample_type}_dist_mat.pickle', 'rb') as handle:
            dist_mat = pickle.load(handle)

    best_res_param, best_s_scr = -1, -1
    step = 0.10

    silh_param_scores = []
    # perform clustering, Rank genes for characterizing groups, plot top 5 genes
    for l_res_param in np.arange(0.1, 1.01, step):
    # for l_param in [0.30]:

        print(f"Creating clusters with Leiden resolution param: {l_res_param:.2f}")
        sc.tl.leiden(adata, resolution = l_res_param, key_added = f"leiden_{l_res_param:.2f}") # default resolution in 1.0
        silh_scr = silhouette_score(dist_mat, np.array(adata.obs[f"leiden_{l_res_param:.2f}"]), metric='precomputed')
        print(f"Clustering param: {l_res_param:.2f}\tSilhoutte score: {silh_scr:.3f}")
        silh_param_scores.append((l_res_param, silh_scr))

    for l_res_param, s_scr in silh_param_scores:
        if s_scr > best_s_scr:
            best_res_param, best_s_scr = l_res_param, s_scr

    adata.uns["leiden_best_silh_param"] = [best_res_param, best_s_scr]

    l_param, _ = adata.uns["leiden_best_silh_param"]
    l_param = f"{l_param:.2f}"

else:
    print("Computing neighbourhood graph... ")
    sc.pp.neighbors(adata)
    print("Leiden")
    sc.tl.leiden(adata, resolution=res_param, key_added = f"leiden_{res_param:.2f}")
    adata.uns["leiden_best_silh_param"] = res_param

sc.tl.rank_genes_groups(adata, groupby=f"leiden_{res_param:.2f}", method='wilcoxon', key_added = f"wilcoxon_{res_param}")

print(f"Saving the object... {sample_type}_integrated_clustered.h5ad...")
# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_integrated_clustered.h5ad'))


# python sc_cluster.py -i ../data/out_data/atlas_integrated.h5ad -o ../data/out_data -st atlas  -an atlas_cluster
# python sc_cluster_annotate.py -i ../data/out_data/atlas_integrated_clustered.h5ad -o ../data/out_data -an atlas_cluster -st atlas
