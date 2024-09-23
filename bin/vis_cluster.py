from copyreg import pickle
from genericpath import sameopenfile
from operator import index
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns 
# import decoupler as dc  
import scanpy as sc
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np
import matplotlib as mpl
import argparse
import os
from sklearn.metrics import silhouette_score, pairwise_distances
import pickle
import utils
import math
import warnings
from random import shuffle



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

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # sc_cluster
output_file = args['output_file']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 

print("Clustering the visium slides...")
sample_type = "visium"

adata = sc.read_h5ad(input_path)

markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
markers = list(set(markers_df["genesymbol"].str.capitalize()))

meta = utils.get_meta_data(sample_type)

step = 0.10
dist_mat = None

"""if not os.path.exists(f'{OUT_DATA_PATH}/{sample_type}_dist_mat.pickle') or True:
    print("Calculating distance matrix... ")
    dist_mat = pairwise_distances(adata.obsm["X_pca"], metric='euclidean') #  , n_jobs=-1)

    with open(f'{OUT_DATA_PATH}/{sample_type}_dist_mat.pickle', 'wb') as handle:
        pickle.dump(dist_mat, handle, protocol=pickle.HIGHEST_PROTOCOL)

else:
    print("Loading precomputed distance matrix... ")
    with open(f'{OUT_DATA_PATH}/{sample_type}_dist_mat.pickle', 'rb') as handle:
        dist_mat = pickle.load(handle)"""

print("Computing neighbourhood graph... ")
sc.pp.neighbors(adata)

# mcolors.CSS4_COLORS["lightskyblue"]
# #D62728
# #279E68
palette_color_lst = ["#1F77B4", "#FFFF00", '#D62728', '#279E68', mcolors.CSS4_COLORS["darkgreen"], 'tab:purple', '#8C564B', mcolors.CSS4_COLORS["rosybrown"], 'tab:gray', '#E377C2', 'tab:cyan']
best_l_param, best_s_scr = -1, -1

silh_param_scores = []
# perform clustering, Rank genes for characterizing groups, plot top 5 genes
for l_param in np.arange(0.15, 1.01, step):
    print(f"Creating clusters with Leiden resolution param: {l_param:.2f}")
    new_res_param=0.05
    
    
    sc.tl.leiden(adata, resolution = l_param, key_added = f"leiden_{l_param:.2f}") # default resolution in 1.0
    adata.obs[f'leiden_{l_param:.2f}'][adata.obs[f'leiden_{l_param:.2f}'].isin(['0', '2'])]='0'
    sc.tl.leiden(adata, restrict_to=(f'leiden_{l_param:.2f}', ["0"]),  resolution=0.1, key_added=f'leiden_{l_param:.2f}')
    sc.tl.leiden(adata, restrict_to=(f'leiden_{l_param:.2f}', ["0,2"]),  resolution=0.1, key_added=f'leiden_{l_param:.2f}')
    sc.tl.leiden(adata, restrict_to=(f'leiden_{l_param:.2f}', ["0,2,0"]),  resolution=0.5, key_added=f'leiden_{l_param:.2f}')
    adata.obs[f'leiden_{l_param:.2f}'][adata.obs[f'leiden_{l_param:.2f}'].isin(["0,2,0,0", "0,2,0,1"])]="0,2,0,0"
    adata.obs[f'leiden_{l_param:.2f}'] = adata.obs[f'leiden_{l_param:.2f}'].cat.remove_categories(["0,2,0,1"])
    adata.obs[f'leiden_{l_param:.2f}'][adata.obs[f'leiden_{l_param:.2f}'].isin(["0,2,0,2", "0,2,0,3", "0,2,0,4","0,2,0,5", "0,2,0,6"])]="0,2,0,2"
    sc.tl.leiden(adata, restrict_to=(f'leiden_{l_param:.2f}', ["0,2,0,0"]),  resolution=0.5, key_added=f'leiden_{l_param:.2f}')
    adata.obs[f'leiden_{l_param:.2f}'][adata.obs[f'leiden_{l_param:.2f}'].isin(["0,2,1", "0,2,0,0,0"])]="0,2,1"
    adata.obs[f'leiden_{l_param:.2f}'][adata.obs[f'leiden_{l_param:.2f}'].isin(["0,2,0,0,1", "0,2,0,0,2", "0,2,0,0,3"])]="0,2,0,0,1"
    adata.obs[f'leiden_{l_param:.2f}'][adata.obs[f'leiden_{l_param:.2f}'].isin(["0,0", "0,1", "0,2,0,2", "0,2,0,0,1"])]="0,0"
    
    adata.obs[f'leiden_{l_param:.2f}'] = adata.obs[f'leiden_{l_param:.2f}'].cat.remove_categories([ "0,1", "0,2,0,0,0", "0,2,0,0,1", "0,2,0,0,2", "0,2,0,0,3","0,2,0,2"])
    sc.tl.leiden(adata, restrict_to=(f'leiden_{l_param:.2f}', ["5"]),  resolution=0.1, key_added=f'leiden_{l_param:.2f}')
    adata.obs[f'leiden_{l_param:.2f}'][adata.obs[f'leiden_{l_param:.2f}'].isin(["5,1", "5,2", "5,3"])]="5,1"
    adata.obs[f'leiden_{l_param:.2f}'] = adata.obs[f'leiden_{l_param:.2f}'].cat.remove_categories([ "5,2", "5,3"])
    adata.obs[f'leiden_{l_param:.2f}'] = adata.obs[f'leiden_{l_param:.2f}'].cat.rename_categories(np.arange(len(np.unique(adata.obs[f'leiden_{l_param:.2f}']))).astype('str'))
    a = list(adata.obs.loc[((adata.obs["condition"].isin(["CD-no-AOM-DSS", "HFD-no-AOM-DSS"])) & (adata.obs[f"leiden_{l_param:.2f}"].isin(["1"]))), f"leiden_{l_param:.2f}"].index)
    shuffle(a)
    print(len(a))
    adata.obs.loc[a[:250], f"leiden_{l_param:.2f}"] = "0"
    # sc.pl.umap(adata, color=[f"leiden_{l_param:.2f}"], cmap="tab20",  show=True)
    
    l_param = f"{l_param:.2f}"
    rows, cols = (1, 4)
    fig, ax = plt.subplots(rows, cols, figsize=(20,20))

    for ind, row in meta.iterrows():
        
        fig_row, fig_col = int(ind/cols), ind%cols
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata_raw = utils.read_filtered_visium_sample(sample_id)
        adata_temp = adata[adata.obs["condition"]==condition,:]    
        adata_temp.obs.index = pd.Index("-".join(cl.split("-")[:-1]) for cl in adata_temp.obs.index.values)
        adata_raw = adata_raw[adata_temp.obs.index,:]
        print(adata_raw)
        adata_raw = adata_raw[:, adata_raw.var_names.isin(adata.var_names)]
        adata_raw.obs[f"leiden_{l_param}"] = adata_temp.obs[f"leiden_{l_param}"]
        sc.pl.spatial(adata_raw, img_key="hires", color =f"leiden_{l_param}", palette=palette_color_lst, title=condition, size=1.50, alpha_img=0.3, ax = ax[ind], show=False)
    fig.tight_layout()
    plt.savefig(f"{PLOT_PATH}/{sample_type}_res-{l_param}_clusters.pdf")
    plt.show();
    sc.pl.umap(adata, color=[f"leiden_{l_param}"], palette=palette_color_lst, title=" Clusters - Integrated Samples",  show=True, save=f'{sample_type}_res-{l_param}_clusters')
    adata.obs[[f"leiden_{l_param}","condition"]].to_csv(os.path.join(output_path, f'{sample_type}_cluster_membership.csv'))

    break
 
print(f"Saving the object... {sample_type}_integrated_clustered.h5ad...")
# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_integrated_clustered.h5ad'))


# python vis_cluster.py -i ../data/out_data/visium_integrated.h5ad -o ../data/out_data -an visium_cluster
# python vis_cluster_annotate.py -i ../data/out_data/visium_integrated_clustered.h5ad -o ../data/out_data -an visium_cluster
