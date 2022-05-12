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
import matplotlib as mpl
import argparse
import os
from sklearn.metrics import silhouette_score, pairwise_distances
import pickle
import utils
import math
import warnings


warnings.simplefilter(action='ignore')

# Read integrated object
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Run Clustering and annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']

###############################

sample_type = "visium"
meta = utils.get_meta_data("visium")

S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
PLOT_PATH =  os.path.join(S_PATH, "../plots", "cluster")

Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = PLOT_PATH
sc.set_figure_params(scanpy=True, facecolor="white", dpi=50, dpi_save=50)
adata = sc.read_h5ad(input_path)

markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
markers = list(set(markers_df["genesymbol"].str.capitalize()))

best_l_param, best_s_scr = -1, -1

step = 0.05
dist_mat = None

print("Calculating distance matrix... ")
dist_mat = pairwise_distances(adata.obsm["X_pca"], metric='euclidean') #  , n_jobs=-1)

with open(f'{OUT_DATA_PATH}/{sample_type}_dist_mat.pickle', 'wb') as handle:
        pickle.dump(dist_mat, handle, protocol=pickle.HIGHEST_PROTOCOL)

"""if not os.path.exists(f'{OUT_DATA_PATH}/{sample_type}_dist_mat.pickle'):
    print("Calculating distance matrix... ")
    dist_mat = pairwise_distances(adata.obsm["X_pca"], metric='euclidean') #  , n_jobs=-1)

    with open(f'{OUT_DATA_PATH}/{sample_type}_dist_mat.pickle', 'wb') as handle:
        pickle.dump(dist_mat, handle, protocol=pickle.HIGHEST_PROTOCOL)

else:
    print("Loading precomputed distance matrix... ")
    with open(f'{OUT_DATA_PATH}/{sample_type}_dist_mat.pickle', 'rb') as handle:
        dist_mat = pickle.load(handle) """


print("Computing neighbourhood graph... ")
sc.pp.neighbors(adata)



silh_param_scores = []
# perform clustering, Rank genes for characterizing groups, plot top 5 genes
for l_param in np.arange(0.1, 1.01, 0.05):

    print(f"Creating clusters with Leiden resolution param: {l_param:.2f}")
    sc.tl.leiden(adata, resolution = l_param, key_added = f"leiden_{l_param:.2f}") # default resolution in 1.0
    silh_scr = silhouette_score(dist_mat, np.array(adata.obs[f"leiden_{l_param:.2f}"]), metric='precomputed')
    # print(f"Clustering param: {l_param:.2f}\tSilhoutte score: {silh_scr:.3f}")
    silh_param_scores.append((l_param, silh_scr))
    # sc.pl.umap(adata, color=[f"leiden_{l_param:.2f}"], cmap="tab20",  show=True, save=f'{sample_type}_res-{l_param:.2f}_silh-{silh_scr:.3f}_clusters')

    """l_param = f"{l_param:.2f}"
    rows, cols = (1, 6)
    fig, ax = plt.subplots(rows, cols, figsize=(20,20))

    for ind, row in meta.iterrows():
        
        fig_row, fig_col = int(ind/cols), ind%cols
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata_raw = utils.read_raw_visium_sample(sample_id)
        adata_temp = adata[adata.obs["condition"]==condition,:]    
        adata_temp.obs.index = pd.Index("-".join(cl.split("-")[:-1]) for cl in adata_temp.obs.index.values)
        adata_raw = adata_raw[adata_temp.obs.index,:]
        adata_raw = adata_raw[:, list(adata.var_names)]
        adata_raw.obs[f"leiden_{l_param}"] = adata_temp.obs[f"leiden_{l_param}"]
        sc.pl.spatial(adata_raw, img_key="hires", color =f"leiden_{l_param}",  title=condition, size=1.25, alpha_img=0.3, ax = ax[ind], show=False)
    fig.tight_layout()
    plt.show();"""

for l_param, s_scr in silh_param_scores:
    if s_scr > best_s_scr:
        best_l_param, best_s_scr = l_param, s_scr

adata.uns["leiden_best_silh_param"] = [best_l_param, best_s_scr]

l_param, _ = adata.uns["leiden_best_silh_param"]

l_param = f"{l_param:.2f}"

sc.pl.umap(adata, color=[f"leiden_{l_param}"], title=" Clusters - Integrated Samples", palette=sc.pl.palettes.default_20, show=True, save=f'{sample_type}_res-{l_param}_clusters')

rows, cols = (1, 6)
fig, ax = plt.subplots(rows, cols, figsize=(20,20))

for ind, row in meta.iterrows():
    
    fig_row, fig_col = int(ind/cols), ind%cols
    sample_id = row["sample_id"]
    condition = row["condition"]
    adata_raw = utils.read_raw_visium_sample(sample_id)
    adata_temp = adata[adata.obs["condition"]==condition,:]    
    adata_temp.obs.index = pd.Index("-".join(cl.split("-")[:-1]) for cl in adata_temp.obs.index.values)
    adata_raw = adata_raw[adata_temp.obs.index,:]
    adata_raw = adata_raw[:, list(adata.var_names)]
    mpl.rcParams["legend.fontsize"]  = 'xx-small'
    adata_raw.obs[f"leiden_{l_param}"] = adata_temp.obs[f"leiden_{l_param}"]
    sc.pl.spatial(adata_raw, img_key="hires", color =f"leiden_{l_param}", title=condition, size=1.25, alpha_img=0.5, ax = ax[ind], show=False)
    # plt.tight_layout(pad=3.0)

    #sc.pl.violin(adata, list(set(adata.var.index) & set(markers)), show=True, groupby=f"leiden_{l_param}")

    
    #mpl.rcParams["image.cmap"]= plt.cm.magma_r
    #sc.pl.spatial(adata_raw, img_key="hires", color =list(set(adata.var.index) & set(markers)),  size=1.0, alpha_img=0.5, wspace = 0.3)
plt.tight_layout()
plt.show();

print(f"Saving the object... {sample_type}_integrated_clustered.h5ad...")
# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_integrated_clustered.h5ad'))



    

