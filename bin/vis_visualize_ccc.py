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


############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=50, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Visualize the mean activities top CCCs')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-cf', '--ccc_file', help='Input ccc file (liana_res)', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
ccc_file = args['ccc_file']
analysis_name = args['analysis_name']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 

print("Clustering the visium slides...")
sample_type = "visium_aomdss"
meta = utils.get_meta_data(sample_type)

df_liana_res = pd.read_csv(ccc_file)
print(df_liana_res)

uniq_pairs = []
for ind, row in df_liana_res.iterrows():
    str_temp = row["ligand_complex"]+"_"+row["receptor_complex"]
    if str_temp not in uniq_pairs:
        uniq_pairs.append(str_temp.split("_"))
    if len(uniq_pairs)==25:
        break
print(uniq_pairs)

for ind, row in meta.iterrows():
    sample_id = row["sample_id"]
    condition = row["condition"]
    row = meta[meta["sample_id"]==sample_id]
    adata = sc.read_h5ad(os.path.join(input_path,f"{sample_id}_filtered.h5ad"))
    adata.var.index = pd.Index(gen.upper() for gen in adata.var.index.values)
    adata.X = adata.layers["log1p_transformed"].copy()
    mpl.rcParams['axes.titlesize'] = 12
    fig, axs = plt.subplots(5, 5, figsize=[20, 20])
    for ind, pair in enumerate(uniq_pairs):
        fig_row, fig_col = int(ind/5), ind%5
        
        if len(pair)==len(set(adata.var.index)&set(pair)):
            adata.obs['my_score2'] = adata[:,pair].X.mean(1)
            sc.pl.spatial(adata, img_key="hires", color = ["my_score2"], title="_".join(pair), cmap="magma_r", size=1.25, alpha_img=0.5, show=False, ax=axs[fig_row][fig_col])
        else:
            sc.pl.spatial(adata, img_key="hires", title="_".join(pair), cmap="magma_r", size=1.25, alpha_img=0.5, show=False, ax=axs[fig_row][fig_col])

        # axs[fig_row][fig_col].plot(cond_list, data[ind,:], marker='o', mec = 'grey', mfc = 'grey', markersize=12, linewidth=5, color=cmap.colors[ind], label=c_t)
        # axs[fig_row][fig_col].set_title(c_t)

        # sc.tl.score_genes(adata, pair, score_name='my_score')
        
    
        # sc.pl.spatial(adata, img_key="hires", color = ["my_score", "my_score2"], cmap="magma_r", size=1.25, alpha_img=0.5, show=True)
    plt.savefig(f"{PLOT_PATH}/{sample_id}_ccc.pdf")
    # plt.show()
    
    # cbar = ax[ind2].collections[0].colorbar
    # cbar.set_ticks([])
    # GNAI2->S1PR1

# adata.obs['my_score'] = adata[:,my_genes].X.sum(1)

"""mpl.rcParams["image.cmap"]= plt.cm.magma_r
            mpl.rcParams['axes.titlesize'] = 12
            # colorbar_loc=None,
            sc.pl.spatial(adata_raw, img_key="hires", color =marker, title=f"{marker} : {condition}", size=1.25, alpha_img=0.5, ax = ax[ind2], show=False)
            cbar = ax[ind2].collections[0].colorbar
            cbar.set_ticks([])"""


# python vis_visualize_ccc.py  -i ../data/out_data -o ../data/out_data -cf ../data/out_data/atlas_bcell1_neutrophils_igap_ccc.csv -an visium_ccc