from reprlib import aRepr
from pathlib import Path
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os
import warnings
import utils

warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0

# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run marker visualization')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

sample_type ="visium"
input_path = args['input_path']
output_path = args['output_dir']
###############################

S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
PLOT_PATH =  os.path.join(S_PATH, "../plots", "visium_visualize_markers")

Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = PLOT_PATH

sc.set_figure_params(scanpy=True, facecolor="white", dpi=80) # , dpi_save=150

adata_integ_clust = sc.read_h5ad(input_path)


meta = utils.get_meta_data("visium")
# adata = sc.read_h5ad(input_path)

adata = []
for ind, row in meta.iterrows():
    sample_id = row["sample_id"]
    condition = row["condition"]
    # print(sample_id)
    tmp = sc.read_h5ad(os.path.join(OUT_DATA_PATH,f"{sample_id}_filtered.h5ad"))
    # Fetch sample metadata
    m = meta[meta['sample_id'] == sample_id]
    # Add metadata to adata
    for col in m.columns:
        tmp.obs[col] = m[col].values[0]
    # Append
    adata.append(tmp)
    del tmp

# Merge objects and delete list
adata = adata[0].concatenate(adata[1:], join='outer')


adata.obsm["X_umap"] = adata_integ_clust.obsm["X_umap"]

markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
markers = list(set(markers_df["genesymbol"].str.capitalize()))

print("Plotting the activities of marker genes on the slides...\n")

marker_intersect = list(set(adata.var.index) & set(markers))
print(f"Number of intersected marker genes: {len(marker_intersect)}")

for ind, marker in enumerate(marker_intersect):
    rows, cols = (len(marker_intersect), 6)
    fig, ax = plt.subplots(1, 6, figsize=(20,2))
    # fig.tight_layout()
    print(f"Plotting marker: {marker}")
    for ind2, row in meta.iterrows():
    
        # fig_row, fig_col = int(ind/cols), ind%cols
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata_raw = utils.read_raw_visium_sample(sample_id)
        adata_temp = adata[adata.obs["condition"]==condition,:]    
        adata_temp.obs.index = pd.Index("-".join(cl.split("-")[:-1]) for cl in adata_temp.obs.index.values)
        adata_raw = adata_raw[adata_temp.obs.index,:]
        adata_raw = adata_raw[:, list(adata.var_names)]

        mpl.rcParams["image.cmap"]= plt.cm.magma_r
        mpl.rcParams['axes.titlesize'] = 12
        # colorbar_loc=None,
        sc.pl.spatial(adata_raw, img_key="hires", color =marker, title=f"{marker} : {condition}", size=1.25, alpha_img=0.5, ax = ax[ind2], show=False)
        cbar = ax[ind2].collections[0].colorbar
        cbar.set_ticks([])
        #plt.tight_layout(h_pad=1)
    
        # sc.pl.violin(adata, list(set(adata.var.index) & set(markers)), show=True, groupby=f"leiden_{l_param}")
    plt.show();

    sc.pl.umap(adata, color=marker)

#  python vis_visualize_cluster.py -i ../data/out_data/visium_integrated_clustered.h5ad -o ../data/out_data