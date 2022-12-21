from reprlib import aRepr
from pathlib import Path
from imageio import save
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

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run marker visualization')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # "visium_visualize_markers"
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ###############################

sample_type ="visium"
adata_integ_clust = sc.read_h5ad(input_path)


meta = utils.get_meta_data(sample_type)
condition = np.unique(meta['condition'])
# adata = sc.read_h5ad(input_path)

adata = utils.get_filtered_concat_data(sample_type)

# filter out the cells missing in adata_integ_clust
adata = adata[adata_integ_clust.obs_names,:]

adata.obsm["X_umap"] = adata_integ_clust.obsm["X_umap"]
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

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
    plt.savefig(f'{PLOT_PATH}/{sample_type}_{marker}.pdf');
    plt.show();
    

    sc.pl.umap(adata, color=marker, show=True, save=f"{sample_type}_{marker}.pdf");
