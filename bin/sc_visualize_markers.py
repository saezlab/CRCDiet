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

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=150)
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run marker visualization')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths("sc_visualize_markers")
############################### BOOOORIING STUFF ABOVE ###############################

sample_type ="sc"
adata_integ_clust = sc.read_h5ad(input_path)

meta = utils.get_meta_data(sample_type)
# adata = sc.read_h5ad(input_path)

adata = utils.get_filtered_concat_data(sample_type)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

adata.obsm["X_umap"] = adata_integ_clust.obsm["X_umap"]

markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
markers = list(set(markers_df["genesymbol"].str.capitalize()))

print("Plotting the activities of marker genes on the slides...\n")

marker_intersect = list(set(adata.var.index) & set(markers))
print(f"Number of marker genes: {len(marker_intersect)}")

marker_ind = 0
while marker_ind<len(marker_intersect):
    mrk_str = ",".join(marker_intersect[marker_ind:marker_ind+4])
    print(f"Plotting markers: {mrk_str}")
    sc.pl.umap(adata, color=marker_intersect[marker_ind:marker_ind+4], ncols=len(marker_intersect[marker_ind:marker_ind+4]), save=f'{sample_type}_marker_{mrk_str}')
    marker_ind += 4
    
    
    # fig.tight_layout()
    
 

#  python sc_visualize_markers.py -i ../data/out_data/sc_integrated_clustered.h5ad -o ../data/out_data