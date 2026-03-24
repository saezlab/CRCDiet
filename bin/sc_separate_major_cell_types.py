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
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Create new object')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
parser.add_argument('-ct', '--cell_type', help='Major cell type', required=True)
args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # subset
m_ct = args['cell_type']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name, plot=False)
############################### BOOOORIING STUFF ABOVE ###############################

# save the cells coming from major cell type as a seperate object
adata= sc.read_h5ad(input_path)

adata = adata[adata.obs["major_cell_types"]==m_ct]

adata_raw = utils.get_filtered_concat_data("sc")
adata_raw = adata_raw[adata.obs_names,:]
adata_raw.obs["cell_type_0.20"] = adata.obs["cell_type_0.20"]
adata_raw.obsm["X_umap"] = adata.obsm["X_umap"]
adata_raw.obsm["X_pca"] = adata.obsm["X_pca"]


for ind in range(6):
    adata_raw.var[f'mt-{ind}'] = adata_raw.var[f'mt-{ind}'].astype(str)
    adata_raw.var[f'rp-{ind}'] = adata_raw.var[f'mt-{ind}'].astype(str)

adata_raw.write(f"{output_path}/sc_{m_ct}.h5ad")
