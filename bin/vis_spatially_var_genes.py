import os
import utils
import pickle
import numpy as np
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
import argparse
import warnings
from pathlib import Path
import matplotlib as mpl
import math
import matplotlib.pyplot as plt
from utils import printmd
from utils import get_meta_data
import squidpy as sq

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run spatially variable genes')
parser.add_argument('-i', '--input_path', help='Input path toobject', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ###############################

sample_type="visium"

meta = get_meta_data(sample_type)
# https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html#Spatially-variable-genes

for ind, row in meta.iterrows():
    

    sample_id = row["sample_id"]
    condition = row["condition"]
    # print(sample_id)
    adata = sc.read_h5ad(os.path.join(OUT_DATA_PATH,f"{sample_id}_filtered.h5ad"))
    # Fetch sample metadata
    m = meta[meta['sample_id'] == sample_id]
    # Add metadata to adata
    for col in m.columns:
        adata.obs[col] = m[col].values[0]
    # Append
    print("burada")
    sq.gr.spatial_neighbors(adata, key_added='spatial')
    print("burada2")
    sc.experimental.pp.highly_variable_genes(
    adata, flavor="pearson_residuals", n_top_genes=3000
    )

    genes = adata[:, adata.var.highly_variable].var_names.values
    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        genes=genes,
        n_perms=100,
        n_jobs=1,
        show_progress_bar=False,
        # backend="multiprocessing"
    )
    print(adata.uns["moranI"].head(10))

    sc.pl.spatial(adata, color=["Guca2a", "Myl9", "Csrp1", "Acta2"], show=True)


    
#python vis_spatially_var_genes.py -i ../data/out_data/ -o ../data/out_data -an "vis_spatially_var_genes"