
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os
import utils
import warnings
import scvi

'''
Integrate the merged samples using Harmony and save the AnnData object
'''
############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run intergration by Harmony')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)
args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # sc_integrate
sample_type = args['sample_type']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 

print("Reading merged object...")
# Read merged object
adata = sc.read_h5ad(input_path)
print(f"Number of cells: {adata.shape[0]}")

# scvi.model.SCVI.setup_anndata(adata,layer="counts", batch_key = 'batch')
scvi.model.SCVI.setup_anndata(adata,layer="counts", batch_key = 'batch')
# scvi.data.view_anndata_setup(adata)

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train(max_epochs=200)

# model, trials = scvi.inference.autotune.auto_tune_scvi_model("cortex", adata)

latent = model.get_latent_representation()

adata.obsm["X_scVI"] = latent
print("calculating neighbors...")
# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI")
print("calculating neighbors...")
sc.tl.umap(adata, min_dist=0.3)
print("performing leiden clustering neighbors...")
sc.tl.leiden(adata, key_added="leiden_scVI", resolution=0.5)




plt.rcParams['figure.dpi']= 300
plt.rcParams["figure.figsize"] = (10,10)
plt.rcParams["legend.fontsize"]  = 'xx-small'
plt.rcParams["legend.loc"]  = "upper right"
plt.rcParams['axes.facecolor'] = "white"

"""# the number of genes expressed in the count matrix
sc.pl.umap(
    adata, color=["condition", "n_genes_by_counts"], color_map =plt.cm.afmhot, 
    title= ["Condition", "Num of exp. genes"], s=10, frameon=False, ncols=2,  show=True, save=f"{sample_type}_all_condition_harmony"
)"""

plt.rcParams['figure.dpi']= 300
plt.rcParams['figure.figsize']= (45, 30)

sc.pl.umap(
    adata, color="condition",
    title= "Condition", size=10, frameon=False, show=True, save=f"{sample_type}_all_condition_scvi"
)

sc.pl.umap(
    adata,
    color=["leiden_scVI"],
    title= "Cluster", size=10, frameon=False, show=True, save=f"{sample_type}_all_condition_cluster_scvi",
)

print("Saving the integrated object...")
# Write to file
# adata.write(os.path.join(output_path, f'{sample_type}_integrated_scvi.h5ad'))

#  python sc_integrate_scvi.py -i ../data/out_data/atlas_merged.h5ad -o ../data/out_data -st atlas -an atlas_integrate  e
