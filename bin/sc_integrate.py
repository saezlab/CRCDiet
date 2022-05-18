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

warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0

# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run intergration by Harmony')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

sample_type ="sc"
input_path = args['input_path']
output_path = args['output_dir']
###############################

S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
PLOT_PATH =  os.path.join(S_PATH, "../plots", "integrate")

Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = PLOT_PATH

sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=150)

print("Reading merged object...")
# Read merged object
adata = sc.read_h5ad(input_path)

print("Running harmony ...")
# Run harmony
sce.pp.harmony_integrate(adata, 'batch', adjusted_basis='X_pca', max_iter_harmony=30)

print("Computing neighbours ...")
# Run umap with updated connectivity
sc.pp.neighbors(adata)
sc.tl.umap(adata)

mpl.rcParams['figure.dpi']= 150
mpl.rcParams["figure.figsize"] = (10,10)
mpl.rcParams["legend.fontsize"]  = 'xx-small'
mpl.rcParams["legend.loc"]  = "upper right"
mpl.rcParams['axes.facecolor'] = "white"



# the number of genes expressed in the count matrix
sc.pl.umap(
    adata, color=["condition", "n_genes_by_counts"], color_map =plt.cm.afmhot, 
    title= ["Condition", "Num of exp. genes"], s=10, frameon=False, ncols=2,  show=True, save=f"{sample_type}_all_condition_harmony"
)

"""rows = 2
columns = 2
grid = plt.GridSpec( rows, columns, wspace = .4, hspace = .4)
plot_list = ["condition", "doublet_score", "n_genes_by_counts", "pct_counts_mt"]
for i in range(rows * columns):
    plt.subplot(grid[i])
    c_ax = plt.gca()
    sc.pl.umap(adata, ax=c_ax, color=[plot_list[i]], color_map =plt.cm.afmhot, frameon=False, show=True)"""

print("Saving the integrated object...")
# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_integrated.h5ad'))

#  python integrate.py -i ../data/out_data/sc_merged.h5ad -o ../data/out_data
