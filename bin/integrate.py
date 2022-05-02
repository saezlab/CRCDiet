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

sc.set_figure_params(scanpy=True, dpi=150, dpi_save=300)


# Read merged object
adata = sc.read_h5ad(input_path)
print(adata)

rows = 2
columns = 2
grid = plt.GridSpec( rows, columns, wspace = .4, hspace = .4)
# plt.figure(figsize=(100, 100))
# plt.subplot(grid[0, 0:-1])

plot_list = ["condition", "doublet_score", "n_genes_by_counts", "pct_counts_mt"]
for i in range(rows * columns):
    mpl.rcParams['figure.dpi']= 300
    mpl.rcParams["figure.figsize"] = (20,20)
    mpl.rcParams["legend.fontsize"]  = 'xx-small'
    mpl.rcParams["legend.loc"]  = "upper right"
    
    plt.subplot(grid[i])
    # plt.legend(title="deneme", loc="upper left")   
    # plt.annotate('grid '+ str(i), xy = (.5, .5), ha = 'center', 
    #              va = 'center')
    
    c_ax = plt.gca()
    # legend = c_ax.legend(titşeloc='upper center', shadow=True, fontsize='x-small')
    # plt.xlabel('randx')
    # plt.ylabel('randy')
    # legnd = c_ax.get_legend()
    # legnd.title_fontsize = 1.0
    sc.pl.umap(adata, ax=c_ax, color=[plot_list[i]], color_map =plt.cm.afmhot, frameon=False, show=True)
    
    # c_ax.legend(loc="upper center")
    #plt.annotate('grid '+ str(i), xy = (.5, .5), ha = 'center', 
    #             va = 'center')
# fig, axs = plt.subplots(2, 2, constrained_layout=True)
# for ax in axs.flat:


    

#sc.pl.umap(adata, color=["condition"], palette=sc.pl.palettes.default_20, frameon=False, show=True, legend_loc='lower center')
# sc.pl.umap(adata, color=["condition"], palette=sc.pl.palettes.default_20, frameon=False, show=True, legend_loc='upper center')



print("Running harmony ...")
# Run harmony
sce.pp.harmony_integrate(adata, 'batch', adjusted_basis='X_pca', max_iter_harmony=30)

print("Computing neighbors ...")
# Run umap with updated connectivity
sc.pp.neighbors(adata)
sc.tl.umap(adata)


# Plot HVG filtering QC plots
#  fig = plt.figure(figsize=(12,6), dpi=150, tight_layout=True, facecolor='white')
# fig.suptitle('HVG filtering QC plots', fontsize=11)
# gs = fig.add_gridspec(2, 2)

# ax = fig.add_subplot(gs[0,0])

sc.pl.umap(
    adata, color=["condition", "doublet_score", "n_genes_by_counts", "pct_counts_mt"], color_map =plt.cm.afmhot, frameon=False,  show=True, save=f"{sample_type}_all_condition_harmony"
)


sc.pl.umap(adata, color='doublet_score', frameon=False, show=True, save=f"{sample_type}_doubletscore_harmony"
)
# ax = fig.add_subplot(gs[0,1])
sc.pl.umap(
    adata, color=["condition"], palette=sc.pl.palettes.plasma,  show=True, save=f"{sample_type}_condition_harmony"
)


# ax = fig.add_subplot(gs[1,0])
sc.pl.umap(adata, color='n_genes_by_counts', show=True, save=f"{sample_type}_ngenesbycounts_harmony"
)

# ax = fig.add_subplot(gs[1,1])
sc.pl.umap(adata, color='pct_counts_mt',  frameon=False, show=True, save=f"{sample_type}_pctcountsmt_harmony"
)

# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_integrated.h5ad'))

#  python integrate.py -i ../data/out_data/sc_merged.h5ad -o ../data/out_data
