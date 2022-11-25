import matplotlib.pyplot as plt
import scanpy as sc
import argparse
from sklearn.metrics import silhouette_score, pairwise_distances
import sys
import warnings
import utils
import os
from utils import printmd
import matplotlib as mpl

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Atlas cell type annotation')
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



l_param = 0.40
cluster_to_cell = [
"Epithelial cells (SI & colon)",
"Epithelial cells in (SI Enterocytes)",
"Fibroblasts",
"Epithelial cells (colon?)",
"Cells with high mitochondrial activity",
"IECs in colon (Goblet)",
"Epithelial cells (colon?)",
"Epithelial cells (SI & colon)",
"Epithelial cells (SI & colon)",
"B cells",
"T cells (and B cells)",
"Epithelial??",
"Immune??",
"Epithelial cells (goblet)",
"Epithelial cells (colon)",
"Goblet cells",
"Antigen presenting cells (Macrophages?)",
"T cells",
"neuronal cells (enteric neurons)?",
"Immune Cells?",
"Epithelial cells (SI & colon)",
"Plasma cells",
"Endothelial cells?",
"Endothelial cells?",
"neuronal cells (enteric neurons?)",
"??",
"epithelial?neuronal?",
"Enterocytes",
"Paneth cells",
"Tuft cells",
"T cells/NK cells",
"neuronal cells (enteric neurons?)",
"Macrophages",
"Enteroendocine cells",
"Smooth muscle cells",
"Enteroendocine cells",
"Epithelial cells (SI & colon)",
"Immune cells: Granulocytes",
"Epithelial cells (colon)",
"??",
"epithelial cells (goblet)",
"Epithelial cells (colon)",
"neuronal cells (enteric neurons?)",
"epithelial cells (goblet)",
"B cells (Plasma cells)",
"??",
"Keratinocytes (distal colonic transitional epithelium??)",
"B cells (Plasma cells)",
"Immune cells",
"B cells (Plasma cells)",
"Paneth cells (secretory IEC)",
"Immune cells + other cells??",
"Paneth cells (secretory IEC)",
"??",
"??"]

print(len(cluster_to_cell))
adata = sc.read_h5ad(input_path)
annotation_dict = dict()
for ind, val in enumerate(cluster_to_cell):
    annotation_dict[str(ind)] = val
print(annotation_dict)

print(adata.obs)
print([clust for clust in adata.obs[f'leiden_{l_param:.2f}']])
print(annotation_dict["54"])
adata.obs['cell_type'] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{l_param:.2f}']]
plt.rcParams['figure.dpi']= 300
plt.rcParams['figure.figsize']= (45, 30)
sc.pl.umap(adata, color=f'cell_type', title= ["Cell type annotation"], show=True, s=10, legend_loc="on data", legend_fontsize="xx-small",  save=f'{sample_type}_cell_type_annot_{l_param:.2f}')

# python atlas_cell_type_annotation.py -i ../data/out_data/atlas_integrated_clustered.h5ad -o ../data/out_data -st atlas  -an atlas_cluster