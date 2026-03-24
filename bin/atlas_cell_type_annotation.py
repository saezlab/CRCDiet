import matplotlib.pyplot as plt
import scanpy as sc
import argparse
from sklearn.metrics import silhouette_score, pairwise_distances
import sys
import warnings
import utils
import os
from utils import printmd
import decoupler as dc
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
parser.add_argument('-st', '--sample_type', default="atlas", help='Sample type', required=False)
args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # sc_integrate
sample_type = args['sample_type']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 



l_param = 0.40


adata_ref = utils.get_filtered_concat_data(sample_type)

for col in adata_ref.obs.columns:
    if col.startswith(("mt", "rp", "mean", "pct", "log", "n_cells", "log1p", "gene", "total", "feature")):
        del adata_ref.obs[col]

for col in adata_ref.var.columns:
    if col.startswith(("mt", "rp", "mean", "pct", "log", "n_cells", "log1p", "gene", "total", "feature")):
        del adata_ref.var[col]

adata_ref.X = adata_ref.layers["counts"]
del adata_ref.layers["counts"]
del adata_ref.layers["raw"]
del adata_ref.layers["sqrt_norm"]

adata = sc.read_h5ad(input_path)


adata_ref.obs["leiden_0.40"] = adata.obs["leiden_0.40"]
adata_ref.obsm["X_umap"] = adata.obsm["X_umap"]


annotation_dict_level0 = {"0": "Epithelial cells-colon", "1": "Epithelial cells-SI enterocytes", "2": "Epithelial cells-SI enterocytes", "3": "Goblet cells", "4": "Epithelial cells-colon", "5": "Fibroblasts", "6": "", "6,0": "Prolif. cells", "6,1": "Prolif. cells", "6,2": "Prolif. cells", "6,3": "Prolif. cells", "6,4": "Prolif. cells", "7": "", "7,0": "Endothelial cells", "7,1": "T cells", "7,2": "T cells-NK cells", "7,3": "ILC2", "7,4": "Epithelial cells", "7,5": "T cells", "7,6": "T cells-NK cells", "8": "Epithelial cells", "9": "Goblet cells", "10": "Goblet cells", "11": "Epithelial cells", "12": "B cells", "13": "", "13,0": "Goblet cells", "13,1": "Goblet cells", "13,2": "Goblet cells", "13,3": "Goblet cells", "14": "Fibroblasts", "15": "T cells", "16": "Neuronal cells", "17": "", "17,0": "Macrophages", "17,1": "Macrophages", "17,2": "Macrophages", "17,3": "Macrophages", "17,4": "IgA plasma cells", "17,5": "B cells", "17,6": "Epithelial cells?", "17,7": "B cells", "17,8": "B cells", "17,9": "Macrophages-DCs", "18": "Epithelial cells", "19": "Endothelial cells", "20": "Smooth muscle cells", "21": "Epithelial cells", "22": "IgA plasma cells", "23": "Epithelial cells", "24": "Epithelial cells-SI enterocytes", "25": "Paneth cells", "26": "B cells", "27": "Goblet cells", "28": "", "28,0": "Endothelial cells", "28,1": "Endothelial cells", "28,2": "Endothelial cells", "28,3": "Endothelial cells", "28,4": "Fibroblasts", "29": "T cells", "30": "", "30,0": "Tuft cells", "30,1": "Tuft cells", "30,2": "Tuft cells", "30,3": "Tuft cells", "30,4": "Tuft cells", "30,5": "Tuft cells", "31": "Goblet cells", "32": "B cells", "33": "Enteroendocrine cells", "34": "Neuronal cells", "35": "Smooth muscle cells", "36": "", "36,0": "Mast cells", "36,1": "Neutrophils", "36,2": "Mast cells", "36,3": "Mast cells", "36,4": "Neutrophils", "36,5": "Neutrophils", "36,6": "Neutrophils", "37": "Epithelial cells", "38": "", "38,0": "Epithelial cells-colon", "38,1": "B cells", "38,2": "Epithelial cells", "38,3": "Epithelial cells", "39": "Paneth-Goblet cells", "40": "Goblet cells", "41": "Paneth cells", "42": "Fibroblasts", "43": "Fibroblasts", "44": "IgA plasma cells", "45": "Epithelial cells", "46": "Keratinocytes", "47": "??", "48": "B cells", "49": "??", "50": "Paneth-Goblet cells"}

annotation_dict_level1 = {"0": "Epithelial cells-colon-1", "1": "Epithelial cells-SI enterocytes-1", "2": "Epithelial cells-SI enterocytes-2", "3": "Goblet cells-1", "4": "Epithelial cells-colon-2", "5": "Fibroblasts-1", "6": "", "6,0": "Prolif. cells", "6,1": "Prolif. cells", "6,2": "Prolif. cells", "6,3": "Prolif. cells", "6,4": "Prolif. cells", "7": "", "7,0": "Endothelial cells-1", "7,1": "T cells", "7,2": "T cells-NK cells", "7,3": "ILC2", "7,4": "Epithelial cells-1", "7,5": "T cells", "7,6": "T cells-NK cells", "8": "Epithelial cells-2", "9": "Goblet cells-2", "10": "Goblet cells-3", "11": "Epithelial cells-3", "12": "B cells-1", "13": "", "13,0": "Goblet cells-4", "13,1": "Goblet cells-4", "13,2": "Goblet cells-4", "13,3": "Goblet cells-4", "14": "Fibroblasts-2", "15": "T cells", "16": "Neuronal cells-1", "17": "", "17,0": "Macrophages", "17,1": "Macrophages", "17,2": "Macrophages", "17,3": "Macrophages", "17,4": "IgA plasma cells-1", "17,5": "B cells-2", "17,6": "Epithelial cells?", "17,7": "B cells-2", "17,8": "B cells-2", "17,9": "Macrophages-DCs", "18": "Epithelial cells-4", "19": "Endothelial cells-2", "20": "Smooth muscle cells-1", "21": "Epithelial cells-5", "22": "IgA plasma cells-2", "23": "Epithelial cells-5", "24": "Epithelial cells-SI enterocytes-3", "25": "Paneth cells-1", "26": "B cells-3", "27": "Goblet cells-5", "28": "", "28,0": "Endothelial cells-3", "28,1": "Endothelial cells-3", "28,2": "Endothelial cells-3", "28,3": "Endothelial cells-3", "28,4": "Fibroblasts-3", "29": "T cells", "30": "", "30,0": "Tuft cells", "30,1": "Tuft cells", "30,2": "Tuft cells", "30,3": "Tuft cells", "30,4": "Tuft cells", "30,5": "Tuft cells", "31": "Goblet cells-6", "32": "B cells-2", "33": "Enteroendocrine cells", "34": "Neuronal cells-2", "35": "Smooth muscle cells-2", "36": "", "36,0": "Mast cells", "36,1": "Neutrophils", "36,2": "Mast cells", "36,3": "Mast cells", "36,4": "Neutrophils", "36,5": "Neutrophils", "36,6": "Neutrophils", "37": "Epithelial cells-6", "38": "", "38,0": "Epithelial cells-colon-3", "38,1": "B cells-5", "38,2": "Epithelial cells", "38,3": "Epithelial cells", "39": "Paneth-Goblet cells", "40": "Goblet cells-7", "41": "Paneth cells-2", "42": "Fibroblasts-4", "43": "Fibroblasts-5", "44": "IgA plasma cells-2", "45": "Epithelial cells", "46": "Keratinocytes", "47": "??", "48": "B cells-2", "49": "??", "50": "Paneth-Goblet cells"}

plt.rcParams['figure.dpi']= 300
plt.rcParams['figure.figsize']= (45, 30)

exluded_clusters = ["17,9", "38,2", "38,3", "45", "47", "49", "50"]
for clust in exluded_clusters:
    del annotation_dict_level0[clust]
    del annotation_dict_level1[clust]


adata = adata[~adata.obs["leiden_0.40"].isin(exluded_clusters)]
adata_ref = adata_ref[adata.obs_names,:]


adata.obs['cell_type_level0'] = [annotation_dict_level0[str(clust)] for clust in adata.obs[f'leiden_0.40']]
adata.obs['cell_type_level1'] = [annotation_dict_level1[str(clust)] for clust in adata.obs[f'leiden_0.40']]

adata_ref.obs['cell_type_level0'] = adata.obs['cell_type_level0']
adata_ref.obs['cell_type_level1'] = adata.obs['cell_type_level1']

sc.pl.umap(adata, color=f'cell_type_level0', title= ["Cell type annotation"], show=True, s=10, legend_loc="on data", legend_fontsize="small",  save=f'atlas_cell_type_annot_level_0_ondata')
sc.pl.umap(adata, color=f'cell_type_level1', title= ["Cell type annotation"], show=True, s=10, legend_loc="on data", legend_fontsize="small",  save=f'atlas_cell_type_annot_level_1_ondata')
sc.pl.umap(adata, color=f'cell_type_level0', title= ["Cell type annotation"], show=True, s=10, legend_fontsize="small",  save=f'atlas_cell_type_annot_level_0')
sc.pl.umap(adata, color=f'cell_type_level1', title= ["Cell type annotation"], show=True, s=10, legend_fontsize="small",  save=f'atlas_cell_type_annot_level_1')

adata_ref.write(os.path.join(OUT_DATA_PATH, f"{sample_type}_cell_type_annot_light_weight.h5ad"))


# python atlas_cell_type_annotation.py -i ../data/out_data/atlas_integrated_subclustered.h5ad -o ../data/out_data -st atlas  -an atlas_cell_type_annot