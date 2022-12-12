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

adata = sc.read_h5ad(input_path)

adata_concat = utils.get_filtered_concat_data(sample_type)
"""for item in adata_conccat.var_names:
    if item.startswith("Pd"):
        print(item)"""
adata_concat = adata_concat[adata.obs_names,:]

sc.pp.normalize_total(adata_concat, target_sum=1e6)
sc.pp.log1p(adata_concat)
# https://github.com/scverse/scanpy/issues/2239
# if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
# KeyError: 'base'
adata.uns['log1p']["base"] = None


marker_genes = dc.get_resource('PanglaoDB')
marker_genes = marker_genes[(marker_genes['mouse']=='True')&(marker_genes['canonical_marker']=='True')]
marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]


marker_genes["genesymbol"] = marker_genes["genesymbol"].str.capitalize()
marker_genes["cell_type"] = marker_genes["cell_type"].str.capitalize()
marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]
adata.obsm["X_umap"] = adata_integ_clust.obsm["X_umap"]


dc.run_ora(mat=adata, net=marker_genes, use_raw=False, source='cell_type', target='genesymbol', min_n=3, verbose=True)
acts = dc.get_acts(adata, obsm_key='ora_estimate')

for c_t in list(adata.obsm["ora_estimate"].columns):
    adata.obs[c_t] = adata.obsm["ora_estimate"][c_t]


dict_mean_enr = dict()

mean_enr = dc.summarize_acts(acts, groupby=f'leiden_{l_param}')
# print(mean_enr)

annotation_dict = dc.assign_groups(mean_enr)
# print(annotation_dict)
# Manual annotation
adata.obs[f'cell_type_{l_param}'] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{l_param}']]



plt.rcParams['figure.dpi']= 300
plt.rcParams['figure.figsize']= (15, 10)
# the number of genes expressed in the count matrix
sc.pl.umap(adata, color=f'cell_type_{l_param}', title= ["Cell type annotation"], show=True, s=10, legend_loc="on data", legend_fontsize="xx-small",  save=f'{sample_type}_cell_type_annot_{l_param}')
sc.pl.umap(adata, color='major_cell_types', title= ["Major cell types"], show=True, s=10, legend_loc="on data", legend_fontsize="xx-small",  save=f'{sample_type}_major_cell_types')
adata_integ_clust.obs[f'cell_type_{l_param}'] = adata.obs[f'cell_type_{l_param}']
adata_integ_clust.obs['major_cell_types'] = adata.obs['major_cell_types']



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