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
parser = argparse.ArgumentParser(prog='qc', description='Run marker visualization')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
parser.add_argument('-ct', '--cell_type', help='Major cell type', required=True)
args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # "sc_cluster_annotate"
m_ct = args['cell_type']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ###############################

sample_type ="sc"

epithelial_samples = ["CD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN"]
immune_samples = ["CD-AOM-DSS-Immune", "HFD-AOM-DSS-Immune", "LFD-AOM-DSS-Immune"]

mct_str = "_".join(m_ct.split(" ")).lower()
adata_integ_clust = sc.read_h5ad(input_path)

"""if "epithelial" in mct_str:
   adata_integ_clust = adata_integ_clust[adata_integ_clust.obs["condition"].isin(epithelial_samples)]
elif "immune" in mct_str:
    adata_integ_clust = adata_integ_clust[adata_integ_clust.obs["condition"].isin(immune_samples)]"""

meta = utils.get_meta_data(sample_type)
# adata = sc.read_h5ad(input_path)

adata = utils.get_filtered_concat_data(sample_type)

adata = adata[adata_integ_clust.obs_names,:]
# keep raw counts in layers
adata.layers['counts'] = adata.X.copy()
# adata.layers["sqrt_norm"] = np.sqrt(
#    sc.pp.normalize_total(adata, inplace=False)["X"]).copy()

adata.layers["log1p_transformed"] = sc.pp.normalize_total(adata, inplace=False, target_sum=1e6)["X"]
sc.pp.log1p(adata, layer="log1p_transformed")

markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
markers = list(set(markers_df["genesymbol"].str.capitalize()))

# for m_ct in set(adata_integ_clust.obs["major_cell_types"]):
adata_tmp = None
if m_ct in ["Epithelial cells", "Immune cells", "Stroma cells"]:
    adata_tmp = adata[adata_integ_clust.obs["major_cell_types"]==m_ct,:].copy()
else:
    adata_tmp = adata[adata_integ_clust.obs["cell_type_0.20"]==m_ct,:].copy()



# adata_tmp = adata[adata_integ_clust.obs["major_cell_types"]==m_ct,:].copy()
del adata
del adata_integ_clust
sc.pp.calculate_qc_metrics(adata_tmp, inplace=True)

for ind in range(6):
    adata_tmp.var[f'mt-{ind}'] = adata_tmp.var[f'mt-{ind}'].astype(str)
    adata_tmp.var[f'rp-{ind}'] = adata_tmp.var[f'mt-{ind}'].astype(str) 

sc.experimental.pp.highly_variable_genes(
    adata_tmp, batch_key='batch', flavor="pearson_residuals", n_top_genes=3000
)

adata_tmp.var = adata_tmp.var[['highly_variable','highly_variable_nbatches']]

# Filter by HVG
num_hvg_genes = 4000
batch_msk = np.array(adata_tmp.var.highly_variable_nbatches > 1)
hvg = adata_tmp.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata_tmp.var['highly_variable'] = [g in hvg for g in adata_tmp.var.index]
adata_tmp.var = adata_tmp.var[['highly_variable','highly_variable_nbatches']]

print("Performing analytic Pearson residual normalization...")
sc.experimental.pp.normalize_pearson_residuals(adata_tmp)

# Run PCA
sc.tl.pca(adata_tmp, svd_solver='arpack', random_state=0)

print("Running harmony ...")
# Run harmony
sce.pp.harmony_integrate(adata_tmp, 'batch', adjusted_basis='X_pca', max_iter_harmony=30)

print("Computing neighbours ...")
# Run umap with updated connectivity
sc.pp.neighbors(adata_tmp)
sc.tl.umap(adata_tmp)

mpl.rcParams['figure.dpi']= 300
mpl.rcParams["figure.figsize"] = (10,10)
mpl.rcParams["legend.fontsize"]  = 'xx-small'
mpl.rcParams["legend.loc"]  = "upper right"
mpl.rcParams['axes.facecolor'] = "white"

# the number of genes expressed in the count matrix
sc.pl.umap(
    adata_tmp, color=["condition", "n_genes_by_counts"], color_map =plt.cm.afmhot, 
    title= ["Condition", "Num of exp. genes"], s=15, frameon=False, ncols=2,  show=False, save=f"{sample_type}_{mct_str}_all_condition_harmony"
)



# Write to file
adata_tmp.write(os.path.join(output_path, f'{sample_type}_{mct_str}_integrated.h5ad'))



"""



python sc_major_cell_type_integration.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an "sc_epithelial_cells_integrate" -ct "Epithelial cells"
python sc_major_cell_type_integration.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an "sc_immune_cells_integrate" -ct "Immune cells"
python sc_major_cell_type_integration.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an "sc_stroma_cells_integrate" -ct "Stroma cells"


python sc_cluster.py -i ../data/out_data/sc_epithelial_cells_integrated.h5ad -o ../data/out_data -an sc_epithelial_cells_cluster -of epithelial_cells
python sc_cluster.py -i ../data/out_data/sc_immune_cells_integrated.h5ad -o ../data/out_data -an sc_immune_cells_cluster -of immune_cells
python sc_cluster.py -i ../data/out_data/sc_stroma_cells_integrated.h5ad -o ../data/out_data -an sc_stroma_cells_cluster -of stroma_cells

python sc_cluster_annotate.py -i ../data/out_data/sc_epithelial_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_epithelial_cells_cluster -of epithelial_cells
python sc_cluster_annotate.py -i ../data/out_data/sc_epithelial_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_immune_cells_cluster -of immune_cells
python sc_cluster_annotate.py -i ../data/out_data/sc_epithelial_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_stroma_cells_cluster -of stroma_cells

python sc_visualize_markers.py -i ../data/out_data/sc_epithelial_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_epithelial_cells_visualize_markers 
python sc_visualize_markers.py -i ../data/out_data/sc_immune_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_immune_cells_visualize_markers 
python sc_visualize_markers.py -i ../data/out_data/sc_stroma_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_stroma_cells_visualize_markers 

python sc_major_cell_type_integration.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an "sc_bcells_integrate" -ct "B cells"
"""
#  python sc_major_cell_type_integration.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data