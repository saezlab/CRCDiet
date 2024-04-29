from copyreg import pickle
from genericpath import sameopenfile
from operator import index
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns 
# import decoupler as dc  
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import os
from sklearn.metrics import silhouette_score, pairwise_distances
import pickle
import warnings
import utils
import matplotlib as mpl

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=50, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Run Clustering and annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
parser.add_argument('-of', '--output_file', help='Output file name', required=False)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)
parser.add_argument('-cc', '--cluster_column', default="cell_type", help='Cluster column', required=True)
parser.add_argument('-res', '--resolution_param', default="0.2", help='Resolution parameter', required=True)
parser.add_argument("-cg", "--cluster_groups", help="Cluster groups to be subclustered", required=True)
parser.add_argument("-save", "--save_adata", action=argparse.BooleanOptionalAction, default=False, help="Save anndata with new clusters?", required=True)



args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # sc_cluster
output_file = args['output_file']
sample_type = args['sample_type']
cluster_column = args['cluster_column']
resolution_param = args['resolution_param']
cluster_groups= args['cluster_groups']
save_adata = args['save_adata']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 

# plt.rcParams['figure.dpi']= 300
# plt.rcParams['figure.figsize']= (45, 30)

resolution_param = resolution_param.split(",")
if len(resolution_param) == 1:
    resolution_param = float(resolution_param[0])
else:
    resolution_param = [float(res) for res in resolution_param]

print(resolution_param)
nc_genes_df = pd.read_csv("/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/musmusculus_non_coding_genes_ncbi.txt", sep="\t")

all_gene_symbols = []
for ind, row in nc_genes_df.iterrows():
    gene_symbol = row["Symbol"]
    aliases = row["Aliases"]
    if pd.isna(aliases):
        continue
    else:
        aliases = aliases.split(", ")
        aliases = [alias.upper() for alias in aliases]
        all_gene_symbols.extend(aliases)
    if gene_symbol is not None:
        all_gene_symbols.append(gene_symbol.upper())

adata = sc.read_h5ad(input_path)
print(adata)
adata.X = adata.layers["log1p_transformed"]
adata.var.index = pd.Index(gen.upper() for gen in adata.var.index.values)


adata = adata[:,~adata.var.index.isin(all_gene_symbols)]
print(adata)
"""
print(adata.var.index.values)
print(adata.layers["raw"][:10])"""

adata.uns['log1p']["base"] = None

to_be_subcluster = cluster_groups.split(",")
print(to_be_subcluster)
print(adata.obs[cluster_column].cat.categories)
for ind, cluster_number in enumerate(to_be_subcluster):
    if len(resolution_param) == 1:  
        sc.tl.leiden(adata, restrict_to=(cluster_column, [str(cluster_number)]),  resolution=resolution_param, key_added=cluster_column)
    else:
        sc.tl.leiden(adata, restrict_to=(cluster_column, [str(cluster_number)]),  resolution=resolution_param[ind], key_added=cluster_column)

# cluster rename dict
renamed_clusters = {
    "B cells-1,0":"B cells-1,0", "B cells-1,3":"B cells-1,0", "B cells-1,1":"B cells-1,1", "B cells-1,2":"B cells-1,2",
    'B cells-2,0':'B cells-2,0', 'B cells-2,1':'B cells-2,0', 'B cells-2,2':'B cells-2,0', 'B cells-2,4':'B cells-2,1', 'B cells-2,3':'B cells-2,2', 'B cells-2,5':'B cells-2,3',
    'B cells-3':'B cells-3', 'B cells-4':'B cells-4'}


adata.obs["cell_type_subclustered"] = [renamed_clusters[clust] for clust in adata.obs[cluster_column]]
# adata = adata[adata.obs["cell type subclustered"]!="B cells-2,3",:]
# change the 
cluster_column = "x"
adata.obs[cluster_column]=adata.obs[cluster_column].astype('str').astype('category')

sc.tl.dendrogram(adata, groupby=cluster_column)

sc.tl.rank_genes_groups(adata, groupby=cluster_column, method='wilcoxon', key_added = f"wilcoxon_{cluster_column}")
# mpl.rcParams['axes.titlesize'] = 20
# mpl.rcParams['font.size'] = 40
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key=f"wilcoxon_{cluster_column}", show=False, groupby=cluster_column, save=f'{sample_type}_one_vs_rest_subcluster')#
# mpl.rcParams['axes.titlesize'] = 60
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=f"wilcoxon_{cluster_column}", show=False, groupby=cluster_column, save=f'{sample_type}_subcluster_5_genes')
sc.pl.rank_genes_groups_dotplot(adata, n_genes=15, key=f"wilcoxon_{cluster_column}", show=False, groupby=cluster_column, save=f'{sample_type}_subcluster_15_genes')

plt.rcParams['figure.dpi']= 150
plt.rcParams['figure.figsize']= (30, 20)
sc.pl.umap(adata, color=cluster_column, palette=sc.pl.palettes.default_20, size=45, show=False, legend_loc='on data', save = f'{sample_type}_subclustered_on_data')
sc.pl.umap(adata, color=cluster_column, palette=sc.pl.palettes.default_20, size=45, show=False, save = f'{sample_type}_subclustered')

# Merge the clusters
"""
adata.obs[cluster_column][adata.obs[cluster_column].isin(['2', '3'])]='0'
adata.obs[cluster_column]=adata.obs[cluster_column].astype('str').astype('category')
### Reorder and rename the Leiden
adata.obs['Leiden'].cat.rename_categories(np.arange(len(np.unique(adata.obs['Leiden']))).astype('str'), inplace=True)

"""
for cat in adata.obs[cluster_column].cat.categories:
    print(cat, adata[adata.obs[cluster_column]==cat].shape)
"""





"""
if save_adata:
    print(f"Saving the object... {sample_type}_integrated_subclustered.h5ad...")
    # Write to file
    adata.write(os.path.join(output_path, f'{sample_type}_integrated_subclustered.h5ad'))

# python sc_subcluster_bcells.py -i ../data/out_data/atlas_bcell_populations_corrected.h5ad -o ../data/out_data/ -st atlas_bcells_subclustered  -an atlas_bcells_subcluster -cc cell_type -res 0.2,0.1 -cg "B cells-1,B cells-2" --save

