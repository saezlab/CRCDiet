# import imp
import scanpy as sc
import scanpy.external as sce
from utils import get_meta_data
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import utils
import os
import argparse
import warnings

warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0

'''
Open all samples QC processed files, merge them
'''

S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
PLOT_PATH =  os.path.join(S_PATH, "../plots", "merge")

Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = PLOT_PATH

sc.set_figure_params(scanpy=True, dpi=150, dpi_save=300)

# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run Merging')
parser.add_argument('-i', '--input_dir', help='Input directory containing the preprocessed AnnData object ', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the processed object', required=True)

args = vars(parser.parse_args())

input_path = args['input_dir']
output_path = args['output_dir']
###############################

sample_type = "sc"
# Load meta data
meta = utils.get_meta_data("sc")
samples = np.unique(meta['sample_id'])

# put the samples in a list
adata = []
# for sample in os.listdir(input_path):
for sample in samples:

    tmp = sc.read_h5ad(os.path.join(input_path,f"{sample}_filtered.h5ad"))

    # Fetch sample metadata
    m = meta[meta['sample_id'] == sample]
    
    # Add metadata to adata
    for col in m.columns:
        tmp.obs[col] = m[col].values[0]

    # Append
    adata.append(tmp)
    del tmp
    
# Merge objects and delete list
adata = adata[0].concatenate(adata[1:], join='outer')


# Log-normalize expression
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

# Compute HVG
sc.pp.highly_variable_genes(adata, batch_key='batch')
sc.pl.highly_variable_genes(adata, save=f'{sample_type}_merged_hvg.pdf')
plt.show()


adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata.raw = adata


# Filter by HVG
num_hvg_genes = 3000
batch_msk = np.array(adata.var.highly_variable_nbatches > 1)
hvg = adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata.var['highly_variable'] = [g in hvg for g in adata.var.index]
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata = adata[:,hvg]

# Run PCA
sc.pp.scale(adata)

sc.tl.pca(adata, svd_solver='arpack', random_state=0)

# Get loadings for each gene for each PC
df_loadings = pd.DataFrame(adata.varm['PCs'], index=adata.var_names)
# get rank of each loading for each PC
df_rankings = pd.DataFrame((-1 * df_loadings.values).argsort(0).argsort(0), index=df_loadings.index, columns=df_loadings.columns)
# c.f. with df_loadings.apply(scipy.stats.rankdata, axis=0)
# evaluate 
# print("Top loadings for PC1...")
# print(df_loadings[0].sort_values().tail())
# print("Rank of IKZF1 for first 5 PCs...")
# print(df_rankings.loc["IKZF1"].head())

"""sc.pl.pca_overview(adata, color='sample_id', show=False, save=f'{sample_type}_pca_overview.pdf')

sc.pl.pca_loadings(adata, components=[1,2,3,4,5,6,7,8],  show=False, save=f'{sample_type}_pca_loadings.pdf')

sc.pl.pca_variance_ratio(adata, n_pcs = 50,  show=False, save=f'{sample_type}_variance_ratio.pdf')"""

# Run UMAP to see the difference after integration
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["condition"], palette=sc.pl.palettes.default_20, save=f'{sample_type}_merged_condition.pdf');
plt.show();
plt.clf()

# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_merged.h5ad'))


# python merge.py -i ../data/out_data -o ../data/out_data  -st mouse
# python merge.py -i ../data/out_data -o ../data/out_data  -st human