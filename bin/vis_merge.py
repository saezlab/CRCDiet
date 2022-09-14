# import imp
from string import whitespace
from tkinter import font
from turtle import width
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

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", fontsize=8, dpi=80, dpi_save=300)
plt.rcParams['figure.constrained_layout.use'] = True
# Read command line and set args
parser = argparse.ArgumentParser(prog='merge', description='Run Merging')
parser.add_argument('-i', '--input_dir', help='Input directory containing the preprocessed AnnData object ', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the processed object', required=True)
parser.add_argument('-n', '--normalization', default="log1p", help='Normalization technique', required=False)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
args = vars(parser.parse_args())
input_path = args['input_dir']
output_path = args['output_dir']
normalization = args['normalization']
analysis_name = args['analysis_name'] # "visium_merge"
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ###############################

sample_type = "visium"
# Load meta data
meta = utils.get_meta_data(sample_type)
samples = np.unique(meta['sample_id'])

# put the samples in a list
# mpl.rcParams['figure.dpi']= 150

markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
markers = list(set(markers_df["genesymbol"].str.capitalize()))

adatas = []

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
sc.pp.calculate_qc_metrics(adata, inplace=True)


"""
fig, axes = plt.subplots(1, 6, figsize=(12, 6))
for ax, sample in zip(axes, samples):

    adata = sc.read_h5ad(os.path.join(input_path,f"{sample}_filtered.h5ad"))
    # adata.var.index = pd.Index(gen.capitalize() for gen in tmp.var.index.values)
    sc.experimental.pp.highly_variable_genes(
        adata, flavor="pearson_residuals", n_top_genes=4000
    )
    
    # Fetch sample metadata
    m = meta[meta['sample_id'] == sample]    
    # Add metadata to adata
    for col in m.columns:
        adata.obs[col] = m[col].values[0]
    
    hvgs = adata.var["highly_variable"]
    # print(hvgs)
    ax.scatter(
        adata.var["mean_counts"], adata.var["residual_variances"], s=3, edgecolor="none"
    )
    ax.scatter(
        adata.var["mean_counts"][hvgs],
        adata.var["residual_variances"][hvgs],
        c="tab:red",
        label="selected genes",
        s=3,
        edgecolor="none",
    )
    ax.scatter(
        adata.var["mean_counts"][np.isin(adata.var_names, markers)],
        adata.var["residual_variances"][np.isin(adata.var_names, markers)],
        c="k",
        label="known marker genes",
        s=10,
        edgecolor="none",
    )
    ax.set_xscale("log")
    ax.set_xlabel("mean expression")
    ax.set_yscale("log")
    ax.set_ylabel("residual variance")
    # print(adata.obs["condition"][0])
    # ax.set_title(adata.uns["name"])
    ax.set_title(adata.obs["condition"][0])

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    adatas.append(adata)
    del adata

plt.legend()


for adata in adatas:
    adata = adata[:, adata.var["highly_variable"]]
    sc.experimental.pp.normalize_pearson_residuals(adata)


for adata in adatas:
    # Run PCA
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack', random_state=0)
    n_cells = len(adata)
    

for adata in adatas:
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

for adata in adatas:
    print(adata.obs["condition"][0], ":")
    sc.pl.umap(adata, color=["leiden"], cmap="tab20")
    sc.pl.umap(adata, color=list(set(adata.var.index) & set(markers)), layer="sqrt_norm")
"""


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],
            wspace=0.3, jitter=0.4, size=0.5, groupby="condition", rotation=75, show=True, save=f"QC_on_merged_objects_after_filtering_{sample_type}_violin.pdf")

# adata.obs['outlier_total'] = adata.obs.total_counts > 30000
# print('%u cells with large total counts' % (sum(adata.obs['outlier_total'])))
# adata_pbmc3k = adata[~adata.obs['outlier_total'], :]

# keep raw counts in layers
adata.layers['counts'] = adata.X.copy()
adata.layers["sqrt_norm"] = np.sqrt(
    sc.pp.normalize_total(adata, inplace=False)["X"]).copy()

adata.layers["log1p_transformed"] = sc.pp.normalize_total(adata, inplace=False, target_sum=1e6)["X"]
sc.pp.log1p(adata, layer="log1p_transformed")



"""
# TODO: Try other normalization techniques
if normalization == "log1p":

    # Log-normalize expression
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    adata.layers['normalized'] = adata.X
elif normalization == "pearson":
    sc.experimental.pp.highly_variable_genes(
        adata, batch_key='batch', flavor="pearson_residuals", n_top_genes=3000
    )
else:
    # TODO: throw an error
    pass
    # throw NotImplementedError """


# 4000 di
# adata.var.index = pd.Index(gen.capitalize() for gen in tmp.var.index.values)
sc.experimental.pp.highly_variable_genes(
    adata, batch_key='batch', flavor="pearson_residuals", n_top_genes=3000
)



fig, ax = plt.subplots(1, 1, figsize=(12, 6))
hvgs = adata.var["highly_variable"]

ax.scatter(
    adata.var["mean_counts"], adata.var["residual_variances"], s=3, edgecolor="none"
)
ax.scatter(
    adata.var["mean_counts"][hvgs],
    adata.var["residual_variances"][hvgs],
    c="tab:red",
    label="selected genes",
    s=3,
    edgecolor="none",
)
ax.scatter(
    adata.var["mean_counts"][np.isin(adata.var_names, markers)],
    adata.var["residual_variances"][np.isin(adata.var_names, markers)],
    c="k",
    label="known marker genes",
    s=10,
    edgecolor="none",
)
ax.set_xscale("log")
ax.set_xlabel("mean expression")
ax.set_yscale("log")
ax.set_ylabel("residual variance")
# ax.set_title(adata.obs["condition"][0])

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")
plt.legend()

txt="In the figure below, red dots show the selected genes (i.e. highly variable genes) and the black ones represent the marker genes."
plt.figtext(0.5, 1.0, txt, wrap=True, horizontalalignment='left', fontsize=12)
plt.show();
# Compute HVG
# sc.experimental.pp.highly_variable_genes(adata, batch_key='batch', n_top_genes=4000)
# sc.pl.highly_variable_genes(adata, save=f'{sample_type}_merged_hvg.pdf')
# plt.show()

adata.var = adata.var[['highly_variable','highly_variable_nbatches']]

# Filter by HVG
num_hvg_genes = 5000
batch_msk = np.array(adata.var.highly_variable_nbatches > 1)
hvg = adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata.var['highly_variable'] = [g in hvg for g in adata.var.index]
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata = adata[:,hvg]

print("Performing analytic Pearson residual normalization...")
sc.experimental.pp.normalize_pearson_residuals(adata)

# print(adata.X)

# Run PCA
# sc.pp.scale(adata)

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

print("Computing neighbors...")
# Run UMAP to see the difference after integration
sc.pp.neighbors(adata)
print("\nUMAP of merged objects before integration")
sc.tl.umap(adata)
sc.pl.umap(adata, color=["condition"], palette=sc.pl.palettes.default_20, save=f'{sample_type}_merged_condition.pdf');
plt.show();

# plt.clf()
print("Saving the merged object...")
# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_merged.h5ad'))


# python vis_merge.py -i ../data/out_data -o ../data/out_data