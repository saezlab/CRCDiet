import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import warnings
import utils
import os

'''
Open all samples QC processed files, merge, perform HVGs selection and save the AnnData object
'''

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", fontsize=8, dpi=150, dpi_save=300)
plt.rcParams['figure.constrained_layout.use'] = True
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run Merging')
parser.add_argument('-i', '--input_dir', help='Input directory containing the preprocessed AnnData object ', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the processed object', required=True)
parser.add_argument('-n', '--normalization', default="log1p", help='Normalization technique', required=False)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)
parser.add_argument('-an', '--analysis_name', default="an", help='Analysis name', required=False)

args = vars(parser.parse_args())
input_path = args['input_dir']
output_path = args['output_dir']
normalization = args['normalization']
sample_type = args['sample_type']
analysis_name = args['analysis_name']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ###############################

# sample_type = "sc"

# Load meta data
meta = utils.get_meta_data(sample_type)
samples = np.unique(meta['sample_id'])

markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
markers = list(set(markers_df["genesymbol"].str.capitalize()))

adata = utils.get_filtered_concat_data(sample_type)
print("Shape of merged object:", adata.shape)

print("Calculating QC metrics for merged samples...")
sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],
            wspace=0.3, jitter=0.4, size=0.5, groupby="condition", rotation=90, show=True, save=f"QC_on_merged_objects_after_filtering_{sample_type}_violin.pdf")

print("Normalizing merged samples...")
# keep raw counts in layers
adata.layers['counts'] = adata.X.copy()
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
    # throw NotImplementedError 
"""

print("Calculating HVGs...")
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
# plt.figtext(0.5, 1.0, txt, wrap=True, horizontalalignment='left', fontsize=12)
plt.savefig(os.path.join(os.path.join(PLOT_PATH, f"HVG_vs_unselected_{sample_type}.pdf")))
plt.show();

adata.var = adata.var[['highly_variable','highly_variable_nbatches']]

# Filter by HVG
num_hvg_genes = 4000
# for sc
if sample_type=="sc":
    num_hvg_genes = 4000
# for atlas
elif sample_type=="atlas":
    num_hvg_genes = 5000

batch_msk = np.array(adata.var.highly_variable_nbatches > 1)
hvg = adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata.var['highly_variable'] = [g in hvg for g in adata.var.index]
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata = adata[:,hvg]

print("Performing analytic Pearson residual normalization...")
sc.experimental.pp.normalize_pearson_residuals(adata)

# Run PCA
# sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack', random_state=0)

# Get loadings for each gene for each PC
df_loadings = pd.DataFrame(adata.varm['PCs'], index=adata.var_names)
# get rank of each loading for each PC
df_rankings = pd.DataFrame((-1 * df_loadings.values).argsort(0).argsort(0), index=df_loadings.index, columns=df_loadings.columns)

# Run UMAP to see the difference after integration
print("Computing neighbors...")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
print("\nUMAP of merged objects before integration")
plt.rcParams['figure.dpi']= 300
plt.rcParams['figure.figsize']= (45, 30)
sc.pl.umap(adata, color=["condition"], size=10, palette=sc.pl.palettes.default_20, save=f'{sample_type}_merged_condition.pdf');

# plt.clf()
print("Saving the merged object...")
# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_merged.h5ad'))


# python sc_merge.py -i ../data/out_data -o ../data/out_data