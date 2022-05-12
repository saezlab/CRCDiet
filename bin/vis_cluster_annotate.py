from genericpath import sameopenfile
from operator import index
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns 
import decoupler as dc  
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import os
import sys
from sklearn.metrics import silhouette_score, pairwise_distances
import sys
import warnings
from utils import printmd
import anndata
import utils
import matplotlib as mpl

sc.settings.verbosity = 0

warnings.simplefilter(action='ignore')

# Read integrated object
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Run annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']

###############################

sample_type = "visium"
meta = utils.get_meta_data("visium")

S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
PLOT_PATH =  os.path.join(S_PATH, "../plots", "visium_annotate")

Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = PLOT_PATH
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80)

adata = sc.read_h5ad(input_path)

adata.raw = anndata.AnnData(adata.layers['counts'], obs=adata.obs, var=adata.var)
# print(adata.layers['counts'])
# https://github.com/scverse/scanpy/issues/2239
# if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
# KeyError: 'base'
# adata.uns['log1p']["base"] = None

l_param, _ = adata.uns["leiden_best_silh_param"]
l_param = f"{l_param:.2f}"

# print(l_param)
# print(adata)

sc.tl.rank_genes_groups(adata, groupby=f"leiden_{l_param}", method='wilcoxon', key_added = f"wilcoxon_{l_param}")
mpl.rcParams['axes.titlesize'] = 20
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key=f"wilcoxon_{l_param}", show=True, groupby=f"leiden_{l_param}", save=f'{sample_type}_one_vs_rest_{l_param}')#
# mpl.rcParams['axes.titlesize'] = 60
# sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=f"wilcoxon_{l_param}", show=True, groupby=f"leiden_{l_param}", save=f'{sample_type}_deg_clusters_dotplot_{l_param}')

wc = sc.get.rank_genes_groups_df(adata, group=None, key=f"wilcoxon_{l_param}", pval_cutoff=0.01, log2fc_min=0)[["group", "names", "scores","logfoldchanges"]]
# print(l_param)
# print(wc.to_csv(os.path.join(output_path, f'{sample_type}_deg_leiden_res_{l_param}.csv'), index=False))

marker_genes = None
own_markers = True
marker_db="PanglaoDB"

if own_markers:
    # adata = adata.raw.to_adata()
    # TODO: make this parametric
    marker_genes = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")

    """# adata_filtered_mouse.var.index = pd.Index(gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.index.values)
    # print([gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.gene_ids.values])

    for ind, row in df_markers.iterrows():

        for cell_type in df_markers.columns:
            if [row[cell_type], cell_type] not in marker_genes:
                marker_genes.append([row[cell_type], cell_type])

    marker_genes = pd.DataFrame(marker_genes, columns=['genesymbol', 'cell_type']) """
    marker_genes["genesymbol"] = marker_genes["genesymbol"].str.capitalize()
    marker_genes['weight'] = 1
    marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]
    
else:
    
    # print(marker_genes)
    # marker_genes = dc.get_resource('PanglaoDB')
    # marker_db="PanglaoDB"
    marker_genes = pd.read_csv(os.path.join(DATA_PATH,"PanglaoDB_markers_27_Mar_2020.tsv"), sep="\t")
    print(marker_genes)
    # Filter by canonical_marker and human
    
    marker_genes = marker_genes[marker_genes["species_official"].str.contains("mouse") & marker_genes["canonical_marker"]==1.0]
    
    print(marker_genes["genesymbol"])
    # Remove duplicated entries
    marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]

dc.run_ora(mat=adata, net=marker_genes, source='cell_type', target='genesymbol', min_n=3, verbose=True)
acts = dc.get_acts(adata, obsm_key='ora_estimate')

"""print("============= ORA p-Vals ============= ")
print(adata.obsm["ora_pvals"])
print("============= ora_estimate ============= ")
print(adata.obsm["ora_estimate"])"""
for c_t in list(adata.obsm["ora_estimate"].columns):
    adata.obs[c_t] = adata.obsm["ora_estimate"][c_t]

for c_t in list(adata.obsm["ora_estimate"].columns):
    rows, cols = (1, 6)
    fig, ax = plt.subplots(rows, cols, figsize=(20,2))

    for ind, row in meta.iterrows():
        
        # fig_row, fig_col = int(ind/cols), ind%cols
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata_raw = utils.read_raw_visium_sample(sample_id)
        adata_temp = adata[adata.obs["condition"]==condition,:]    
        adata_temp.obs.index = pd.Index("-".join(cl.split("-")[:-1]) for cl in adata_temp.obs.index.values)
        adata_raw = adata_raw[adata_temp.obs.index,:]
        adata_raw = adata_raw[:, list(adata.var_names)]
        adata_raw.obs[c_t] = adata_temp.obs[c_t]
        adata_raw.obs[f"leiden_{l_param}"] = adata_temp.obs[f"leiden_{l_param}"]
        mpl.rcParams["image.cmap"]= plt.cm.magma
        mpl.rcParams['axes.titlesize'] = 14
        sc.pl.spatial(adata_raw, img_key="hires", color =c_t, title=f"{c_t}_{condition}", size=1.25, alpha_img=0.5, ax = ax[ind], show=False) # ax[fig_row][fig_col]
        cbar = ax[ind].collections[0].colorbar
        cbar.set_ticks([])
        # plt.tight_layout()
    
    
        # ax[fig_row][fig_col].colorbar(fraction=0.046, pad=0.04)
        # plt.tight_layout(pad=3.0)

        #sc.pl.violin(adata, list(set(adata.var.index) & set(markers)), show=True, groupby=f"leiden_{l_param}")

        
        #mpl.rcParams["image.cmap"]= plt.cm.magma_r
        #sc.pl.spatial(adata_raw, img_key="hires", color =list(set(adata.var.index) & set(markers)),  size=1.0, alpha_img=0.5, wspace = 0.3)
    # plt.tight_layout()
    plt.show();

mpl.rcParams['axes.titlesize'] = 20
cell_types = list(set(marker_genes["cell_type"]))
sc.pl.umap(adata, color=list(adata.obsm["ora_pvals"].columns))

dict_mean_enr = dict()

mean_enr = dc.summarize_acts(acts, groupby=f'leiden_{l_param}')

"""sns.clustermap(mean_enr, xticklabels=mean_enr.columns)
plt.savefig(f"{PLOT_PATH}/{sample_type}_clustermap_res_{l_param}.pdf" if own_markers else f"{PLOT_PATH}/{sample_type}_clustermap_{marker_db}_res_{l_param}.pdf") 


annotation_dict = dc.assign_groups(mean_enr)
print(annotation_dict)
# annotation_dict06 = dc.assign_groups(mean_enr06)
# annotation_dict04 = dc.assign_groups(mean_enr04)

if own_markers:
    # Add cell type column based on annotation
    adata.obs[f'cell_type_{l_param}'] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{l_param}']]
else:
    # Add cell type column based on annotation
    adata.obs[f'cell_type_{l_param}_panglao'] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{l_param}']]
"""   

"""if own_markers:
    # Visualize
    sc.pl.umap(adata, color=f'cell_type_{l_param}', show=True, legend_loc="on data", legend_fontsize="xx-small", save=f'{sample_type}_cell_type_res.pdf' if own_markers else f'{sample_type}_cell_type_PanglaoDB_res.pdf')
else:
    sc.pl.umap(adata, color=f'cell_type_{l_param}_panglao', show=True, legend_loc="on data", legend_fontsize="xx-small", save=f'{sample_type}_cell_type_res.pdf' if own_markers else f'{sample_type}_cell_type_PanglaoDB_res.pdf')
"""



# python vis_cluster_annotate.py -i ../data/out_data/visium_integrated.h5ad -o ../data/out_data
