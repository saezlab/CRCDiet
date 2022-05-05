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


# Read integrated object
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Run annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']

###############################

sample_type = "sc"

S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
PLOT_PATH =  os.path.join(S_PATH, "../plots", "annotate")

Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = PLOT_PATH

adata = sc.read_h5ad(input_path)

# https://github.com/scverse/scanpy/issues/2239
# if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
# KeyError: 'base'
adata.uns['log1p']["base"] = None

l_param, _ = adata.uns["leiden_best_silh_param"]

sample_type = "sc"

# print(l_param)
# print(adata)

sc.tl.rank_genes_groups(adata, groupby=f"leiden_{l_param}", method='wilcoxon', key_added = f"wilcoxon_{l_param}")

sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=f"wilcoxon_{l_param}", show=False, groupby=f"leiden_{l_param}", save=f'{sample_type}_deg_clusters_dotplot_{l_param}')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key=f"wilcoxon_{l_param}", show=False, groupby=f"leiden_{l_param}", save=f'{sample_type}_one_vs_rest_{l_param}')

wc = sc.get.rank_genes_groups_df(adata, group=None, key=f"wilcoxon_{l_param}", pval_cutoff=0.01, log2fc_min=0)[["group", "names", "scores","logfoldchanges"]]
print(l_param)
print(wc.to_csv(os.path.join(output_path, f'{sample_type}_deg_leiden_res_{l_param}.csv'), index=False))



sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon", show=False, groupby="leiden", show_gene_labels=True, save=f'{sample_type}_deg_clusters_heatmap')

sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=f"wilcoxon_{l_param}", show=False, groupby=f"leiden_{l_param}", save=f'{sample_type}_deg_clusters_dotplot_{l_param}')
"""
marker_genes = None
own_markers = False
marker_db="PanglaoDB"

if own_markers:
    # adata = adata.raw.to_adata()
    # TODO: make this parametric
    df_markers = pd.read_csv("..", sep=";")
    marker_genes = []

    # adata_filtered_mouse.var.index = pd.Index(gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.index.values)
    # print([gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.gene_ids.values])

    for ind, row in df_markers.iterrows():

        for cell_type in df_markers.columns:
            if [row[cell_type], cell_type] not in marker_genes:
                marker_genes.append([row[cell_type], cell_type])

    marker_genes = pd.DataFrame(marker_genes, columns=['genesymbol', 'cell_type']) 
    marker_genes['weight'] = 1
else:
    
    # print(marker_genes)
    marker_genes = dc.get_resource('PanglaoDB')
    marker_db="PanglaoDB"
    # Filter by canonical_marker and human
    print(marker_genes)
    marker_genes = marker_genes[(marker_genes["mouse"]=='True')&(marker_genes['canonical_marker']=='True')]

    # Remove duplicated entries
    marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]

print(marker_genes)

#marker_genes['genesymbol'] = marker_genes['genesymbol'].str.upper()

dc.run_ora(mat=adata, net=marker_genes, source='cell_type', target='genesymbol', min_n=3)

# print(adata)
# dc.run_ora(mat=adata, net=markers, source='cell_type', target='genesymbol', min_n=3, verbose=True)
acts = dc.get_acts(adata, obsm_key='ora_estimate')
"""
"""print(adata.obsm)
print("============= ORA p-Vals ============= ")
print(adata.obsm["ora_pvals"])
print("============= ora_estimate ============= ")
print(adata.obsm["ora_estimate"])"""
"""
dict_mean_enr = dict()



mean_enr = dc.summarize_acts(acts, groupby=f'leiden_{l_param}')

sns.clustermap(mean_enr, xticklabels=mean_enr.columns)
plt.savefig(f"{PLOT_PATH}/{sample_type}_clustermap_res_{l_param}.pdf" if own_markers else f"{PLOT_PATH}/{sample_type}_clustermap_{marker_db}_res_{l_param}.pdf") 


annotation_dict = dc.assign_groups(mean_enr)
# annotation_dict06 = dc.assign_groups(mean_enr06)
# annotation_dict04 = dc.assign_groups(mean_enr04)

if own_markers:
    # Add cell type column based on annotation
    adata.obs[f'cell_type_{l_param}'] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{l_param}']]
else:
    # Add cell type column based on annotation
    adata.obs[f'cell_type_{l_param}_panglao'] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{l_param}']]
        

if own_markers:
    # Visualize
    sc.pl.umap(adata, color=f'cell_type_{l_param}', show=True, legend_loc="on data", legend_fontsize="xx-small", save=f'{sample_type}_cell_type_res.pdf' if own_markers else f'{sample_type}_cell_type_PanglaoDB_res.pdf')
else:
    sc.pl.umap(adata, color=f'cell_type_{l_param}_panglao', show=True, legend_loc="on data", legend_fontsize="xx-small", save=f'{sample_type}_cell_type_res.pdf' if own_markers else f'{sample_type}_cell_type_PanglaoDB_res.pdf')

"""


# python cluster_annotate.py -i ../data/out_data/mouse_integrated.h5ad -o ../data/out_data -st mouse
# python cluster_annotate.py -i ../data/out_data/human_integrated.h5ad -o ../data/out_data -st human
