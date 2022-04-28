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


# Read integrated object
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run Clustering and annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-st', '--sample_type', help='Human, mouse or tumor', required=True)
parser.add_argument('-res', '--resolution',  help='Resolution param for leiden', type=float, default=1.0)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']
sample_type = args['sample_type']

###############################

plot_path="../plots/cluster_annotate"
Path(plot_path).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = plot_path

adata = sc.read_h5ad(input_path)

leiden_res_params = [0.1, 0.2, 0.5, 0.7, 1.0]


sc.pp.neighbors(adata)

# perform clustering, Rank genes for characterizing groups, plot top 5 genes
for l_param in leiden_res_params:
    sc.tl.leiden(adata, resolution = l_param, key_added = f"leiden_{l_param}") # default resolution in 1.0
    print(type(adata.obs[f"leiden_{l_param}"]))
    print(list(adata.obs[f"leiden_{l_param}"].cat.categories))
    
    sc.tl.rank_genes_groups(adata, groupby=f"leiden_{l_param}", method='wilcoxon', key_added = f"wilcoxon_{l_param}")
    
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key=f"wilcoxon_{l_param}", show=False, groupby=f"leiden_{l_param}", save=f'{sample_type}_deg_clusters_dotplot_{l_param}')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key=f"wilcoxon_{l_param}", show=False, groupby=f"leiden_{l_param}", save=f'{sample_type}_one_vs_rest_{l_param}')

    wc = sc.get.rank_genes_groups_df(adata, group=None, key=f"wilcoxon_{l_param}", pval_cutoff=0.01, log2fc_min=0)[["group", "names", "scores","logfoldchanges"]]
    print(l_param)
    print(wc.to_csv(os.path.join(output_path, f'{sample_type}_deg_leiden_res_{l_param}.csv'), index=False))






"""
sc.tl.rank_genes_groups(adata, 'leiden_1.0', method='wilcoxon', key_added = "wilcoxon_1.0")
sc.tl.rank_genes_groups(adata, 'leiden_0.6', method='wilcoxon', key_added = "wilcoxon_0.6")
sc.tl.rank_genes_groups(adata, 'leiden_0.4', method='wilcoxon', key_added = "wilcoxon_0.4")
# sc.pl.rank_genes_groups(adata)

# sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon", show=False, groupby="leiden", show_gene_labels=True, save=f'{sample_type}_deg_clusters_heatmap')

sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon_1.0", show=False, groupby="leiden_1.0", save=f'{sample_type}_deg_clusters_dotplot_1.0')
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon_0.6", show=False, groupby="leiden_0.6", save=f'{sample_type}_deg_clusters_dotplot_0.6')
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon_0.4", show=False, groupby="leiden_0.4", save=f'{sample_type}_deg_clusters_dotplot_0.4')
"""

# sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=5, key="wilcoxon", show=False, groupby="leiden", save=f'{sample_type}_ranked_genes_stackedviolin')


sc.pl.umap(adata, color=[f"leiden_{l_param}" for l_param in leiden_res_params], legend_loc="on data", legend_fontsize="xx-small", palette=sc.pl.palettes.default_20, show=False, save=f'{sample_type}_clusters_diff_res')
"""sc.pl.umap(
    adata, color=["leiden"], palette=sc.pl.palettes.default_20, show=False, save=f'{sample_type}_clusters'
)"""

sc.pl.umap(
    adata, color=["sample_id"], palette=sc.pl.palettes.default_20, show=False, save=f'{sample_type}_samples'
)

marker_genes = None
own_markers = True
marker_db="PanglaoDB"

if own_markers:
    # adata = adata.raw.to_adata()
    df_markers = pd.read_csv("../data/Stroma_oldclusters_DEgenes_Venice_top100xlsx.csv", sep=";")
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
    marker_genes = marker_genes[(marker_genes[sample_type]=='True')&(marker_genes['canonical_marker']=='True')]

    # Remove duplicated entries
    marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]

print(marker_genes)

#marker_genes['genesymbol'] = marker_genes['genesymbol'].str.upper()

dc.run_ora(mat=adata, net=marker_genes, source='cell_type', target='genesymbol', min_n=3)

# print(adata)
# dc.run_ora(mat=adata, net=markers, source='cell_type', target='genesymbol', min_n=3, verbose=True)
acts = dc.get_acts(adata, obsm_key='ora_estimate')

print(adata.obsm)
print("============= ORA p-Vals ============= ")
print(adata.obsm["ora_pvals"])
print("============= ora_estimate ============= ")
print(adata.obsm["ora_estimate"])

dict_mean_enr = dict()

for l_param in leiden_res_params:

    mean_enr = dc.summarize_acts(acts, groupby=f'leiden_{l_param}')
    
    sns.clustermap(mean_enr, xticklabels=mean_enr.columns)
    plt.savefig(f"{plot_path}/{sample_type}_clustermap_res_{l_param}.pdf" if own_markers else f"{plot_path}/{sample_type}_clustermap_{marker_db}_res_{l_param}.pdf") 
    """sns.clustermap(mean_enr06, xticklabels=mean_enr06.columns)
    plt.savefig(f"{plot_path}/{sample_type}_clustermap_res06.pdf" if own_markers else f"{plot_path}/{sample_type}_clustermap_{marker_db}_res06.pdf") 
    sns.clustermap(mean_enr04, xticklabels=mean_enr04.columns)
    plt.savefig(f"{plot_path}/{sample_type}_clustermap_res04.pdf" if own_markers else f"{plot_path}/{sample_type}_clustermap_{marker_db}_res04.pdf")"""

    annotation_dict = dc.assign_groups(mean_enr)
    """annotation_dict06 = dc.assign_groups(mean_enr06)
    annotation_dict04 = dc.assign_groups(mean_enr04)"""

    if own_markers:
        # Add cell type column based on annotation
        adata.obs[f'cell_type_{l_param}'] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{l_param}']]
    else:
        # Add cell type column based on annotation
        adata.obs[f'cell_type_{l_param}_panglao'] = [annotation_dict[clust] for clust in adata.obs[f'leiden_{l_param}']]
        

if own_markers:
    # Visualize
    sc.pl.umap(adata, color=[f'cell_type_{l_param}' for l_param in leiden_res_params], show=False, legend_loc="on data", legend_fontsize="xx-small", save=f'{sample_type}_cell_type_res.pdf' if own_markers else f'{sample_type}_cell_type_PanglaoDB_res.pdf')
else:
    sc.pl.umap(adata, color=[f'cell_type_{l_param}_panglao'for l_param in leiden_res_params], show=False, legend_loc="on data", legend_fontsize="xx-small", save=f'{sample_type}_cell_type_res.pdf' if own_markers else f'{sample_type}_cell_type_PanglaoDB_res.pdf')




# sc.pl.umap(tmp, color=list(tmp.var.index)+['leiden']+['sample_id'],  show=False, save=f"{sample_type}_cell_type")

# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_integrated.h5ad'))


# python cluster_annotate.py -i ../data/out_data/mouse_integrated.h5ad -o ../data/out_data -st mouse
# python cluster_annotate.py -i ../data/out_data/human_integrated.h5ad -o ../data/out_data -st human
