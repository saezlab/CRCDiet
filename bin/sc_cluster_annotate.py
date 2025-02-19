import matplotlib.pyplot as plt
import scanpy as sc
import argparse
from sklearn.metrics import silhouette_score, pairwise_distances
import sys
import warnings
import utils
import os
import pandas as pd
from utils import printmd
import matplotlib as mpl

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Run annotation')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
parser.add_argument('-of', '--output_file', help='Output file name', required=False)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)
parser.add_argument('-oc', '--obs_column', default="oc", help='Obs column', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # "sc_cluster_annotate"
output_file = args['output_file']
sample_type = args['sample_type']
obs_column  = args['obs_column']

# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
###############################



adata = sc.read_h5ad(input_path)
print(adata)
"""print(adata)
for col in adata.obs.columns:
    print(col)
    if col not in ['condition', 'batch', 'sample_id']:
        del adata.obs[col]"""
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
# adata.uns['log1p']["base"] = None
# l_param = adata.uns["leiden_best_silh_param"]
# # l_param = f"{l_param:.2f}"



adata_concat.var.index = pd.Index(gen.upper() for gen in adata_concat.var.index.values)
adata_concat.var_names_make_unique()
markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
markers = list(set(markers_df["genesymbol"].str.upper()))

marker_intersect = list(set(adata_concat.var_names) & set(markers))

print(marker_intersect)
broad_markers =["Ptprc", "EpCam", "Pdgfra", "Pdpn", "Myh11", "Ms4a1", "Cd3e", "Itgax", "CD14", "S100a9", "Cpa3", "Mzb1", "Igha"]
broad_markers =["Ptprc", "EpCam", "Pdgfra", "Myh11", "Ms4a1", "Cd3e", "Itgax", "CD14", "S100a9", "Cpa3", "Mzb1", "Igha"]
broad_markers = [mrk.upper() for mrk in broad_markers]

l_param_list = [0.30] # for major cell types
l_param_list = [0.30] # for major cell types
l_param_list = [0.20] # for sc and epicells
"""if sample_type=="atlas":
    l_param_list = [0.40] # for Atlas data
elif sample_type=="sc":
    l_param_list = [0.20] # for SC data"""


step = 0.10
for l_param in l_param_list:
#for l_param in np.arange(0.1, 1.01, step):
    
    l_param = f"{l_param:.2f}"
    group_by = f"leiden_{l_param}"
    if obs_column:
        group_by=obs_column
    printmd(f"## Clusters with resolution param: {l_param} <a class='anchor' id='seventh-bullet-1'></a>")
    
    if obs_column:
        adata_concat.obs[obs_column] = adata.obs[obs_column]
    else:

        adata_concat.obs[f"leiden_{l_param}"] = adata.obs[f"leiden_{l_param}"]
    adata_concat.obsm["X_umap"] = adata.obsm["X_umap"]
    del adata
    
    
    # do not  del adata for mor e than one value
    # del adata
    # mpl.rcParams['figure.dpi']= 300
    #mpl.rcParams["figure.figsize"] = (40,30)
    #mpl.rcParams["legend.fontsize"]  = 30
    
    mpl.rcParams["legend.loc"]  = "upper right"
    mpl.rcParams['axes.facecolor'] = "white"

    sc.tl.rank_genes_groups(adata_concat, method="wilcoxon", groupby=group_by, show=False, key_added = f"wilcoxon_{obs_column}")
    mpl.rcParams['axes.titlesize'] = 20
    sc.pl.rank_genes_groups(adata_concat, n_genes=25, sharey=False, standard_scale='var', key=f"wilcoxon_{obs_column}", show=False, groupby=group_by, save=f'{sample_type}_one_vs_rest_{obs_column}_25.pdf')
    wc = sc.get.rank_genes_groups_df(adata_concat, group=None, key=f"wilcoxon_{obs_column}", pval_cutoff=0.05)# [["group", "names", "scores","logfoldchanges"]]
    wc.to_csv(os.path.join(PLOT_PATH, f"{analysis_name}_rank_genes_groups_df.csv"))
    
    # sc.pl.rank_genes_groups_dotplot(adata_concat, key=f"wilcoxon_{l_param}", standard_scale='var', show=False, groupby=f"leiden_{l_param}", save=f'{sample_type}_deg_clusters_dotplot_{l_param}_default')
    """sc.pl.rank_genes_groups_dotplot(adata_concat, n_genes=5, key=f"wilcoxon_{obs_column}", standard_scale='var',  show=False, groupby=group_by, save=f'{sample_type}_deg_clusters_dotplot_{obs_column}_default')
    sc.pl.rank_genes_groups_dotplot(adata_concat, n_genes=10, key=f"wilcoxon_{obs_column}", standard_scale='var',  show=False, groupby=group_by, save=f'{sample_type}_deg_clusters_dotplot_{obs_column}_10')
    """
    # standard scale false for b cell populations
    sc.pl.rank_genes_groups_dotplot(adata_concat, n_genes=5, key=f"wilcoxon_{obs_column}", standard_scale='var',  show=False, groupby=group_by, save=f'{sample_type}_deg_clusters_dotplot_{obs_column}_default')
    sc.pl.rank_genes_groups_dotplot(adata_concat, n_genes=10, key=f"wilcoxon_{obs_column}", standard_scale='var',  show=False, groupby=group_by, save=f'{sample_type}_deg_clusters_dotplot_{obs_column}_10')

    sc.pl.rank_genes_groups_dotplot(adata_concat, n_genes=10, key=f"wilcoxon_{obs_column}", standard_scale='var',  show=False, groupby=group_by)
    plt.savefig('destination_path.eps', format='eps')
    #adata_concat = adata_concat[:,marker_intersect]
    # mpl.rcParams["figure.figsize"] = 5,5)
    # sc.pl.dotplot(adata_concat, broad_markers, groupby=f'leiden_{l_param}', swap_axes=True, dendrogram=False,  show=True, save=f'{sample_type}_clusters_broad_markers_{l_param}_dotplot')
    # sc.pl.stacked_violin(adata_concat, broad_markers, groupby=f'leiden_{l_param}', dendrogram=False, save=f'{sample_type}_clusters_broad_markers_{l_param}_stacked_violin')
    # sc.pl.stacked_violin(adata_concat, broad_markers, groupby=f'leiden_{l_param}', dendrogram=False, swap_axes=True,  save=f'{sample_type}_clusters_broad_markers_{l_param}_stacked_violin_axis_swapped')
    # sc.pl.violin(adata_concat2, broad_markers, groupby=f'leiden_{l_param}', save=f'{sample_type}_clusters_broad_markers_{l_param}_violin')

    # sc.pl.dotplot(adata_concat, marker_intersect, groupby=f'leiden_{l_param}', swap_axes=True, dendrogram=True,  show=True, save=f'{sample_type}_clusters_all_marker_{l_param}_dotplot_dendogram')
    #sc.pl.dotplot(adata_concat, marker_intersect, groupby=f'leiden_{l_param}', swap_axes=True, dendrogram=False,  show=True, save=f'{sample_type}_clusters_all_marker_{l_param}_dotplot')
    #sc.pl.stacked_violin(adata_concat, marker_intersect, groupby=f'leiden_{l_param}', dendrogram=False, save=f'{sample_type}_clusters_all_marker_{l_param}_stacked_violin')
    # sc.pl.violin(adata_concat, marker_intersect, groupby=f'leiden_{l_param}', save=f'{sample_type}_clusters_all_marker_{l_param}_violin')

    sc.pl.umap(adata_concat, color=obs_column, palette=sc.pl.palettes.default_20, size=8 , show=False, legend_loc='on data', save=f'{sample_type}_leiden_{obs_column}_ondata')
    sc.pl.umap(adata_concat, color=obs_column, palette=sc.pl.palettes.default_20, size=8 , show=False, save=f'{sample_type}_leiden_{obs_column}_umap')
    plt.savefig('destination_path.eps', format='eps')
    """
    # "Pdgrfra",
    markers_dot_plot = markers_dot_plot = ["Epcam", "Agr2", "Fabp2", "Krt14", "Pdgfra", "Myh11", "Ano1", "Lyve1", "Esam", "Ptprc", "Itgax", "Cd3g", "Mzb1", "Jchain", "Il17rb", "Cpa3", "S100a9", "Mki67"]
    sc.pl.dotplot(adata, markers_dot_plot, groupby=f'leiden_{l_param}', swap_axes=True, dendrogram=True,  show=True, save=f'{sample_type}_clusters_marker_{l_param}_dotplot_dendogram')
    sc.pl.dotplot(adata, markers_dot_plot, groupby=f'leiden_{l_param}', swap_axes=True, dendrogram=False,  show=True, save=f'{sample_type}_clusters_marker_{l_param}_dotplot')
    """
    # change below anndata objects to "anndata" to run on only HVGs
    """mpl.rcParams['figure.dpi']= 300
    mpl.rcParams["figure.figsize"] = (5,5)
    print("DEGs per cluster!")
    sc.tl.rank_genes_groups(adata_concat, groupby=f"leiden_{l_param}", method='wilcoxon', key_added = f"wilcoxon_{l_param}")
    adata.uns[f"wilcoxon_{l_param}"] = adata_concat.uns[f"wilcoxon_{l_param}"]
    mpl.rcParams['axes.titlesize'] = 20
    sc.pl.rank_genes_groups(adata_concat, n_genes=25, sharey=False, key=f"wilcoxon_{l_param}", show=True, groupby=f"leiden_{l_param}", save=f'{sample_type}_one_vs_rest_{l_param}')#
    mpl.rcParams['axes.titlesize'] = 60
    sc.pl.rank_genes_groups_dotplot(adata_concat, n_genes=5, key=f"wilcoxon_{l_param}", show=True, groupby=f"leiden_{l_param}", save=f'{sample_type}_deg_clusters_dotplot_{l_param}')"""


    print(f"Saving the object... {sample_type}_integrated_clustered.h5ad...")
    # Write to file
    # adata.write(os.path.join(output_path, f'{sample_type}_integrated_clustered.h5ad'))

    """
    sc.tl.rank_genes_groups(adata_concat, groupby=f"condition", method='wilcoxon', key_added = f"wilcoxon_condition")
    sc.pl.rank_genes_groups(adata_concat, n_genes=25, sharey=False, key=f"wilcoxon_condition", show=True, groupby="condition", save=f'{sample_type}_one_vs_rest_condition')#
    """

    """mpl.rcParams['figure.dpi']= 300
    mpl.rcParams["figure.figsize"] = (5,5)
    mpl.rcParams['axes.titlesize'] = 15
    # mpl.rcParams["font.size"]  = 50
    adata_concat_immune = adata_concat[adata_concat.obs["condition"].str.contains("Immune"), : ].copy()
    sc.tl.rank_genes_groups(adata_concat_immune, groupby=f"condition", method='wilcoxon', key_added = f"wilcoxon_condition_immune")
    sc.pl.rank_genes_groups(adata_concat_immune, n_genes=25, sharey=False, key=f"wilcoxon_condition_immune", show=True, groupby="condition", save=f'{sample_type}_immune_one_vs_rest_condition')
    


    adata_concat_epi = adata_concat[adata_concat.obs["condition"].str.contains("Epi_"), : ].copy()
    sc.tl.rank_genes_groups(adata_concat_epi, groupby=f"condition", method='wilcoxon', key_added = f"wilcoxon_condition_epithelial")
    sc.pl.rank_genes_groups(adata_concat_epi, n_genes=25, sharey=False, key=f"wilcoxon_condition_epithelial", show=True, groupby="condition", save=f'{sample_type}_epithelial_one_vs_rest_condition')
    """
    # wc = sc.get.rank_genes_groups_df(adata_concat, group=None, key=f"wilcoxon_{l_param}", pval_cutoff=0.01, log2fc_min=0)[["group", "names", "scores","logfoldchanges"]]
    # print(l_param)
    # print(wc.to_csv(os.path.join(output_path, f'{sample_type}_deg_leiden_res_{l_param}.csv'), index=False))


# python sc_cluster_annotate.py -i ../data/out_data/sc_epicells_integrated_clustered.h5ad -o ../data/out_data -st sc_epicells_aom_noaom -an sc_epicells_aom_noaom_cluster_3 -oc leiden_0.20

# python sc_cluster_annotate.py -i ../data/out_data/sc_epicells_aom_noaom_concatenated_celltype_annot.h5ad -o ../data/out_data -st sc_epicells_aom_noaom -an sc_epicells_aom_noaom_cluster -oc cluster 
# python sc_cluster_annotate.py -i ../data/out_data/sc_epicells_aom_noaom_concatenated_celltype_annot.h5ad -o ../data/out_data -st sc_epicells_aom_noaom -an sc_epicells_aom_noaom_cluster -oc cell_type