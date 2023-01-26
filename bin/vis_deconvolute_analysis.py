
import seaborn as sns 
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import argparse
import os
import pickle
from matplotlib import rcParams
import utils
import math
import warnings
import argparse


############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
# sc.set_figure_params(scanpy=True, facecolor="white", dpi=50, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Run deconvolution')
parser.add_argument('-vp', '--vis_path', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True) # vis_deconvolution


args = vars(parser.parse_args())
vis_path = args['vis_path']
analysis_name = args['analysis_name'] # vis_deconvolute

# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 


#############################################
################## Part 2 ################### 
#############################################

ref_run_name = f'{OUT_DATA_PATH}/cell2location_{analysis_name}'
run_name = f'{OUT_DATA_PATH}/cell2location_map_{analysis_name}'

lst_samples = ["ext_L24854_CD-AOM-DSS-colon-d81-visium_filtered_deconv_15_20",
"ext_L24854_CD-no-AOM-DSS-colon-d81-visium_filtered_deconv_15_20",
"ext_L24854_HFD-AOM-DSS-colon-d81-visium_filtered_deconv_15_20",
"ext_L24854_HFD-no-AOM-DSS-colon-d81-visium_filtered_deconv_15_20",
"ext_L24854_LFD-AOM-DSS-colon-d81-visium_filtered_deconv_15_20",
"ext_L24854_LFD-no-AOM-DSS-colon-d81-visium_filtered_deconv_15_20"]

cell_types = None
dict_cluster_to_count = dict()
for sample in lst_samples:
    print(f"Processing sample {sample}")
    adata_vis = sc.read_h5ad(f"../data/out_data/cell2location_map/{sample}.h5ad")
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    if not cell_types:
        cell_types = list(adata_vis.uns['mod']['factor_names'])
    
    assert cell_types == list(adata_vis.uns['mod']['factor_names'])
    
    condition = adata_vis.obs["condition"].cat.categories[0]
    adata_vis_integ = sc.read_h5ad("../data/out_data/visium_integrated_clustered.h5ad")
    adata_vis_integ = adata_vis_integ[adata_vis_integ.obs["condition"]==condition]

    adata_vis.obs["leiden_0.10"] = adata_vis_integ.obs["leiden_0.10"].values
    
    for clust_id in  adata_vis.obs["leiden_0.10"].cat.categories:
        
        if clust_id not in dict_cluster_to_count:
            dict_cluster_to_count[clust_id] = []

    for cat in adata_vis.obs["leiden_0.10"].cat.categories:
        adata_temp = adata_vis[adata_vis.obs["leiden_0.10"]==cat]
        # all cell type per spot will be accumulated here
        lst_cell_type_prop_clust = []
        for ind, val in enumerate(list(adata_temp.obs_names)):
            # print(ind, val)
            cell_abundance_list = adata_temp[ind,:].obs[cell_types].values[0]
            norm = [float(i)/max(cell_abundance_list) for i in cell_abundance_list]
            cell_type_prop = [val/sum(norm) for val in norm]
            # append raw for cell abundances
            # append cell_type_prop for proportions
            lst_cell_type_prop_clust.append(cell_abundance_list)
        
        df = pd.DataFrame(lst_cell_type_prop_clust, index=list(adata_temp.obs_names), columns=cell_types)
        # calculate log2(count+1)
        row_sum = np.log2(np.array(df.sum())+1)
        row_sum_norm = [float(i)/max(row_sum) for i in row_sum]
        overall_cell_type_prop = [val/sum(row_sum_norm) for val in row_sum_norm]
        # all_cluster_vs_cell_type_prop.append(overall_cell_type_prop)
        # all_cluster_vs_cell_type_prop.append(list(row_sum))
        dict_cluster_to_count[cat].append(list(row_sum))
    del adata_vis

cluster_numbers = sorted(dict_cluster_to_count.keys())
all_cluster_vs_cell_type_prop  = []
print("cluster categories:", cluster_numbers)
for clust in cluster_numbers:
    all_cluster_vs_cell_type_prop.append(np.sum(dict_cluster_to_count[clust], axis=0))


correlation_matrix  = np.corrcoef(all_cluster_vs_cell_type_prop) 
df_overall = pd.DataFrame(all_cluster_vs_cell_type_prop, columns=cell_types)
print(df_overall)
# rcParams['figure.figsize'] = (30,10)

heap_map_plot = sns.heatmap(df_overall, xticklabels=True, yticklabels=True)

fig = heap_map_plot.get_figure()
fig.savefig("../plots/vis_deconvolution/cluster_vs_deconv_cell_abundance_hmap.pdf", dpi=300)

c_map_plot = sns.clustermap(df_overall, xticklabels=True, yticklabels=True)
c_map_plot.savefig("../plots/vis_deconvolution/cluster_vs_deconv_cell_abundance_cmap.pdf", dpi=300)


# calculate proportion of each cell type per cluster
# filter by cluster, sum cell_type scores and divide by number of clusters create a dataframe

# python vis_deconvolute_analysis.py 