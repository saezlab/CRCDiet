
import seaborn as sns 
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import argparse
import os
import pickle
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


adata_vis = sc.read_h5ad("/net/data.isilon/ag-saez/bq_arifaioglu/home/Projects/CRCDiet/data/out_data/cell2location_map/ext_L24854_CD-AOM-DSS-colon-d81-visium_filtered_deconv_15_20.h5ad")

adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

cell_types = list(adata_vis.uns['mod']['factor_names'])
condition = adata_vis.obs["condition"].cat.categories[0]

adata_vis_integ = sc.read_h5ad("/net/data.isilon/ag-saez/bq_arifaioglu/home/Projects/CRCDiet/data/out_data/visium_integrated_clustered.h5ad")

adata_vis_integ = adata_vis_integ[adata_vis_integ.obs["condition"]==condition]

adata_vis.obs["leiden_0.10"] = adata_vis_integ.obs["leiden_0.10"]
adata_vis.obs["leiden_0.10"] = adata_vis_integ.obs["leiden_0.10"].values



all_cluster_vs_cell_type_prop = []
for cat in adata_vis.obs["leiden_0.10"].cat.categories:
    adata_temp = adata_vis[adata_vis.obs["leiden_0.10"]==cat]
    lst_cell_type_prop_clust = []
    for ind, val in enumerate(list(adata_temp.obs_names)):
        # print(ind, val)
        raw = adata_temp[ind,:].obs[cell_types].values[0]
        norm = [float(i)/max(raw) for i in raw]
        cell_type_prop = [val/sum(norm) for val in norm]
        lst_cell_type_prop_clust.append(cell_type_prop)
    
    df = pd.DataFrame(lst_cell_type_prop_clust, index=list(adata_temp.obs_names), columns=cell_types)
    row_sum = df.sum()
    row_sum_norm = [float(i)/max(row_sum) for i in row_sum]
    overall_cell_type_prop = [val/sum(row_sum_norm) for val in row_sum_norm]
    all_cluster_vs_cell_type_prop.append(overall_cell_type_prop)
    # all_cluster_vs_cell_type_prop.append(list(row_sum))


correlation_matrix  = np.corrcoef(all_cluster_vs_cell_type_prop) 
print(correlation_matrix)
df_overall = pd.DataFrame(all_cluster_vs_cell_type_prop, index=list(adata_vis.obs["leiden_0.10"].cat.categories), columns=cell_types)
print(df_overall)
heap_map_plot = sns.heatmap(df_overall, xticklabels=True, yticklabels=True)

fig = heap_map_plot.get_figure()

fig.savefig("output.pdf", dpi=300)

c_map_plot = sns.clustermap(df_overall, xticklabels=True, yticklabels=True)

# fig = c_map_plot.get_figure()

c_map_plot.savefig("output_cmap.pdf", dpi=300)


# calculate proportion of each cell type per cluster
# filter by cluster, sum cell_type scores and divide by number of clusters create a dataframe