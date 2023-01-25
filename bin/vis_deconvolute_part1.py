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
import matplotlib as mpl
import argparse
import os
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes
from sklearn.metrics import silhouette_score, pairwise_distances
import pickle
import utils
import math
import warnings
import cell2location
import scvi
import argparse


############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
# sc.set_figure_params(scanpy=True, facecolor="white", dpi=50, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='cluster', description='Run deconvolution')
parser.add_argument('-rp', '--input_path', help='path to referecence atlas', required=True)
parser.add_argument('-vp', '--vis_path', help='Output directory where to store the object', required=False)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True) # vis_deconvolution
parser.add_argument('-of', '--output_file', help='Output file name', required=False)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
vis_path = args['vis_path']
analysis_name = args['analysis_name'] # vis_deconvolute
output_file = args['output_file']
sample_type = args['sample_type']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 




#############################################
################## Part 1 ################### 
#############################################

ref_run_name = f'{OUT_DATA_PATH}/cell2location_{analysis_name}'
# run_name = f'{OUT_DATA_PATH}/cell2location_map_{analysis_name}'

annotated_sc_data = sc.read_h5ad(input_path)
annotated_sc_data.var_names_make_unique()

# calculate qc metrics
annotated_sc_data.var["mt"] = annotated_sc_data.var_names.str.contains("^MT-")


# remove mitochondrial genes
annotated_sc_data = annotated_sc_data[:,~annotated_sc_data.var["mt"]]
selected = filter_genes(annotated_sc_data, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
print("adata_ref" ,annotated_sc_data)


# filter the object
annotated_sc_data = annotated_sc_data[annotated_sc_data.obs_names, selected].copy()

# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=annotated_sc_data,
                        # 10X reaction / sample / batch
                        batch_key='batch',
                        # cell type, covariate used for constructing signatures
                        labels_key='cell_type_0.20',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        # categorical_covariate_keys=['technology', "study"]
                       )

# create the regression model
mod = RegressionModel(annotated_sc_data)

# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(max_epochs=250, use_gpu=True) # check qc and change max_epochs 
plt.clf()
mod.plot_history()
plt.savefig(f"{ref_run_name}/regression_mod.png",
                   bbox_inches='tight')
plt.clf()
# In this section, we export the estimated cell abundance (summary of the posterior distribution).


# Save model
mod.save(f"{ref_run_name}", overwrite=True)

#  mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", annotated_sc_data)

annotated_sc_data = mod.export_posterior(
    annotated_sc_data, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)
print(annotated_sc_data.varm["means_per_cluster_mu_fg"])
plt.clf()
mod.plot_QC()
plt.clf()
# Save anndata object with results
adata_file = f"{ref_run_name}/cell2loc.h5ad"
annotated_sc_data.write(adata_file)
adata_file
plt.clf()
mod.plot_QC()

plt.savefig(f"{ref_run_name}/c2l_qc.png",
                   bbox_inches='tight')

plt.clf()


# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in annotated_sc_data.varm.keys():
    inf_aver = annotated_sc_data.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in annotated_sc_data.uns['mod']['factor_names']]].copy()
else:
    inf_aver = annotated_sc_data.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in annotated_sc_data.uns['mod']['factor_names']]].copy()
inf_aver.columns = annotated_sc_data.uns['mod']['factor_names']

inf_aver.index = inf_aver.index.str.upper()
# find shared genes and subset both anndata and reference signatures

inf_aver = inf_aver[~inf_aver.index.duplicated(keep='first')]
inf_aver.to_csv(f"{ref_run_name}/inf_aver.csv")

# python vis_deconvolute_part1.py -rp ../data/out_data/atlas_cell_type_annot_light_weight.h5ad -vp ../data/out_data/ext_L24854_CD-AOM-DSS-colon-d81-visium_filtered.h5ad -an  vis_deconvolution -of path_to_out

# python vis_deconvolute_part1.py -rp ../data/out_data/atlas_cell_type_annot_light_weight.h5ad -vp ../data/out_data/ext_L24854_CD-AOM-DSS-colon-d81-visium_filtered.h5ad -an  vis_deconvolution -of path_to_out 

# python vis_deconvolute_part1.py -rp "../data/out_data/sc_Immune cells.h5ad" -an vis_immune_deconvolution
"""
import cv2
import matplotlib.pyplot as plt
import pandas as pd
cv2.startWindowThread()


img = cv2.imread('TileScan 3_Region3_HFD_No_AOMDSS_Merged_RAW_ch00.tif')
columns = ["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres", "pxl_col_in_fullres"]
df_locs = pd.read_csv("tissue_positions_list.csv", header=None)
df_locs.columns = columns


for ind, row in df_locs.iterrows():
    # print(row)
    x_coor = row["pxl_row_in_fullres"]
    y_coor = row["pxl_col_in_fullres"]
    cv2.circle(img,(x_coor,y_coor), 33, (0,0,255), -1)
    


cv2.imshow('imag2e',img)
cv2.waitKey(1000)
"""