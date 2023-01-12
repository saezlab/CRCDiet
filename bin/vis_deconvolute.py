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
parser.add_argument('-vp', '--vis_path', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True) # vis_deconvolution
parser.add_argument('-of', '--output_file', help='Output file name', required=False)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['vis_path']
analysis_name = args['analysis_name'] # vis_deconvolute
output_file = args['output_file']
sample_type = args['sample_type']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 




#############################################
################## Part 1 ################### 
#############################################

ref_run_name = f'{OUT_DATA_PATH}/cell2location'
run_name = f'{OUT_DATA_PATH}/cell2location_map'

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
                        labels_key='cell_type_level1',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['technology', "study"]
                       )

# create the regression model
mod = RegressionModel(annotated_sc_data)

# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(max_epochs=100, use_gpu=True)
plt.clf()
mod.plot_history()
plt.savefig(f"{ref_run_name}/regression_mod.png",
                   bbox_inches='tight')
plt.clf()
# In this section, we export the estimated cell abundance (summary of the posterior distribution).


# Save model
mod.save(f"{ref_run_name}", overwrite=True)

mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", annotated_sc_data)

annotated_sc_data = mod.export_posterior(
    annotated_sc_data, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)
print(annotated_sc_data.varm["means_per_cluster_mu_fg"])
plt.clf()
mod.plot_QC()
plt.clf()
# Save anndata object with results
adata_file = f"{ref_run_name}/atlas_cell2loc.h5ad"
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

#############################################
################## Part 2 ################### 
#############################################

adata_vis = sc.read_h5ad("/net/data.isilon/ag-saez/bq_arifaioglu/home/Projects/CRCDiet/data/out_data/ext_L24854_HFD-AOM-DSS-colon-d81-visium_filtered.h5ad")
adata_vis.var.index = pd.Index(gen.upper() for gen in adata_vis.var.index.values)
adata_vis.var_names_make_unique()

intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
print("Intersect:", intersect)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

inf_aver = inf_aver[~inf_aver.index.duplicated(keep='first')]

print("adata_vis", adata_vis)
print("inf_aver",inf_aver)
inf_aver.to_csv(f"{ref_run_name}/inf_aver.csv")


# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="condition")


# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver, 
    # the expected average cell abundance: tissue-dependent 
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=8,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
) 
mod.view_anndata_setup()


mod.train(max_epochs=30000, 
          # train using full data (batch_size=None)
          batch_size=None, 
          # use all data points in training because 
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training'])

plt.savefig(f"{ref_run_name}/training_ELBO_history_minus1k.png",
                   bbox_inches='tight')


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file



# python vis_deconvolute.py -rp ../data/out_data/atlas_cell_type_annot_light_weight.h5ad -vp path_to_slide -an  vis_deconvolution -of path_to_out 
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