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
from cell2location.utils import select_slide


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
parser.add_argument('-ia', '--infav_path', help='Path to inferred signatures', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True) # vis_deconvolution
parser.add_argument('-nc', '--num_of_cells', default=8, help='Number of cells per location', required=False)
parser.add_argument('-da', '--detection_alpha', default=20, help='Detection alpha', required=False)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False) # visium, visium_helminth, visium_microbiota

args = vars(parser.parse_args())
input_path = args['input_path']
vis_path = args['vis_path']
inf_av_path = args['infav_path']
analysis_name = args['analysis_name'] # vis_deconvolute
n_cells_per_location = int(args['num_of_cells'])
detection_alpha = int(args['detection_alpha'])
sample_type = args['sample_type'] # visium, visium_helminth
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 


#############################################
################## Part 2 ################### 
#############################################

ref_run_name = f'{OUT_DATA_PATH}/cell2location_{analysis_name}'
Path(ref_run_name).mkdir(parents=True, exist_ok=True)

# sample_type = "visium_helminth"

# "../data/out_data/cell2location_atlas/inf_aver.csv"
# "../data/out_data/cell2location_vis_immune_deconvolution/inf_aver.csv"
inf_aver = pd.read_csv(inf_av_path, index_col=0)
print(inf_aver)
# out_fl_name = vis_path.split("/")[-1].split(".")[0]+"_deconv_"+str(n_cells_per_location)+"_"+str(detection_alpha)


adata_vis = sc.read_h5ad(vis_path)
out_fl_name = "all_sample_deconv_"+adata_vis.obs["condition"].cat.categories[0]+"_"+str(n_cells_per_location)+"_"+str(detection_alpha)

# adata_vis = utils.get_filtered_concat_data("visium_aom")
# adata_vis = utils.get_filtered_concat_data("visium_noaom")
# adata_vis = utils.get_filtered_concat_data(sample_type)
# adata_vis.var.index = pd.Index(gen.upper() for gen in adata_vis.var.index.values)

if "korean" in sample_type:
    adata_vis.var_names = adata_vis.var["SYMBOL"]
    
adata_vis.var_names_make_unique()

print("adata_vis.var_names", adata_vis.var_names)
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
print("Intersect:", intersect)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="condition")

print(type(n_cells_per_location))
# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver, 
    # the expected average cell abundance: tissue-dependent 
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=n_cells_per_location,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=detection_alpha
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

plt.savefig(f"{ref_run_name}/{out_fl_name}_training_ELBO_history_minus1k.png",
                   bbox_inches='tight')


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)



# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

for col in adata_vis.var.columns:
    if col.startswith("mt-") or col.startswith("rp-"):
        adata_vis.var[col] = adata_vis.var[col].astype(str)

# Save anndata object with results

adata_file = f"{ref_run_name}/{out_fl_name}.h5ad"
adata_vis.write(adata_file)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

samples = list(adata_vis.obs["condition"].cat.categories)

    
# slide = select_slide(adata_vis,sample, "condition")

lst_cell_types = list(adata_vis.uns['mod']['factor_names'])

sc.pl.spatial(adata_vis, cmap='magma', color=lst_cell_types, ncols=4, size=1.3, img_key='hires', vmin=0, vmax='p99.2', save=f"{out_fl_name}_cell2loc.pdf")

"""meta = utils.get_meta_data(sample_type)

for ind, row in meta.iterrows():
        
    # fig_row, fig_col = int(ind/cols), ind%cols
    sample_id = row["sample_id"]
    condition = row["condition"]
    if condition in samples:
        
        slide = utils.read_filtered_visium_sample(sample_id)
        tmp_adata_vis = adata_vis[adata_vis.obs["condition"]==condition,:]
        sample_id_t = slide.obs["condition"].cat.categories[0]
        out_fl_name = f"{condition}_{n_cells_per_location}_{detection_alpha}_{sample_id_t}"
        
        for c_type in lst_cell_types:
            slide.obs[c_type] = tmp_adata_vis.obs[c_type].values
        with mpl.rc_context({'axes.facecolor':  'black', 'figure.figsize': [4.5, 5]}):
            sc.pl.spatial(slide, cmap='magma', color=lst_cell_types, ncols=4, size=1.3, img_key='hires', vmin=0, vmax='p99.2', library_id=slide.obs["sample"].cat.categories[0], save=f"{out_fl_name}_cell2loc.pdf")


        # lst_cell_types = ["B cells-1", "B cells-2", "B cells-3", "B cells-5", "Endothelial cells-1", "Endothelial cells-2", "Endothelial cells-3", "Enteroendocrine cells", "Epithelial cells-1", "Epithelial cells-2", "Epithelial cells-3", "Epithelial cells-4", "Epithelial cells-5", "Epithelial cells-6", "Epithelial cells-colon-1", "Epithelial cells-colon-2", "Epithelial cells-colon-3", "Epithelial cells-SI enterocytes-1", "Epithelial cells-SI enterocytes-2", "Epithelial cells-SI enterocytes-3", "Epithelial cells?", "Fibroblasts-1", "Fibroblasts-2", "Fibroblasts-3", "Fibroblasts-4", "Fibroblasts-5", "Goblet cells-1", "Goblet cells-2", "Goblet cells-3", "Goblet cells-4", "Goblet cells-5", "Goblet cells-6", "Goblet cells-7", "IgA plasma cells-1", "IgA plasma cells-2", "ILC2", "Keratinocytes", "Macrophages", "Mast cells", "Neuronal cells-1", "Neuronal cells-2", "Neutrophils", "Paneth cells-1", "Paneth cells-2", "Paneth-Goblet cells", "Prolif. cells", "Smooth muscle cells-1", "Smooth muscle cells-2", "T cells", "T cells-NK cells", "Tuft cells"]
        # lst_cell_types = ["B cells", "Dendritic cells", "ILC2", "Mast cells", "Myeloid cells", "Neutrophils", "Plasma cells", "T cells"]

        # plot in spatial coordinates
        


"""


# to adata.obs with nice names for plotting




# python vis_deconvolute_part2.py -rp ../data/out_data/atlas_cell_type_annot_light_weight.h5ad -vp ../data/out_data/ext_L24854_CD-AOM-DSS-colon-d81-visium_filtered.h5ad -an  vis_deconvolution 
"""
import cv2
import matplotlib.pyplot as plt
import pandas as pd
import json
cv2.startWindowThread()


img = cv2.imread('TileScan 3_Region3_HFD_No_AOMDSS_Merged_RAW_ch00.tif')
# img = cv2.imread('V1_Adult_Mouse_Brain_image.tif')
columns = ["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres", "pxl_col_in_fullres"]
df_locs = pd.read_csv("tissue_positions_list.csv", header=None)
df_locs.columns = columns
json_file = open("scalefactors_json.json", "r")
scale_dict = data = json.loads(json_file.read())
radius = round(scale_dict["spot_diameter_fullres"]/2)+1


for ind, row in df_locs.iterrows():
    # print(row)
    if row["in_tissue"]==1:
    
    	x_coor = row["pxl_col_in_fullres"]
    	y_coor = row["pxl_row_in_fullres"]
    	cv2.circle(img,(x_coor,y_coor), radius, (0,0,255), -1)
    	barcode= row["barcode"]

    	cropped = img[y_coor-(radius+5):y_coor+(radius+5),  x_coor-(radius+5):x_coor+(radius+5)]
    	cv2.imwrite(f"{barcode}.tif", cropped)

		
    


# cv2.imshow('imag2e',img)
cv2.imwrite("deneme.tif", img)
cv2.waitKey(1000)
"""

# sbatch --job-name=deconv_12-01-2023_1 -p gpu --gres=gpu:1 --mem=30g  -n 1 --time=7-00:00:00 --output=deconv.out "deconv.sh"
# sbatch --job-name=deconv_08-07-2024_2 -p gpu --gres=gpu:1 --mem=20 -n 1 --time=7-00:00:00 --output=vis_deconvolution_helminth_stdiet_15_1.out "vis_deconvolution_helminth_stdiet_15_1.sh"
# sbatch --job-name=deconv_08-07-2024_1 -p gpu --gres=gpu:1 --mem=20 -n 1 --time=7-00:00:00 --output=vis_deconvolution_helminth_stdiet_15_2.out "vis_deconvolution_helminth_stdiet_15_2.sh"
# sbatch --job-name=deconv_08-07-2024_3 -p gpusaez --gres=gpu:1 --mem=20 -n 1 --time=7-00:00:00 --output=vis_deconvolution_helminth_stdiet_10_1.out "vis_deconvolution_helminth_stdiet_10_1.sh"
# sbatch --job-name=deconv_08-07-2024_5 -p gpusaez --gres=gpu:1 --mem=20 -n 1 --time=7-00:00:00 --output=vis_deconvolution_helminth_stdiet_10_2.out "vis_deconvolution_helminth_stdiet_10_2.sh"