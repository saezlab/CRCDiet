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
parser.add_argument('-vp', '--vis_path', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True) # vis_deconvolution
parser.add_argument('-nc', '--num_of_cells', default=8, help='Number of cells per location', required=False)
parser.add_argument('-da', '--detection_alpha', default=20, help='Detection alpha', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
vis_path = args['vis_path']
analysis_name = args['analysis_name'] # vis_deconvolute
n_cells_per_location = int(args['num_of_cells'])
detection_alpha = int(args['detection_alpha'])
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 


#############################################
################## Part 3 ################### 
#############################################

ref_run_name = f'{OUT_DATA_PATH}/cell2location_{analysis_name}'
run_name = f'{OUT_DATA_PATH}/cell2location_map_{analysis_name}'



inf_aver = pd.read_csv(f"{ref_run_name}/inf_aver.csv", index_col=0)
print(inf_aver)
out_fl_name = vis_path.split("/")[-1].split(".")[0]+"_deconv_"+str(n_cells_per_location)+"_"+str(detection_alpha)


adata_vis = utils.get_filtered_concat_data("visium")
# adata_vis = sc.read_h5ad(vis_path)
adata_vis.var.index = pd.Index(gen.upper() for gen in adata_vis.var.index.values)
adata_vis.var_names_make_unique()

intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
print("Intersect:", intersect)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()


print("adata_vis", adata_vis)
print("inf_aver",inf_aver)



# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="condition")


sc.settings.set_figure_params(dpi = 100, color_map = 'viridis', dpi_save = 100,
                              vector_friendly = True, format = 'pdf',
                              facecolor='white')

r = cell2location.run_cell2location(

      # Single cell reference signatures as pd.DataFrame
      # (could also be data as anndata object for estimating signatures
      #  as cluster average expression - `sc_data=adata_snrna_raw`)
      sc_data=inf_aver,
      # Spatial data as anndata object
      sp_data=adata_vis,

      # the column in sc_data.obs that gives cluster idenitity of each cell
      summ_sc_data_args={'cluster_col': "annotation_1",
                        },

      train_args={'use_raw': True, # By default uses raw slots in both of the input datasets.
                  'n_iter': 40000, # Increase the number of iterations if needed (see QC below)

                  # Whe analysing the data that contains multiple experiments,
                  # cell2location automatically enters the mode which pools information across experiments
                  'sample_name_col': 'condition'}, # Column in sp_data.obs with experiment ID (see above)


      export_args={'path': results_folder, # path where to save results
                   'run_name_suffix': '' # optinal suffix to modify the name the run
                  },

      model_kwargs={ # Prior on the number of cells, cell types and co-located groups

                    'cell_number_prior': {
                        # - N - the expected number of cells per location:
                        'cells_per_spot': 8, # < - change this
                        # - A - the expected number of cell types per location (use default):
                        'factors_per_spot': 7,
                        # - Y - the expected number of co-located cell type groups per location (use default):
                        'combs_per_spot': 7
                    },

                     # Prior beliefs on the sensitivity of spatial technology:
                    'gene_level_prior':{
                        # Prior on the mean
                        'mean': 1/2,
                        # Prior on standard deviation,
                        # a good choice of this value should be at least 2 times lower that the mean
                        'sd': 1/4
                    }
      }
)



# python vis_deconvolute_part3.py -rp ../data/out_data/atlas_cell_type_annot_light_weight.h5ad -vp ../data/out_data/ext_L24854_CD-AOM-DSS-colon-d81-visium_filtered.h5ad -an vis_deconvolution 
