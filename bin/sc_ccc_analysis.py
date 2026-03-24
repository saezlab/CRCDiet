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
import sys
import warnings
from utils import printmd
import anndata
import utils
import plotnine as p9
import matplotlib as mpl
from copy import deepcopy
import liana as li
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean
from plotnine.scales import scale_x_continuous, scale_x_discrete
from plotnine import ggplot, geom_point, aes, theme, element_text, facet_grid

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
sc.settings.verbosity = 0
warnings.simplefilter(action='ignore')
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Parse arguments
parser = argparse.ArgumentParser(prog='Liana', description='Cell cell communication')
parser.add_argument('-i', '--input_path', help='Path of the anndata object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
parser.add_argument('-gb', '--group_by', help='Group by cell type, condition etc. for comparison', required=False)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # "sc_pse_func_analys"  
group_by = args['group_by']
sample_type = args['sample_type']

# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ############################### 

meta = utils.get_meta_data(sample_type)

adata = sc.read_h5ad(input_path)

print(adata.var_names)
adata.var.index = pd.Index(gen.upper() for gen in adata.var.index.values)


# Run rank_aggregate
li.mt.rank_aggregate(adata, groupby='cell_type', use_raw=False, layer="log1p_transformed", expr_prop=0.1, verbose=True)

print(adata.uns["liana_res"].columns)

my_p = li.pl.dotplot(adata = adata,
                        # colour='magnitude_rank', # "propportion" mran #take the mean olrs faf
                        # colour='magnitude_rank', # "propportion" mran #take the mean olrs faf
                        # inverse_colour=True,
                        colour='lr_means', # "propportion" mran #take the mean olrs faf

                        # size='spec_weight',
                        # size='lr_means',
                        size='magnitude_rank',
                        inverse_size=True,
                        # inverse_size=True, # we inverse sign since we want small p-values to have large sizes
                        source_labels=adata.obs["cell_type"].cat.categories,
                        target_labels=adata.obs["cell_type"].cat.categories,
                        top_n=25,
                        orderby="magnitude_rank",
                        orderby_ascending=True,
                        # finally, since cpdbv2 suggests using a filter to FPs
                        # we filter the pvals column to <= 0.05
                        filterby='cellphone_pvals',
                        filter_lambda=lambda x: x <= 0.05,
                        figure_size=(8, 7)
                        )

adata.uns["liana_res"].to_csv(f"{OUT_DATA_PATH}/{analysis_name}.csv")

"""uniq_pairs = []
for ind, row in adata.uns["liana_res"].iterrows():
    str_temp = row["ligand_complex"]+"->"+row["receptor_complex"]
    if str_temp not in uniq_pairs:
        uniq_pairs.append(str_temp)
        print(str_temp)
    if len(uniq_pairs)==25:
        break"""


my_p = (my_p +
    facet_grid('~source', scales="free_x", space='free_x') +
    p9.theme(strip_text = p9.element_text(size = 6, face="bold", colour = "gray"), 
        axis_text_x=element_text(size=7, face="bold", angle=90))
    )

my_p.save(f"{PLOT_PATH}/ccc.pdf", limitsize = False)

# python sc_ccc_analysis.py -i ../data/out_data/atlas_bcell1_neutrophils_igaplasma.h5ad -o ../data/out_data -st atlas -an atlas_bcell1_neutrophils_igap_ccc -gb cell_type