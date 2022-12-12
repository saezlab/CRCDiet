from reprlib import aRepr
from pathlib import Path
from imageio import save
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os
import warnings
import utils

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run marker visualization')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-ml', '--marker_list', help='Markers seperated by commmas', required=False)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
parser.add_argument('-st', '--sample_type', default="sc", help='Sample type', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
marker_list = args['marker_list']
analysis_name = args['analysis_name'] # "sc_        "
sample_type = args['sample_type']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ###############################

adata_integ_clust = sc.read_h5ad(input_path)

meta = utils.get_meta_data(sample_type)
condition = np.unique(meta['condition'])

# adata = sc.read_h5ad(input_path)

adata = utils.get_filtered_concat_data(sample_type)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

# filter out the cells missing in adata_integ_clust
adata = adata[adata_integ_clust.obs_names,:]

adata.obsm["X_umap"] = adata_integ_clust.obsm["X_umap"]

adata.var.index = pd.Index(gen.upper() for gen in adata.var.index.values)

print("Plotting marker genes on UMAPs...\n")

if marker_list is None:

    markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
    markers = list(set(markers_df["genesymbol"].str.upper()))
    # print(set(markers))


    marker_intersect = list(set(adata.var.index) & set(markers))
    print(f"Number of marker genes: {len(marker_intersect)}")

    marker_ind = 0
    while marker_ind<len(marker_intersect):
        mrk_str = ",".join(marker_intersect[marker_ind:marker_ind+4])
        print(f"Plotting markers: {mrk_str}")
        sc.pl.umap(adata, color=marker_intersect[marker_ind:marker_ind+4], ncols=len(marker_intersect[marker_ind:marker_ind+4]), save=f'{sample_type}_marker_{mrk_str}')
        marker_ind += 4
        
        
        # fig.tight_layout()

else:
    
    marker_list = [mrk.upper() for mrk in marker_list.split(",")]
    
    for mrk in marker_list:
        
        if mrk in adata.var_names.str.upper():
            print(mrk)

            sc.pl.umap(adata, color=mrk.upper(), size=5, title=f"Immune: {mrk}", show=False, save=f"{sample_type}_marker_{mrk}")
            
            fig, axs = plt.subplots(1, 4, figsize=(40, 10));
            sc.pl.umap(adata, color=mrk.upper(), size=30, title=f"{mrk}: All conditions", show=False, ax=axs[0])
            
            ind = 1
            for cond in condition:
                if "Immune" in cond:
                    adata_temp = adata[adata.obs["condition"]==cond,:]        
                    sc.pl.umap(adata_temp,  title=f"{mrk}: {cond}", color=mrk.upper(), show=False, ax=axs[ind])
                    ind += 1
            
            fig.savefig(os.path.join(PLOT_PATH, f"{sample_type}_marker_{mrk}_conditions"));
            plt.show();

# Uba52,Gm10076,Ubb,Wdr89,Bloc1s1,Tmsb10,Fau,H3f3a

#  python sc_visualize_markers.py -i ../data/out_data/sc_integrated_clustered.h5ad -o ../data/out_data
# python sc_visualize_markers.py -i ../data/out_data/atlas_integrated_clustered.h5ad -o ../data/out_data -an atlas_visualize_markers -st atlas