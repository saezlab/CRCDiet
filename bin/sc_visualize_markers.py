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
parser.add_argument('-fn', '--file_name', default="marker_plot", help='File name', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
marker_list = args['marker_list']
analysis_name = args['analysis_name'] # "sc_        "
sample_type = args['sample_type']
file_name = args['file_name']
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
        sc.pl.umap(adata, color=marker_intersect[marker_ind:marker_ind+4], show=False, ncols=len(marker_intersect[marker_ind:marker_ind+4]), save=f'{sample_type}_marker_{mrk_str}')
        marker_ind += 4
        
        
        # fig.tight_layout()

else:
    
    marker_list = [mrk.upper() for mrk in marker_list.split(",")]
    if file_name:
        marker_intersect = list(set(adata.var.index) & set(marker_list))
        sc.pl.umap(adata, color=marker_intersect, size=10, frameon=False, ncols=4, show=False, save=f"{sample_type}_marker_{file_name}")

        plt.clf()
        # print(marker_intersect)
        # plt.rcParams['figure.figsize']=(6,4)
        # ax = sc.pl.umap(adata, color=["WIF1"], color_map=mpl.cm.Reds, size=10, show=False, vmax=8, frameon=False)
        # ax = sc.pl.umap(adata, color=["NKD1"], color_map=mpl.cm.Blues, size=10, show=True, vmax=8, frameon=False)
        # ax = sc.pl.umap(adata, color=["NKD1"], color_map=mpl.cm.Blues, size=10, show=False, vmax=8, frameon=False, ax=ax)
        # ax = sc.pl.umap(adata, color=["AXIN2"], color_map=mpl.cm.Greens, size=10, show=False, vmax=8, frameon=False, ax=ax)
        # ax = sc.pl.umap(adata, color=["NOTUM"], color_map=mpl.cm.Reds, size=10, show=False, vmax=8, frameon=False, ax=ax)
        # ax = sc.pl.umap(adata, color=["PROX1"], color_map=mpl.cm.Reds, size=10, show=False, vmax=8, frameon=False, ax=ax)
        # ax = sc.pl.umap(adata, color=["MMP7"], color_map=mpl.cm.Reds, size=10, show=False, vmax=8, frameon=False, ax=ax)
        # ax = sc.pl.umap(adata, color=["SOX4"], color_map=mpl.cm.Reds, size=10, show=True, vmax=8, frameon=False, ax=ax)

        plot_multiple_markers_bool = True
        if plot_multiple_markers_bool:
        # 10 is the threshold
            umi_thr = 5
            plt.rcParams['figure.dpi']= 300
            plt.rcParams['figure.figsize']= (15, 10)
            bool_over_threshold = (adata[:,adata.var_names.isin(marker_intersect)].layers["counts"]>umi_thr).toarray()
            bool_over_threshold = np.logical_or.reduce(bool_over_threshold, 1).astype(np.int32)
            df_over_threshold = pd.DataFrame(bool_over_threshold, index=adata.obs_names, columns=["Over threshold"])
            df_over_threshold.to_csv(f"{OUT_DATA_PATH}/{sample_type}_tumor_markers_combined_thr_{umi_thr}.csv")
            adata.obs["tumor_markers"] = bool_over_threshold
            sc.pl.umap(adata, color=["tumor_markers"], color_map=mpl.cm.Reds, title="Markers: "+", ".join(marker_intersect), size=40, show=False, frameon=False, save=f"{sample_type}_tumor_markers_combined_thr_{umi_thr}.pdf")

    else:

        for mrk in marker_list:
            
            if mrk in adata.var_names.str.upper():
                print(mrk)
                plt.rcParams['figure.dpi']= 300
                plt.rcParams['figure.figsize']= (15, 10)
                sc.pl.umap(adata, color=mrk.upper(), size=10, title=f"{mrk}", show=False, save=f"{sample_type}_marker_{mrk}")
                
                # fig, axs = plt.subplots(1, 4, figsize=(40, 10));
                # sc.pl.umap(adata, color=mrk.upper(), size=30, title=f"{mrk}: All conditions", show=False, ax=axs[0])
                
                #ind = 1
                #for cond in condition:
                #    if "Immune" in cond:
                #        adata_temp = adata[adata.obs["condition"]==cond,:]        
                #        sc.pl.umap(adata_temp,  title=f"{mrk}: {cond}", color=mrk.upper(), show=False, ax=axs[ind])
                #        ind += 1
                
                # fig.savefig(os.path.join(PLOT_PATH, f"{sample_type}_marker_{mrk}_conditions"));
                # plt.show();

# Uba52,Gm10076,Ubb,Wdr89,Bloc1s1,Tmsb10,Fau,H3f3a

#  python sc_visualize_markers.py -i ../data/out_data/sc_integrated_clustered.h5ad -o ../data/out_data
# python sc_visualize_markers.py -i ../data/out_data/atlas_integrated_clustered.h5ad -o ../data/out_data -an atlas_visualize_markers -st atlas

# python sc_visualize_markers.py -i ../data/out_data/atlas_integrated_clustered.h5ad -o ../data/out_data -an atlas_visualize_markers -st atlas -ml RELN,CD31,PECAM
# python sc_visualize_markers.py -i ../data/out_data/sc_integrated_subclustered_only_tcells.h5ad -o ../data/out_data -an sc_tcells_visualize_markers -st sc
# python sc_visualize_markers.py -i ../data/out_data/sc_epithelial_cells_integrated.h5ad -o ../data/out_data -an sc_epithelial_cells_visualize_markers -st sc -ml CDH2,FN1,VIM,SNAI1,Twist1,Defa,Mptx2,Axin2,Myc,Ccnd1

# python sc_visualize_markers.py -i ../data/out_data/sc_epicells_integrated_clustered.h5ad -o ../data/out_data -an sc_epicells_visualize_markers -st sc_epicells

# python sc_visualize_markers.py -i ../data/out_data/sc_epicells_only_integrated_clustered.h5ad -o ../data/out_data -an sc_epi_cells_aom_noaom_visualize_markers -st sc_epicells -ml Wif1,Nkd1,Axin2,Notum,Prox1,MMP7,Sox4,Ifitm3 -fn tumor_markers
# python sc_visualize_markers.py -i ../data/out_data/sc_epicells_only_integrated_clustered.h5ad -o ../data/out_data -an sc_epi_cells_aom_noaom_visualize_markers -st sc_epicells -ml CDH2,FN1,VIM,SNAI1,Twist1 -fn bonafide_EMT_markers 
# python sc_visualize_markers.py -i ../data/out_data/sc_epicells_only_integrated_clustered.h5ad -o ../data/out_data -an sc_epi_cells_aom_noaom_visualize_markers -st sc_epicells
# python sc_visualize_markers.py -i ../data/out_data/sc_epicells_only_integrated_clustered.h5ad -o ../data/out_data -an sc_epi_cells_aom_noaom_visualize_markers -st sc_epicells -ml Wif1,Nkd1,Axin2,Notum,Prox1,MMP7,Sox4,Ifitm3,CDH2,FN1,VIM,SNAI1,Twist1 -fn all_tumor_markers

# python sc_visualize_markers.py -i ../data/out_data/sc_epicells_integrated_clustered.h5ad -o ../data/out_data -an sc_epi_cells_aom_noaom_visualize_markers -st sc_epicells -ml Wif1,Nkd1,Notum,Prox1,MMP7,FN1 -fn all_tumor_markers_figure

# python sc_visualize_markers.py -i ../data/out_data/sc_epicells_integrated_clustered.h5ad -o ../data/out_data -an sc_epi_cells_aom_noaom_visualize_markers -st sc_epicells -ml Wif1,Nkd1,Axin2,Notum,Prox1,MMP7,Sox4 -fn tumor_markers_together


# Wif1, Nkd1, Axin2, Notum, Prox1, MMP7, Sox4

# python sc_visualize_markers.py -i ../data/out_data/sc_epicells_aom_noaom_integrated_clustered.h5ad -o ../data/out_data -an sc_epi_cells_aom_noaom_visualize_markers_2 -st sc_epicells_aom_noaom -ml WNT6,WNT10A,FZD10,DKK3,WIF1,NKD1,NOTUM,PROX1,MMP7 -fn tumor_markers_together
