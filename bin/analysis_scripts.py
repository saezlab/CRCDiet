import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import argparse
import os
import utils
import seaborn as sns
import matplotlib.pyplot as plt
integrated_sample_path = "../data/out_data/integrated.h5ad"

def check_marker_vs_sample_gene_overlap():
    adata = sc.read_h5ad(integrated_sample_path)
    print(adata.var.index)

    df_markers = pd.read_csv("../data/Stroma_oldclusters_DEgenes_Venice_top100xlsx.csv", sep=";")
    marker_genes = []

    # adata_filtered_mouse.var.index = pd.Index(gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.index.values)
    # print([gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.gene_ids.values])
    sample_gene_ids_set = set()
    for item in adata.var.index.values:
        sample_gene_ids_set.add(item.upper())
    print("Number of gene ids in the filtered data:", len(sample_gene_ids_set))

    all_marker_gene_set = set()
    for ind, row in df_markers.iterrows():

        for cell_type in df_markers.columns:
            marker_genes.append([row[cell_type], cell_type])
            all_marker_gene_set.add(row[cell_type])

    print("Number of gene ids in the marker gene set:", len(all_marker_gene_set))

    print("Overlap:", len(sample_gene_ids_set & all_marker_gene_set))

    marker_genes = pd.DataFrame(marker_genes, columns=['gene', 'cell_type']) 

# check_marker_vs_sample_gene_overlap()

"""

"""

def calculate_correlation_between_factors():
    anndata_

def colocalization_analysis():

    # colocalization analysis performed based on https://www.nature.com/articles/s41467-021-21892-z#Sec8
    sample_type = "visium"
    meta = utils.get_meta_data(sample_type)
    
    dict_colocalization = dict()
    for ind, row in meta.iterrows():
        lst_cell_types = []
        sample_id = row["sample_id"]
        condition = row["condition"]
        df_abundance = pd.read_csv(f"../data/out_data/cell2location_map/cell_type_abundances_{sample_id}_filtered_deconv_15_20.csv", index_col=0)
        for ct in df_abundance.columns:
                ct = ct.split("_")[-1]
                lst_cell_types.append(ct)
        df_abundance.columns = lst_cell_types
        
        corr = df_abundance.corr(method='pearson', min_periods=1, numeric_only=False)
        
        cmap = sns.diverging_palette(230, 20, as_cmap=True)
        sns.clustermap(corr, annot=False, cmap=cmap)
        plt.rcParams['figure.dpi']= 300
        plt.rcParams['figure.figsize']= (90, 90)
        plt.savefig(f"../plots/vis_deconvolution/corr_{sample_id}.pdf")

        
        # print(cc) 

        """
        if len(lst_cell_types)==0:
            for ct in df_abundance.columns[1:]:
                ct = ct.split("_")[-1]
                lst_cell_types.append(ct)
                dict_colocalization[ct] = dict()
        
            for ct in lst_cell_types:
                for ct2 in lst_cell_types:
                    dict_colocalization[ct][ct2] = []
                
                    
        
        
        for ind, row in df_abundance.iterrows():
            # print(df_abundance.iloc[ind,1:].sort_values(ascending=False))
            
            total_abundance = sum(df_abundance.iloc[ind,1:])

            df_proportions = df_abundance.iloc[ind,1:]/total_abundance
            # TODO: Try only top k, that is why I am sorting below
            df_proportions = df_proportions.sort_values(ascending=False)
            print(len(df_proportions))
            for index, values in df_proportions.items():
                print(index, values)
            # df_abundance.iloc[ind,1:].set_values( = )df_abundance.iloc[ind,1:].values/total_abundance
            # print(sum(df_abundance.iloc[ind,1:].sort_values(ascending=False)[:5].values))
            # print(sum(df_abundance.iloc[ind,1:].sort_values(ascending=False).values))
            # row[1].arg_sort())
            break
        """
    print(lst_cell_types)    
# colocalization_analysis()


def extract_cell_type_abundances():
    sample_type = "visium"
    meta = utils.get_meta_data(sample_type)
    for ind, row in meta.iterrows():
        
        
        sample_id = row["sample_id"]
        condition = row["condition"]
        adata = sc.read_h5ad(f"../data/out_data/cell2location_map/{sample_id}_filtered_deconv_15_20.h5ad")
        adata.obsm['q05_cell_abundance_w_sf'].to_csv(f"../data/out_data/cell2location_map/cell_type_abundances_{sample_id}_filtered_deconv_15_20.csv")

# extract_cell_type_abundances()

"""
ext_L24854_CD-AOM-DSS-colon-d81-visium
ext_L24854_CD-no-AOM-DSS-colon-d81-visium
ext_L24854_HFD-AOM-DSS-colon-d81-visium
ext_L24854_HFD-no-AOM-DSS-colon-d81-visium
ext_L24854_LFD-AOM-DSS-colon-d81-visium
ext_L24854_LFD-no-AOM-DSS-colon-d81-visium
"""