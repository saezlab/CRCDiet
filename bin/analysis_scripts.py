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

def calculate_correlation_between_factors(df_path1,df_path2,condition):
    df1 = pd.read_csv(df_path1, index_col=0)
    df2 = pd.read_csv(df_path2, index_col=0)
    df2["condition"] = df1["condition"]
    df1 = df1.loc[df1["condition"].isin([condition]),:]
    df2 = df2.loc[df2["condition"].isin([condition]),:]
    
    
    # mergedStuff = pd.concat([df1, df2], axis=1)
    print(df1)
    print(df2)
    # print(mergedStuff)
    # corr = df1.corrwith(df2)
    # corr = mergedStuff.corr(method='pearson', min_periods=1, numeric_only=False)
    corr_mat = np.ones((len(df1.columns)-1, len(df2.columns)-1))
    print(corr_mat)
    for ind, fac1 in enumerate(df1.columns[:-1]):
        for ind2, pw in enumerate(df2.columns[:-1]):
            corr_mat[ind][ind2] = df1[fac1].corr(df2[pw])
    print(corr_mat)
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    df_corr = pd.DataFrame(corr_mat, index=df1.columns[:-1], columns=df2.columns[:-1])
    print(df_corr)
    sns.clustermap(df_corr, annot=False, cmap=cmap)
    # sns.heatmap(df_corr, annot=False, cmap=cmap)
    plt.rcParams['figure.dpi']= 300
    plt.rcParams['figure.figsize']= (90, 90)
    plt.savefig(f"../plots/corr_nnmf_pathway_{condition}.pdf")
    plt.clf()
    # print(corr)
    

"""calculate_correlation_between_factors("../data/out_data/sc_bcells_nnmf_factors.csv","../data/out_data/sc_bcells_pathway_act_est_mlm_estimate.csv", "LFD-AOM-DSS-Immune")
calculate_correlation_between_factors("../data/out_data/sc_bcells_nnmf_factors.csv","../data/out_data/sc_bcells_pathway_act_est_mlm_estimate.csv", "HFD-AOM-DSS-Immune")
calculate_correlation_between_factors("../data/out_data/sc_bcells_nnmf_factors.csv","../data/out_data/sc_bcells_pathway_act_est_mlm_estimate.csv", "CD-AOM-DSS-Immune")"""

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

def extract_nnmf_weights(adata_path, n_of_factors):
    adata = sc.read_h5ad(adata_path)

    print(adata.obs["condition"])
    weights_arr  = []
    lst_columns = []
    for fact_num in range(1,n_of_factors+1):
        wght = list(adata.obs[f"W{n_of_factors}_{fact_num}"].values)

        weights_arr.append(wght)
        lst_columns.append(f"W{n_of_factors}_{fact_num}")
    
    weights_arr = np.array(weights_arr).reshape(n_of_factors, adata.shape[0]).T
    df_weights = pd.DataFrame(data = weights_arr, 
                  index = adata.obs_names, 
                  columns = lst_columns)
    df_weights["condition"] = adata.obs["condition"]
    
    df_weights.to_csv("../data/out_data/sc_bcells_nnmf_factors.csv")


# extract_nnmf_weights("/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/sc_bcells_nnmf_42.h5ad", 20)  
# adata = sc.read_h5ad("/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/adata_sc_bcells_nnmf_42.h5ad")
# extract_cell_type_abundances()

"""
ext_L24854_CD-AOM-DSS-colon-d81-visium
ext_L24854_CD-no-AOM-DSS-colon-d81-visium
ext_L24854_HFD-AOM-DSS-colon-d81-visium
ext_L24854_HFD-no-AOM-DSS-colon-d81-visium
ext_L24854_LFD-AOM-DSS-colon-d81-visium
ext_L24854_LFD-no-AOM-DSS-colon-d81-visium
"""


def create_adata_file_with_only_epi_cells():
    df_meta_data = pd.read_csv("/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/Parigi_GSE163638_GSM4983265_IEC-Stroma-control_Epi_MetaData.txt")
    lst_barcodes = list(df_meta_data[df_meta_data["orig.ident"]=="control"]["Unnamed: 0"])
    
    for ind, barcode in enumerate(lst_barcodes):
        lst_barcodes[ind] = barcode.split("_")[1]+"-1"
    


    print(lst_barcodes)
    adata_stroma_bcells_ctrl = sc.read_h5ad("/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/Parigi_GSE163638_GSM4983265_IEC-Stroma-control_filtered.h5ad")
    adata_stroma_bcells_ctrl = adata_stroma_bcells_ctrl[adata_stroma_bcells_ctrl.obs_names.isin(lst_barcodes),:]
    adata_stroma_bcells_ctrl.write_h5ad(f"../data/out_data/Parigi_GSE163638_GSM4983265_only-IEC-control_filtered.h5ad")

    


create_adata_file_with_only_epi_cells()