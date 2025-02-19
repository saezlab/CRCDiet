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
        sns.clustermap(corr, annot=False, cmap=cmap, xticklabels=True, yticklabels=True)
        plt.rcParams['figure.dpi']= 300
        plt.rcParams['figure.figsize']= (120, 120)
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


def extract_cell_type_abundances_adata(adata_path, sample_type, n_of_cells):
    # sample_type = "visium"
    meta = utils.get_meta_data(sample_type)
    for ind, row in meta.iterrows():
        sample_id = row["sample_id"]
        adata = sc.read_h5ad(os.path.join(adata_path, f"all_sample_deconv_{sample_id}_{n_of_cells}_20.h5ad"))
        adata.obsm['q05_cell_abundance_w_sf'].to_csv(f"{adata_path}/cell_type_abundances_{sample_id}_{n_of_cells}_20.csv")


def colocalization_analysis_adata(adata_path, sample_type, exp_name, n_of_cells):

    # colocalization analysis performed based on https://www.nature.com/articles/s41467-021-21892-z#Sec8
    meta = utils.get_meta_data(sample_type)
    
    dict_colocalization = dict()
    for ind, row in meta.iterrows():
        lst_cell_types = []
        sample_id = row["sample_id"]
        condition = row["condition"]
        df_abundance = pd.read_csv(f"{adata_path}/cell_type_abundances_{sample_id}_{n_of_cells}_20.csv", index_col=0)
        for ct in df_abundance.columns:
                ct = ct.split("_")[-1]
                lst_cell_types.append(ct)
        df_abundance.columns = lst_cell_types
        
        corr = df_abundance.corr(method='pearson', min_periods=1, numeric_only=False)
        
        cmap = sns.diverging_palette(230, 20, as_cmap=True)
        sns.clustermap(corr, annot=False, cmap=cmap, xticklabels=True, yticklabels=True)
        plt.rcParams['figure.dpi']= 300
        plt.rcParams['figure.figsize']= (120, 120)
        plt.savefig(f"../plots/{exp_name}/corr_{sample_id}_{n_of_cells}.pdf")

        
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

def colocalization_analysis(sample_id_list):

    # colocalization analysis performed based on https://www.nature.com/articles/s41467-021-21892-z#Sec8
    sample_type = "visium"
    meta = utils.get_meta_data(sample_type)
    
    dict_colocalization = dict()
    #for ind, row in meta.iterrows():
    for condition in sample_id_list:

        lst_cell_types = []
        # sample_id = row["sample_id"]
        # condition = row["condition"]
        df_abundance = pd.read_csv(f"/net/data.isilon/ag-saez/bq_arifaioglu/home/Projects/CRCDiet/data/out_data/cell2location_vis_deconvolution/cell_type_abundances_{condition}_15_20.csv", index_col=0)
        for ct in df_abundance.columns:
                ct = ct.split("_")[-1]
                lst_cell_types.append(ct)
        df_abundance.columns = lst_cell_types
        
        corr = df_abundance.corr(method='pearson', min_periods=1, numeric_only=False)
        plt.rcParams['figure.dpi']= 300
        plt.rcParams['figure.figsize']= (300, 300)
        cmap = sns.diverging_palette(230, 20, as_cmap=True)
        ax = sns.clustermap(corr, annot=False, cmap=cmap, xticklabels=True, yticklabels=True)
        ax.ax_heatmap.set_xticklabels(ax.ax_heatmap.get_xmajorticklabels(), fontsize = 8)
        ax.ax_heatmap.set_yticklabels(ax.ax_heatmap.get_ymajorticklabels(), fontsize = 8)
        # sns.set_theme(font_scale=0.8)
        # ax.ax_cbar.set_ylabel("log2-fold change",size=25);        
        plt.savefig(f"../plots/visium_correlation/corr_{condition}.pdf")

        
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
colocalization_analysis(sample_id_list = ["CD-AOM-DSS", "HFD-AOM-DSS"])

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


# create_adata_file_with_only_epi_cells()

def create_adata_file_with_only_epi_cells():
    # since the epithelial samples includes both epi and stromal i exluded the epithelial cells 
    # with this script. 
    # this is done after the annotation that we were done using both epi plus dn samples. 
    sample_type = "sc_epicells"
    meta = utils.get_meta_data(sample_type)
    samples = np.unique(meta['sample_id'])
    lst_samples = ["CD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN"]
    adata_all_integrated = sc.read_h5ad("../data/out_data/sc_integrated_clustered.h5ad")
    adata_epi_plus_dn = adata_all_integrated[adata_all_integrated.obs["condition"].isin(lst_samples)]
    print(adata_epi_plus_dn)
    adata_epi_plus_dn = adata_epi_plus_dn[adata_epi_plus_dn.obs["major_cell_types"]=="Epithelial cells",:]
    print(adata_epi_plus_dn.obs["condition"])
    adata_immun_epi_cells = sc.read_h5ad("../data/out_data/Parigi_GSE163638_GSM4983265_only-IEC-control_filtered.h5ad")
    # print(adata_immun_epi_cells)
    lst_obs_names = []
    for ind, samp in enumerate(lst_samples):
        # print(samp)
        adata_tmp = adata_epi_plus_dn[adata_epi_plus_dn.obs["condition"]==samp,:]
        tmp_lst_obs_names = list(adata_tmp.obs_names)
        for ind2, o_n in enumerate(tmp_lst_obs_names):
            # print(ind2, o_n)
            obs_id = o_n.split("-")[0] 
            tmp_lst_obs_names[ind2] = f"{obs_id}-1-{ind}"

        lst_obs_names.extend(tmp_lst_obs_names)

    lst_obs_names.extend([f"{obs_n}-3" for obs_n in list(adata_immun_epi_cells.obs_names)])
    print(lst_obs_names) # 
    adata_only_epi_cells = sc.read_h5ad("../data/out_data/sc_epi_plus_DN_cells_integrated_clustered.h5ad")
    adata_only_epi_cells = adata_only_epi_cells[adata_only_epi_cells.obs_names.isin(lst_obs_names),:]
    clusters = ["2", "5", "7", "8", "14", "16"]
    adata_only_epi_cells = adata_only_epi_cells[adata_only_epi_cells.obs["leiden_0.20"].isin(clusters)]
    print(adata_only_epi_cells.obs["condition"])
    adata_only_epi_cells.write("../data/out_data/sc_epicells_only_integrated_clustered.h5ad")

    sc.pl.umap(adata_only_epi_cells, color="leiden_0.20", legend_loc='on data')
    # print(adata_only_epi_cells)

# this is to correct the epithelial samples 
# create_adata_file_with_only_epi_cells()

    
def create_concat_epi_cells():
    sample_type = "sc_epicells"
    meta = utils.get_meta_data(sample_type)
    samples = np.unique(meta['sample_id'])

    # this anndata includes only epithelial cells but the integration is done with DN cells too
    adata_sc_epi_integ_clust = sc.read_h5ad("../data/out_data/sc_epicells_integrated_clustered.h5ad")

    for cond in adata_sc_epi_integ_clust.obs["condition"].cat.categories:
        print(cond, adata_sc_epi_integ_clust[adata_sc_epi_integ_clust.obs["condition"]==cond,:].shape[0]/adata_sc_epi_integ_clust.shape[0])
    
    # this anndata includes both epi and DN
    adata_concat = utils.get_filtered_concat_data(sample_type)
    # remove the DN samples
    adata_concat = adata_concat[adata_sc_epi_integ_clust.obs_names,:]

    for ind in range(4):
        adata_concat.var[f'mt-{ind}'] = adata_concat.var[f'mt-{ind}'].astype(str)
        adata_concat.var[f'rp-{ind}'] = adata_concat.var[f'mt-{ind}'].astype(str)

    adata_concat.write("../data/out_data/sc_epicells_aom_noaom_concatenated.h5ad")

# create_concat_epi_cells()


def create_concat_epi_cells_with_annotations():
    sample_type = "sc_epicells_aom_noaom"
    meta = utils.get_meta_data(sample_type)
    samples = np.unique(meta['sample_id'])

    # this anndata includes only epithelial cells but the integration is done with DN cells too
    adata_sc_epi_integ_clust = sc.read_h5ad("../data/out_data/sc_epicells_integrated_clustered.h5ad")

    for cond in adata_sc_epi_integ_clust.obs["condition"].cat.categories:
        print(cond, adata_sc_epi_integ_clust[adata_sc_epi_integ_clust.obs["condition"]==cond,:].shape[0]/adata_sc_epi_integ_clust.shape[0])
    

    res_param=0.2
    new_res_param=0.2
    sc.tl.leiden(adata_sc_epi_integ_clust, restrict_to=(f'leiden_{res_param:.2f}', ["2"]),  resolution=new_res_param, key_added=f'leiden_{res_param:.2f}')
    sc.pl.umap(adata_sc_epi_integ_clust, color=f'leiden_{res_param:.2f}', palette=sc.pl.palettes.default_20, size=50)
    res_param=0.2
    new_res_param=0.2
    sc.tl.leiden(adata_sc_epi_integ_clust, restrict_to=(f'leiden_{res_param:.2f}', ["2,1"]),  resolution=new_res_param, key_added=f'leiden_{res_param:.2f}')
    sc.pl.umap(adata_sc_epi_integ_clust, color=f'leiden_{res_param:.2f}', palette=sc.pl.palettes.default_20, size=50)
    res_param=0.2
    new_res_param=0.1
    sc.tl.leiden(adata_sc_epi_integ_clust, restrict_to=(f'leiden_{res_param:.2f}', ["2,1,2"]),  resolution=new_res_param, key_added=f'leiden_{res_param:.2f}')
    sc.pl.umap(adata_sc_epi_integ_clust, color=f'leiden_{res_param:.2f}', palette=sc.pl.palettes.default_20, size=50)
    # this anndata includes both epi and DN
    adata_concat = utils.get_filtered_concat_data(sample_type)
    # remove the DN samples
    adata_concat = adata_concat[adata_sc_epi_integ_clust.obs_names,:]

    for ind in range(4):
        adata_concat.var[f'mt-{ind}'] = adata_concat.var[f'mt-{ind}'].astype(str)
        adata_concat.var[f'rp-{ind}'] = adata_concat.var[f'mt-{ind}'].astype(str)
    
    adata_concat.obs["leiden_0.20"] = adata_sc_epi_integ_clust.obs["leiden_0.20"]
    adata_concat.obsm["X_umap"] = adata_sc_epi_integ_clust.obsm["X_umap"]
    annotation_dict = {"2,0":"Proliferating Enterocytes", "2,1,2,1":"Proliferating Enterocytes", "2,1,0":"Differentiated Enterocytes", "2,1,1":"Differentiated Enterocytes", "2,1,2,0":"Differentiated Enterocytes", "2,2":"Differentiated Enterocytes", "2,3":"Differentiated Enterocytes", "2,4":"Differentiated Enterocytes", "2,5":"Differentiated Enterocytes", "2,6":"Differentiated Enterocytes", "5": "Goblet cells", "7":"Proliferating ISCs", "8":"Enteroendocrine cell", "14":"Tuft cells", "16":"EC Tumor" }

    new_cluster_dict = {"2,0":"1", "2,1,2,1":"1", "2,1,0":"0", "2,1,1":"0", "2,1,2,0":"0", "2,2":"0", "2,3":"0", "2,4":"0", "2,5":"0", "2,6":"0", "5": "4", "7":"2", "8":"3", "14":"5", "16":"6" }
    adata_concat.obs[f'cell_type'] = [annotation_dict[clust] for clust in adata_concat.obs[f'leiden_0.20']]
    adata_concat.obs[f'cluster'] = [new_cluster_dict[clust] for clust in adata_concat.obs[f'leiden_0.20']]
    sc.pl.umap(adata_concat, color="cell_type")
    sc.pl.umap(adata_concat, color="cluster")
    adata_concat.write("../data/out_data/sc_epicells_aom_noaom_concatenated_celltype_annot.h5ad")

# create_concat_epi_cells_with_annotations()

