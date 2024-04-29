from sklearn.preprocessing import MinMaxScaler
from reprlib import aRepr
from pathlib import Path
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.decomposition import NMF
from matplotlib import cm
import argparse
import math
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
parser = argparse.ArgumentParser(prog='nnmfbc', description='NNMF Analysis on B cells')
parser.add_argument('-i', '--input_path', help='Input path to integrated object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
# parser.add_argument('-ct', '--cell_type', help="Cell type of interest", required=False)
args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # subset
# m_ct = args['cell_type']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name, plot=True)
############################### BOOOORIING STUFF ABOVE ###############################
                     
sample_type = "sc"
# Load meta data
"""meta = utils.get_meta_data(sample_type)
markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
markers = list(set(markers_df["genesymbol"].str.capitalize()))"""

def run_NMF(adataX, n_components=2, random_state=0):
    """Run NMF on adata

    Run NMF on adata given number of components and random state

    Args:
        adataX (np.array): num of cells x num of genes
        n_components (int): number of components for NMF
        random_state (int): random seed
    
    Returns:
        Metagenes and factors

    """
    model = NMF(n_components, init='random', random_state=random_state)
    W = model.fit_transform(adataX)
    H = model.components_

    return W, H


def analyse_nmf_results(random_state, adata_nnmf_merged, n_of_factors=20):

    processed_sample_dict = dict()
    # adata_nnmf_merged = sc.read_h5ad(os.path.join(OUT_DATA_PATH, f'adata_{analysis_name}_nnmf_{random_seed}.pckl'))
    # adata_filt_concat = utils.get_filtered_concat_data(sample_type)
    # adata_filt_concat = adata_filt_concat[adata_bcells_merged.obs_names,:]
    sample_num = 0
    is_first = True
    for ind, row in meta.iterrows():
        sample_id = row["sample_id"]
        condition = row["condition"]
        # print(sample_id, condition)
        # print(adata_cc_merged[adata_cc_merged.obs["condition"]==condition,:])
        
        if is_first:
            is_first=False
        else:
            sample_num+=1 
        row = meta[meta["sample_id"]==sample_id]
    
        # adata = sc.read_h5ad(os.path.join(input_path,f"{sample_id}_filtered.h5ad"))
        adata_bcells_cond = adata_nnmf_merged[adata_nnmf_merged.obs["condition"]==condition,:].copy()        
        processed_sample_dict[sample_id] = adata_bcells_cond
            
    
    for factor_ind in range(1,(n_of_factors+1)):
       
        # printmd(f"## Factor {factor_ind} <a class='anchor' id='seventh-bullet-1'></a>")
        fig = plt.figure(constrained_layout=True, figsize=(20, 20))
        subfigs = fig.subfigures(1, 2, hspace=0, wspace=0, width_ratios=[1, 1])

        # plt.rcParams["figure.figsize"] = (8, 8)
        rows, cols = (3, 1)
        # fig, ax = plt.subplots(rows, cols, figsize=(25,2))

        axsLeft = subfigs[0].subplots(rows, cols)

        ind=0
        for _, row in meta.iterrows():
            sample_id = row["sample_id"]
            condition = row["condition"]
            
            # change this condition based on the sample
            # TODO: Refactor this part
            # if "Epi" in condition or True:
            # To visualize all the samples: if "Epi" in condition or True:
            if "Immune" in condition:
            
                fig_row, fig_col = int(ind/cols), ind%cols
                
                magma_r = cm.get_cmap('magma_r', 12)
                plt.set_cmap(magma_r)
                mpl.rcParams['axes.titlesize'] = 30
                # sc.pl.umap(processed_sample_dict[sample_id], title=condition, color=f"W20_{factor_ind}", ax = axsLeft[fig_row][fig_col], show=False)
                # cbar = axsLeft[fig_row][fig_col].collections[0].colorbar
                sc.pl.umap(processed_sample_dict[sample_id], size=30, title=condition, color=f"W{n_of_factors}_{factor_ind}", ax = axsLeft[fig_row], show=False)
                cbar = axsLeft[fig_row].collections[0].colorbar
                cbar.set_ticks([])
                cbar = None
                ind+=1
    
            # fig.tight_layout(pad=1.0)
        # fig.savefig(os.path.join(plot_path, "NMF", str(random_state), f"factor_{factor_ind}_20_all_samples_{mode}.png") , dpi=300)
        # plt.close(fig)
        # fig.tight_layout()

        axsRight = subfigs[1].subplots(1, 1)
        top40_genes = adata_nnmf_merged.var[f"H{n_of_factors}_{factor_ind}"].sort_values(ascending=False)[:40]
        # print(top40_genes)
        genes = list(top40_genes.keys())
        loadings = list(top40_genes.tolist())
        genes.reverse()
        loadings.reverse()
        # print(genes)
        # print(loadings)
        high = math.floor(max(loadings))+1
        low = max(0, min(loadings)-1.01)
        # print(high, low)
        # print([math.ceil(low-0.5*(high-low)), math.ceil(high+0.5*(high-low))])
        plt.xlim([low, high])
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        axsRight.barh(genes, loadings, color='grey')
        plt.savefig(f'{PLOT_PATH}/{sample_type}_{n_of_factors}_factor_{factor_ind}.pdf');


def analyse_nmf_single_sample(random_state, adata_nnmf_merged, n_of_factors=20, major_cell_type=None):
    processed_sample_dict = dict()
    # adata_nnmf_merged = sc.read_h5ad(os.path.join(OUT_DATA_PATH, f'adata_{analysis_name}_nnmf_{random_seed}.pckl'))
    # adata_filt_concat = utils.get_filtered_concat_data(sample_type)
    # adata_filt_concat = adata_filt_concat[adata_bcells_merged.obs_names,:]
    
    for factor_ind in range(1,(n_of_factors+1)):
       
        # printmd(f"## Factor {factor_ind} <a class='anchor' id='seventh-bullet-1'></a>")
        
        print(factor_ind)
        print(f'{PLOT_PATH}/{sample_type}_{n_of_factors}_factor_{factor_ind}.pdf')
        fig, (axsLeft, axsRight)  = plt.subplots(1, 2, figsize=(60, 30), width_ratios=[2, 1])

        # axsLeft = subfigs[0]
        
        magma_r = cm.get_cmap('magma_r', 12)
        plt.set_cmap(magma_r)
        mpl.rcParams['axes.titlesize'] = 30
        # sc.pl.umap(processed_sample_dict[sample_id], title=condition, color=f"W20_{factor_ind}", ax = axsLeft[fig_row][fig_col], show=False)
        # cbar = axsLeft[fig_row][fig_col].collections[0].colorbar
        sc.pl.umap(adata_nnmf_merged, size=40, color=f"W{n_of_factors} {factor_ind}", ax = axsLeft, show=False)
        cbar = axsLeft.collections[0].colorbar
        cbar.set_ticks([])
        cbar = None

        #axsRight = subfigs[1]
        top40_genes = adata_nnmf_merged.var[f"H{n_of_factors} {factor_ind}"].sort_values(ascending=False)[:40]
        # print(top40_genes)
        genes = list(top40_genes.keys())
        loadings = list(top40_genes.tolist())
        genes.reverse()
        loadings.reverse()
        # print(genes)
        # print(loadings)
        high = math.floor(max(loadings))+1
        low = max(0, min(loadings)-1.01)
        # print(high, low)
        # print([math.ceil(low-0.5*(high-low)), math.ceil(high+0.5*(high-low))])
        plt.xlim([low, high])
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        axsRight.barh(genes, loadings, color='grey')
        if major_cell_type:
            plt.savefig(f'{PLOT_PATH}/{sample_type}_{major_cell_type}_{n_of_factors}_factor_{factor_ind}.pdf');
        else:
            plt.savefig(f'{PLOT_PATH}/{sample_type}_{n_of_factors}_factor_{factor_ind}.pdf');

def generate_heatmap(adata, n_of_factors=3, major_cell_type=None):

    from PyComplexHeatmap import ClusterMapPlotter
    # group by cell type, calculate mean of each factor
    immune_cell_types = [ "T cells", "Plasma cells", "Mast cells", "Neutrophils", "ILC2",  "Dendritic cells", "Myeloid cells", "B cells"]
    epithelial_cell_types = ["Tuft cells", "Goblet cells", "Prolif.", "Enteroendocrine", "Keratynocytes", "Prolif. + Mature enterocytes"]
    stroma_cell_types = ["Stroma", "Myofibroblasts", "Endothelial cells"]


    sns.set_style("whitegrid", {'axes.grid' : False})
    

    factor_columns = [f"W{n_of_factors} {i}" for i in range(1, n_of_factors+1)]
    cell_type_summarised_data = adata.obs.groupby(["cell_type_0.20"]).mean()
    cell_type_summarised_data = cell_type_summarised_data[factor_columns]
    factor_columns = [f"Factor {i}" for i in range(1, n_of_factors+1)]
    scaler = MinMaxScaler(feature_range=(-2,2))
    scaler.fit(cell_type_summarised_data)
    scaled_factors = scaler.transform(cell_type_summarised_data)
    # sns.clustermap(scaled_factors.T, cmap="viridis", xticklabels=cell_type_summarised_data.index.values, yticklabels=factor_columns)
    # plt.savefig(f'{PLOT_PATH}/{sample_type}_{n_of_factors}_factors_heatmap_celltype.pdf');

    ct_df = pd.DataFrame(scaled_factors.T, index=factor_columns,   columns=cell_type_summarised_data.index.values)
    
    factor_columns = [f"W{n_of_factors} {i}" for i in range(1, n_of_factors+1)]
    cond_summarised_data = adata.obs.groupby(["condition"]).mean()
    cond_summarised_data = cond_summarised_data[factor_columns]
    factor_columns = [f"Factor {i}" for i in range(1, n_of_factors+1)]
    scaler = MinMaxScaler(feature_range=(-2,2))
    scaler.fit(cond_summarised_data)
    scaled_factors = scaler.transform(cond_summarised_data)
    # sns.clustermap(scaled_factors.T, cmap="viridis", xticklabels=cond_summarised_data.index.values, yticklabels=factor_columns)
    # plt.savefig(f'{PLOT_PATH}/{sample_type}_{n_of_factors}_factors_heatmap_sample.pdf');
    cond_df = pd.DataFrame(scaled_factors.T, index=factor_columns,   columns=cond_summarised_data.index.values)

    merged_df = pd.concat([ct_df, cond_df], axis=1) # pd.merge(ct_df, cell_type_summarised_data, right_index=True) 
    sep_df = pd.DataFrame(['GroupA'] * len(ct_df.columns) + ['GroupB']*len(cond_df.columns), index= merged_df.columns, columns=['AB'])

    plt.figure(figsize=(16, 4))
    if n_of_factors==20:
        plt.figure(figsize=(16, 16))

    cm = ClusterMapPlotter(
        data=merged_df,
        col_cluster=True,row_cluster=True,
        col_split=sep_df.AB,
        col_split_gap=1.5,row_split_gap=0.8,
        label='values',row_dendrogram=True,
        show_rownames=True,show_colnames=True,
        tree_kws={'row_cmap': 'Set1'},verbose=0,legend_gap=5,
        cmap='viridis',xticklabels_kws={'labelrotation':-90,'labelcolor':'black'},
        ylabel="Factors",)
    if major_cell_type:
        
        plt.savefig(f'{PLOT_PATH}/{sample_type}_{major_cell_type}_{n_of_factors}_factors_complex_heatmap.pdf', bbox_inches='tight', dpi=300);
    else:   
        plt.savefig(f'{PLOT_PATH}/{sample_type}_{n_of_factors}_factors_complex_heatmap.pdf', bbox_inches='tight', dpi=300);




def run_NMF_on_adata(adata, n_of_factors=3):

    # save the cells coming from major cell type as a seperate object
    
    # print(adata_integ_clust.obs["cell_type_0.20"])
    print(f"Performing NMF... Number of factors: {n_of_factors}")
    # W = num of cells x factor
    # H = factor x num of genes
    W3, H3 = run_NMF(adata.layers["log1p_transformed"], n_of_factors, random_state=random_seed)

    for factor_ind in range(W3.shape[1]):
        adata.obs[f"W{n_of_factors} {factor_ind+1}"] = W3[: , factor_ind]

    for factor_ind in range(H3.shape[0]):
        adata.var[f"H{n_of_factors} {factor_ind+1}"] = H3[factor_ind , :]

    return adata


random_seed= 42 


major_cell_types = ["epi", "stroma", "immune"]
major_cell_types = ["immune"]
factor_list = [3, 4, 5, 20]

# adata_integ_clust= sc.read_h5ad(input_path)
adata_integ_clust = sc.read_h5ad(os.path.join(OUT_DATA_PATH,f'{analysis_name}_{random_seed}.h5ad'))

for factor in factor_list:
    
    adata_integ_clust = run_NMF_on_adata(adata_integ_clust, n_of_factors=factor)
    analyse_nmf_single_sample(random_seed, adata_integ_clust, factor)
    generate_heatmap(adata_integ_clust, factor)


for mct in major_cell_types:
    adata_integ_clust = sc.read_h5ad(os.path.join(OUT_DATA_PATH,f'{analysis_name}_{random_seed}.h5ad'))
    for factor in factor_list:
        adata_integ_clust = utils.get_subset_adata(adata=adata_integ_clust, cell_types=mct)
        print(adata_integ_clust)
        adata_integ_clust = run_NMF_on_adata(adata_integ_clust, n_of_factors=factor)
        analyse_nmf_single_sample(random_seed, adata_integ_clust, factor, major_cell_type=mct)
        generate_heatmap(adata_integ_clust, factor, major_cell_type=mct)








# python sc_nnmf.py -i ../data/out_data/sc_immune_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_immunecells_nnmf -ct "B cells"
# python sc_nnmf.py -i ../data/out_data/sc_immune_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_bcells_nnmf -ct "B cells"
# python sc_nnmf.py -i ../data/out_data/sc_stroma_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_stromacells_nnmf
# python sc_nnmf.py -i ../data/out_data/atlas_bcell_populations.h5ad -o ../data/out_data -an atlas_bcell_populations_nnmf

# python sc_pseudobulk_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an sc_test -ct "B cells"
# python sc_nnmf.py -i ../data/out_data/atlas_bcells_subclustered_integrated_subclustered.h5ad -o  ../data/out_data -an atlas_bcell_populations_nnmf
# python sc_nnmf.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o  ../data/out_data -an sc_nnmf