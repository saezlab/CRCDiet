from reprlib import aRepr
from pathlib import Path
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.decomposition import NMF
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
meta = utils.get_meta_data(sample_type)
markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
markers = list(set(markers_df["genesymbol"].str.capitalize()))

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
                
                
                mpl.rcParams["image.cmap"]= plt.cm.magma_r
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


sample_type = "sc" 
random_seed= 42                      
# save the cells coming from major cell type as a seperate object
adata_integ_clust= sc.read_h5ad(input_path)

print("Performing NMF... Number of factors: 3")
# W = num of cells x factor
# H = factor x num of genes
W3, H3 = run_NMF(adata_integ_clust.layers["log1p_transformed"], 3, random_state=random_seed)

for factor_ind in range(W3.shape[1]):
    adata_integ_clust.obs[f"W3_{factor_ind+1}"] = W3[: , factor_ind]

for factor_ind in range(H3.shape[0]):
    adata_integ_clust.var[f"H3_{factor_ind+1}"] = H3[factor_ind , :]


print("Performing NMF... Number of factors: 5")
# W = num of cells x factor
# H = factor x num of genes
W5, H5 = run_NMF(adata_integ_clust.layers["log1p_transformed"], 5, random_state=random_seed)

for factor_ind in range(W5.shape[1]):
    adata_integ_clust.obs[f"W5_{factor_ind+1}"] = W5[: , factor_ind]

for factor_ind in range(H5.shape[0]):
    adata_integ_clust.var[f"H5_{factor_ind+1}"] = H5[factor_ind , :]

print("Performing NMF... Number of factors: 20")
W20, H20 = run_NMF(adata_integ_clust.layers["log1p_transformed"], 20, random_state=random_seed)

for factor_ind in range(W20.shape[1]):
    adata_integ_clust.obs[f"W20_{factor_ind+1}"] = W20[: , factor_ind]

for factor_ind in range(H20.shape[0]):
    adata_integ_clust.var[f"H20_{factor_ind+1}"] = H20[factor_ind , :]

# adata_integ_clust.write(os.path.join(OUT_DATA_PATH,f'adata_{analysis_name}_{random_seed}.pckl'))
#utils.write_pickle(os.path.join(OUT_DATA_PATH,f'adata_bcells_merged_nnmf_{random_seed}.pckl'), adata_integ_clust)



analyse_nmf_results(42, adata_integ_clust, 20)
analyse_nmf_results(42, adata_integ_clust, 3)
analyse_nmf_results(42, adata_integ_clust, 5)




# python sc_nnmf.py -i ../data/out_data/sc_immune_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_immunecells_nnmf -ct "B cells"
# python sc_nnmf.py -i ../data/out_data/sc_immune_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_bcells_nnmf -ct "B cells"
# python sc_nnmf.py -i ../data/out_data/sc_stroma_cells_integrated_clustered.h5ad -o ../data/out_data -an sc_stromacells_nnmf
# python sc_pseudobulk_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an sc_test -ct "B cells"