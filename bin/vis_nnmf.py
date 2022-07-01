import os
import utils
import pickle
import numpy as np
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
import argparse
import warnings
from pathlib import Path
import matplotlib as mpl
import math
import matplotlib.pyplot as plt
from utils import printmd

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
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ###############################

sample_type = "visium"
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



def apply_nmf_on_merged_data(random_state):

    adata_cc_merged = get_merged_adata()
    sc.pp.normalize_total(adata_cc_merged, target_sum=1e6)
    sc.pp.log1p(adata_cc_merged)
    
    # print(adata_cc_merged)

    # NMF is applied on the higly variable genes
    # Try without HVG
    # adata_cc_merged = adata_cc_merged[:, adata_cc_merged.var_names.isin(list(merged_object.obsm["pearson_residuals"].columns))] if mode=="hvg" else adata_cc_merged
    # adata_cc_merged.X[adata_cc_merged.X < 0] = 0
    
    # print("NMF 3 Factors!")
    # W = num of cells x factor
    # H = factor x num of genes
    W3, H3 = run_NMF(adata_cc_merged.X, 3, random_state=random_state)
    
    for factor_ind in range(W3.shape[1]):
        adata_cc_merged.obs[f"W3_{factor_ind+1}"] = W3[: , factor_ind]
    
    for factor_ind in range(H3.shape[0]):
        adata_cc_merged.var[f"H3_{factor_ind+1}"] = H3[factor_ind , :]
    
    print("Performing NMF... Number of factors: 20")
    W20, H20 = run_NMF(adata_cc_merged.X, 20, random_state=random_state)

    for factor_ind in range(W20.shape[1]):
        adata_cc_merged.obs[f"W20_{factor_ind+1}"] = W20[: , factor_ind]
    
    for factor_ind in range(H20.shape[0]):
        adata_cc_merged.var[f"H20_{factor_ind+1}"] = H20[factor_ind , :]
    
    utils.write_pickle(os.path.join(OUT_DATA_PATH,f'adata_merged_sct_normalized_nmf_{random_state}.pckl'), adata_cc_merged)
    

def analyse_nmf_results(random_state):

    processed_sample_dict = dict()
    adata_cc_merged = utils.read_pickle(os.path.join(OUT_DATA_PATH,f'adata_merged_sct_normalized_nmf_{random_state}.pckl'))

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
        

        adata = sc.read_h5ad(os.path.join(input_path,f"{sample_id}_filtered.h5ad"))
        adata_samp_cc = adata_cc_merged[adata_cc_merged.obs["condition"]==condition,:].copy()
        adata_samp_cc = adata_samp_cc[:, adata.var_names.intersection(adata_samp_cc.var_names)]
        # print("first", adata_samp_cc) # 3962 × 18308

        adata = adata[:, adata.var_names.intersection(adata_samp_cc.var_names)] #  3962 × 16447 
        # print("second", adata)
        # print(adata_samp_cc)
        adata.X = adata_samp_cc.X
    
        for factor_ind in range(1,4):
            adata.obs[f"W3_{factor_ind}"] = adata_samp_cc.obs[f"W3_{factor_ind}"].values
            adata.var[f"H3_{factor_ind}"] = adata_samp_cc.var[f"H3_{factor_ind}"].values
        
        for factor_ind in range(1,21):
            adata.obs[f"W20_{factor_ind}"] = adata_samp_cc.obs[f"W20_{factor_ind}"].values
            adata.var[f"H20_{factor_ind}"] = adata_samp_cc.var[f"H20_{factor_ind}"].values
        
        processed_sample_dict[sample_id] = adata
            
    
    for factor_ind in range(1,21):
       
        printmd(f"## Factor {factor_ind} <a class='anchor' id='seventh-bullet-1'></a>")
        fig = plt.figure(constrained_layout=True, figsize=(40, 20))
        subfigs = fig.subfigures(1, 2, hspace=0, wspace=0, width_ratios=[2, 1])

        # plt.rcParams["figure.figsize"] = (8, 8)
        rows, cols = (3, 2)
        # fig, ax = plt.subplots(rows, cols, figsize=(25,2))

        axsLeft = subfigs[0].subplots(rows, cols)

        
        for ind, row in meta.iterrows():
            fig_row, fig_col = int(ind/cols), ind%cols
            sample_id = row["sample_id"]
            condition = row["condition"]
            
            mpl.rcParams["image.cmap"]= plt.cm.magma_r
            mpl.rcParams['axes.titlesize'] = 30
            sc.pl.spatial(processed_sample_dict[sample_id], img_key="hires", title=condition, color=f"W20_{factor_ind}", size=1.25, alpha_img=0.3, ax = axsLeft[fig_row][fig_col], show=False)
            cbar = axsLeft[fig_row][fig_col].collections[0].colorbar
            cbar.set_ticks([])
            cbar = None
    
            # fig.tight_layout(pad=1.0)
        # fig.savefig(os.path.join(plot_path, "NMF", str(random_state), f"factor_{factor_ind}_20_all_samples_{mode}.png") , dpi=300)
        # plt.close(fig)
        # fig.tight_layout()

        axsRight = subfigs[1].subplots(1, 1)
        top40_genes = adata_samp_cc.var[f"H20_{factor_ind}"].sort_values(ascending=False)[:40]
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
        plt.savefig(f'{PLOT_PATH}/{sample_type}_factor_{factor_ind}.pdf');
        plt.show();
        print()
        

    """for factor_ind in range(1,4):
        fig, axs = plt.subplots(1, 5, figsize=(20, 2))
        # plt.rcParams["figure.figsize"] = (8, 8)
        samp_num = 0
        for ind, row in meta.iterrows():
            sample_id = row["sample_id"]
            condition = row["condition"]
            
            sc.pl.spatial(processed_sample_dict[sample_id], img_key="hires", title=condition, color=f"W3_{factor_ind}", color_map="bwr", size=1.0, alpha_img=0.5, wspace = 0.1, hspace = 1.0, ax=axs[samp_num], show=True)
            samp_num+=1

        plt.show()

        fig.savefig(os.path.join(PLOT_PATH, "NMF", str(random_state),  f"factor_{factor_ind}_3_all_samples_{mode}.png") , dpi=300)
        plt.close(fig)"""
    

    utils.write_pickle(os.path.join(OUT_DATA_PATH, f'adata_dict_sct_normalized_nmf_{random_state}.pckl'), processed_sample_dict)


def get_merged_adata():
    samples = np.unique(meta['sample_id'])

    adata = []
    # for sample in os.listdir(input_path):
    for sample in samples:

        tmp = sc.read_h5ad(os.path.join(input_path,f"{sample}_filtered.h5ad"))

        # Fetch sample metadata
        m = meta[meta['sample_id'] == sample]
        
        # Add metadata to adata
        for col in m.columns:
            tmp.obs[col] = m[col].values[0]

        # Append
        adata.append(tmp)
        del tmp
        
    # Merge objects and delete list
    adata = adata[0].concatenate(adata[1:], join='outer')
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # keep raw counts in layers
    adata.layers['counts'] = adata.X.copy()
    adata.layers["sqrt_norm"] = np.sqrt(
        sc.pp.normalize_total(adata, inplace=False)["X"]).copy()
    
    # adata.layers["log1p_transformed"] = sc.pp.normalize_total(adata, inplace=False, target_sum=1e6)["X"]
    return adata


apply_nmf_on_merged_data(42)
analyse_nmf_results(42)