import os
import utils
import pickle
import json
import numpy as np
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from sklearn.metrics import mean_squared_error

import argparse
import warnings
from pathlib import Path
import matplotlib as mpl
import math
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from utils import printmd
import seaborn as sns
############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='vis_nnmf', description='Visium NMF analysis')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name']
############################### BOOOORIING STUFF ABOVE ###############################


random_state = 42
sample_type = "sc"
donor_id = "all_donors"
# markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
# markers = list(set(markers_df["genesymbol"].str.capitalize()))
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
    V = W @ H
        
    # Calculate RMSE of original df and new V
    RMSE = np.sqrt(mean_squared_error(adataX, V))




    return W, H, RMSE



def apply_nmf_on_merged_data(adata, factor_list = range(5,21)):

    factor_rmse_dict = dict()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    

    ct_list =['SC_Differentiated Enterocytes',
    'SC_EC Tumor',
    'SC_Enteroendocrine cell',
    'SC_Goblet cells',
    'SC_Proliferating Enterocytes',
    'SC_Proliferating IECs',
    'SC_Tuft cells']
    ct_list = ['B cells-1',
    'B cells-2',
    'B cells-3',
    'B cells-4',
    'Endothelial cells-1',
    'Endothelial cells-2',
    'Endothelial cells-3',
    'Enteroendocrine cells',
    'Epithelial cells-1',
    'Epithelial cells-2',
    'Epithelial cells-3',
    'Epithelial cells-4',
    'Epithelial cells-5',
    'Epithelial cells-6',
    'Epithelial cells-SI enterocytes-1',
    'Epithelial cells-SI enterocytes-2',
    'Epithelial cells-SI enterocytes-3',
    'Epithelial cells-colon-1',
    'Epithelial cells-colon-2',
    'Epithelial cells-colon-3',
    'Epithelial cells?',
    'Fibroblasts-1',
    'Fibroblasts-2',
    'Fibroblasts-3',
    'Fibroblasts-4',
    'Fibroblasts-5',
    'Goblet cells-1',
    'Goblet cells-2',
    'Goblet cells-3',
    'Goblet cells-4',
    'Goblet cells-5',
    'Goblet cells-6',
    'Goblet cells-7',
    'ILC2',
    'IgA plasma cells-1',
    'IgA plasma cells-2',
    'Keratinocytes',
    'Macrophages',
    'Mast cells',
    'Neuronal cells-1',
    'Neuronal cells-2',
    'Neutrophils',
    'Paneth cells-1',
    'Paneth cells-2',
    'Paneth-Goblet cells',
    'Prolif. cells',
    'Smooth muscle cells-1',
    'Smooth muscle cells-2',
    'T cells',
    'T cells-NK cells',
    'Tuft cells']

    for fact_n in factor_list:
        
    
        print(f"Performing NMF... Number of factors: {fact_n}")
        W, H, RMSE = run_NMF(adata.obs[ct_list].to_numpy(), fact_n, random_state=random_state)
        factor_rmse_dict[f"Factor {fact_n}"] = [RMSE]

        w_factor_col_list = []
        for factor_ind in range(W.shape[1]):
            w_factor_col_list.append(f"{analysis_name}_W{factor_ind+1}")
            adata.obs[f"{analysis_name}_W{factor_ind+1}"] = W[: , factor_ind]

        adata.obs[w_factor_col_list+ct_list].to_csv(os.path.join(output_path,f"{analysis_name}_W_{fact_n}.csv" ))
        h_list = []
        for factor_ind in range(H.shape[0]):
            h_list.append(list(H[factor_ind , :]))
            # adata.var[f"H{fact_n}_{factor_ind+1}"] = H[factor_ind , :]
        h_cols = [f"{analysis_name}_F{i}" for i in range(1, fact_n+1)]
        df_factor_h = pd.DataFrame(h_list, columns = ct_list)
        df_factor_h["factor"] = h_cols
        df_factor_h.set_index("factor", inplace=True)

        df_factor_h.to_csv(os.path.join(output_path,f"{analysis_name}_H_{fact_n}.csv" ))
    # print(factor_rmse_dict)
    df_rmse = pd.DataFrame.from_dict(factor_rmse_dict, orient="index", columns=["rmse"])
    df_rmse.to_csv(os.path.join(output_path, f"{analysis_name}_RMSE.csv"))
    # with open(os.path.join(output_path, f"{analysis_name}_{donor_id}_{fact_n}_RMSE.json"), "w") as outfile: 
    #     json.dump(factor_rmse_dict, outfile)

    return adata
    

def analyse_nmf_results(adata_cc_merged, fact_n, random_state):
    top_cts_df = pd.read_csv(os.path.join(output_path,f"{analysis_name}_H_{fact_n}.csv" ), index_col="factor")
    # normalise each tow sums to 1
    row_sums = top_cts_df.sum(axis=1)
    top_cts_df = top_cts_df.div(row_sums, axis=0)

    """for factor_ind in range(fact_n):
       
        print(f"## Factor {factor_ind}")
        
        # subfigs = fig.subfigures(1, 2, hspace=0, wspace=0, width_ratios=[2, 1])
        plt.rcParams['figure.figsize'] = [30, 15]
        fig, (ax1, ax2) = plt.subplots(1, 2)
        
        sc.pl.spatial(adata_cc_merged, color=f"{analysis_name}_W{factor_ind+1}", cmap="magma_r",size=1.3, img_key='hires', vmax='p99.2', ax = ax1, show=False)

    
        top_cts = top_cts_df.loc[f"{analysis_name}_F{factor_ind+1}"].sort_values(ascending=False)[:20]
        
        
        genes = list(top_cts.keys())
        loadings = list(top_cts.tolist())
        genes.reverse()
        loadings.reverse()
        # print(genes)
        # print(loadings)
        high = math.floor(max(loadings))+1
        low = max(0, min(loadings)-1.01)
        # print(high, low)
        # print([math.ceil(low-0.5*(high-low)), math.ceil(high+0.5*(high-low))])
        plt.xlim([low, high])
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        ax2.barh(genes, loadings, color='grey')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f'{analysis_name}_factor_{factor_ind+1}.pdf'));
        # plt.show();
        print()"""
    

# output_path = "/home/rifaioglu/projects/CRCDiet/plots/visium_nnmf_deconv/"
adata = sc.read_h5ad(input_path)

adata = apply_nmf_on_merged_data(adata)
adata.obs.to_csv(os.path.join(output_path, f"{analysis_name}_obs.csv"))
adata.var.to_csv(os.path.join(output_path, f"{analysis_name}_var.csv"))

# TODO: Seperate plots
plt.rcParams['figure.figsize'] = [10, 10]
df_rmse = pd.read_csv(os.path.join(output_path, f"{analysis_name}_RMSE.csv"), index_col=0).plot.line()
plt.xticks(rotation=90)
plt.savefig(os.path.join(output_path, f"{analysis_name}_RMSE.pdf"))

# This is for deconv
k_dict = {"CD":15,
          "HFD":15,
          "STD_D0":15,
          "STD_D22":15,
          "STD_D43":15,
          "STD_D70":15
        }


analyse_nmf_results(adata, k_dict[analysis_name], 42)
adata.write(os.path.join(f"../data/out_data/", f"{analysis_name}_atlas_deconv_nmf.h5ad"))

"""k_dict = {"CD":7,
          "HFD":7
        }
"""
"""
python visium_nnmf_deconv.py -i /home/rifaioglu/projects/CRCDiet/data/out_data/ext_L24854_CD-AOM-DSS-colon-d81-visium_filtered_clustered_deconv.h5ad -o  /home/rifaioglu/projects/CRCDiet/plots/visium_nnmf_deconv -an CD
python visium_nnmf_deconv.py -i /home/rifaioglu/projects/CRCDiet/data/out_data/ext_L24854_HFD-AOM-DSS-colon-d81-visium_filtered_clustered_deconv.h5ad -o  /home/rifaioglu/projects/CRCDiet/plots/visium_nnmf_deconv -an HFD

python visium_nnmf_deconv.py -i /home/rifaioglu/projects/CRCDiet/data/out_data/all_sample_deconv_STD_D0_15_20_deconv.h5ad -o  /home/rifaioglu/projects/CRCDiet/plots/visium_nnmf_deconv -an STD_D0
python visium_nnmf_deconv.py -i /home/rifaioglu/projects/CRCDiet/data/out_data/all_sample_deconv_STD_D22_15_20_deconv.h5ad -o  /home/rifaioglu/projects/CRCDiet/plots/visium_nnmf_deconv -an STD_D22
python visium_nnmf_deconv.py -i /home/rifaioglu/projects/CRCDiet/data/out_data/all_sample_deconv_STD_D43_15_20_deconv.h5ad -o  /home/rifaioglu/projects/CRCDiet/plots/visium_nnmf_deconv -an STD_D43
python visium_nnmf_deconv.py -i /home/rifaioglu/projects/CRCDiet/data/out_data/all_sample_deconv_STD_D70_15_20_deconv.h5ad -o  /home/rifaioglu/projects/CRCDiet/plots/visium_nnmf_deconv -an STD_D70
"""

# this is for deconv_sc
# python visium_nnmf_deconv.py -i /home/rifaioglu/projects/CRCDiet/data/out_data/ext_L24854_CD-AOM-DSS-colon-d81-visium_filtered_clustered_deconv.h5ad -o  /home/rifaioglu/projects/CRCDiet/plots/visium_nnmf_deconv_sc -an CD
# python visium_nnmf_deconv.py -i /home/rifaioglu/projects/CRCDiet/data/out_data/ext_L24854_HFD-AOM-DSS-colon-d81-visium_filtered_clustered_deconv.h5ad -o  /home/rifaioglu/projects/CRCDiet/plots/visium_nnmf_deconv_sc -an HFD