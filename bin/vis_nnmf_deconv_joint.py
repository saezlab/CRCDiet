import os
import utils
import pickle
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
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
############################### BOOOORIING STUFF ABOVE ###############################

sample_type = "visium"
# Load meta data
meta = utils.get_meta_data("visium_diet_kinetics_AOM_DSS")
# markers_df = pd.read_csv(os.path.join(DATA_PATH, "marker_genes.txt"), sep="\t")
# markers = list(set(markers_df["genesymbol"].str.capitalize()))
l_param = 0.15
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



def apply_nmf_on_merged_data(adata, factor_list = range(5,21), random_state=42):

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
            w_factor_col_list.append(f"W{fact_n}_{factor_ind+1}")
            # adata.obs[f"{analysis_name}_W{factor_ind+1}"] = W[: , factor_ind]
            adata.obs[f"W{fact_n}_{factor_ind+1}"] = W[: , factor_ind]

        adata.obs[w_factor_col_list+ct_list].to_csv(os.path.join(output_path,f"{analysis_name}_W_{fact_n}.csv" ))
        h_list = []
        for factor_ind in range(H.shape[0]):
            h_list.append(list(H[factor_ind , :]))
            # adata.var[f"H{fact_n}_{factor_ind+1}"] = H[factor_ind , :]
        # f"W{fact_n}_{factor_ind+1}"
        h_cols = [f"H{fact_n}_{i}" for i in range(1, fact_n+1)]
        print(h_cols)
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
    processed_sample_dict = dict()
    sample_num = 0
    is_first = True
    for ind, row in meta.iterrows():
        sample_id = row["sample_id"]
        condition = row["condition"]
        print(sample_id, condition)
        # print(adata_cc_merged[adata_cc_merged.obs["condition"]==condition,:])
        
        if is_first:
            is_first=False
        else:
            sample_num+=1 
        row = meta[meta["sample_id"]==sample_id]
        

        adata = sc.read_h5ad(os.path.join(input_path,f"{sample_id}.h5ad"))
        adata_samp_cc = adata_cc_merged[adata_cc_merged.obs["condition"]==condition,:].copy()
        adata_samp_cc = adata_samp_cc[:, adata.var_names.intersection(adata_samp_cc.var_names)]
        # print("first", adata_samp_cc) # 3962 × 18308

        adata = adata[:, adata.var_names.intersection(adata_samp_cc.var_names)] #  3962 × 16447 
        # print("second", adata)
        # print(adata_samp_cc)
        # adata.X = adata_samp_cc.X
    
        """for factor_ind in range(1,4):
            adata.obs[f"W3_{factor_ind}"] = adata_samp_cc.obs[f"W3_{factor_ind}"].values
            adata.var[f"H3_{factor_ind}"] = adata_samp_cc.var[f"H3_{factor_ind}"].values

        for factor_ind in range(1,6):
            adata.obs[f"W5_{factor_ind}"] = adata_samp_cc.obs[f"W5_{factor_ind}"].values
            adata.var[f"H5_{factor_ind}"] = adata_samp_cc.var[f"H5_{factor_ind}"].values"""
        
        for factor_ind in range(1,fact_n+1):
            adata.obs[f"W{fact_n}_{factor_ind}"] = adata_samp_cc.obs[f"W{fact_n}_{factor_ind}"].values
            # adata.var[f"H20_{factor_ind}"] = adata_samp_cc.var[f"H20_{factor_ind}"].values
        top_cts_df = pd.read_csv(os.path.join(output_path,f"{analysis_name}_H_{fact_n}.csv" ), index_col="factor")
        # normalise each tow sums to 1
        row_sums = top_cts_df.sum(axis=1)
        top_cts_df = top_cts_df.div(row_sums, axis=0)
        processed_sample_dict[sample_id] = (adata, top_cts_df)
        
    
    for factor_ind in range(1,fact_n+1):
        
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
                
                # mpl.rcParams["image.cmap"]= mpl.colormaps['viridis']# plt.cm.magma_r
                mpl.rcParams['axes.titlesize'] = 30
                sc.pl.spatial(processed_sample_dict[sample_id][0], img_key="hires", title=condition, color=f"W{fact_n}_{factor_ind}", cmap="magma_r", size=1.25, alpha_img=0.3, ax = axsLeft[fig_row][fig_col], show=False)
                cbar = axsLeft[fig_row][fig_col].collections[0].colorbar
                cbar.set_ticks([])
                cbar = None
        
                # fig.tight_layout(pad=1.0)
            # fig.savefig(os.path.join(PLOT_PATH, "NMF", str(random_state), f"factor_{factor_ind}_20_all_samples.pdf") , dpi=300)
            # plt.close(fig)
            # fig.tight_layout()

            axsRight = subfigs[1].subplots(1, 1)
            
            top_cts = processed_sample_dict[sample_id][1].loc[f"H{fact_n}_{factor_ind}"].sort_values(ascending=False)[:20]
        
        
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
            axsRight.barh(genes, loadings, color='grey')
            plt.savefig(f'{PLOT_PATH}/{sample_type}_factor_{factor_ind}.pdf');
            # plt.show();
            print()
    

def get_merged_adata():
    samples = np.unique(meta['sample_id'])

    adata = []
    # for sample in os.listdir(input_path):
    for sample in samples:

        tmp = sc.read_h5ad(os.path.join(input_path,f"{sample}.h5ad"))

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
    
    return adata


def generate_heatmap(adata, n_of_factors=3, major_cell_type=None):

    from PyComplexHeatmap import ClusterMapPlotter
    sns.set_style("whitegrid", {'axes.grid' : False})
    

    factor_columns = [f"W{n_of_factors}_{i}" for i in range(1, n_of_factors+1)]
    
    
    
    # obs_df = adata.obs[factor_columns+[f"leiden_{l_param:.2f}"]].copy()
    #cell_type_summarised_data = obs_df.groupby(f"leiden_{l_param:.2f}").mean()
    # cell_type_summarised_data = cell_type_summarised_data[factor_columns]
    
    # In the original script the above part was used to aggreagate factor loading across clusters
    # Below part was used just top plot aggreatgegion acrross all donor inclusing cohort 1
    obs_df = adata.obs[factor_columns+["condition"]].copy()
    cell_type_summarised_data = obs_df.groupby(f"condition").mean()
    cell_type_summarised_data = cell_type_summarised_data[factor_columns]
    
    
    factor_columns = [f"Factor {i}" for i in range(1, n_of_factors+1)]
    scaler = MinMaxScaler(feature_range=(-1,1))
    scaler.fit(cell_type_summarised_data)
    scaled_factors = scaler.transform(cell_type_summarised_data)
    # sns.clustermap(scaled_factors.T, cmap="viridis", xticklabels=cell_type_summarised_data.index.values, yticklabels=factor_columns)
    # plt.savefig(f'{PLOT_PATH}/{sample_type}_{n_of_factors}_factors_heatmap_celltype.pdf');

    ct_df = pd.DataFrame(scaled_factors.T, index=factor_columns,   columns=cell_type_summarised_data.index.values)
    
    factor_columns = [f"W{n_of_factors}_{i}" for i in range(1, n_of_factors+1)]
    obs_df = adata.obs[factor_columns+["condition"]].copy()
    cond_summarised_data = obs_df.groupby("condition").mean()
    # cond_summarised_data = adata.obs.groupby(["condition"]).mean()
    cond_summarised_data = cond_summarised_data[factor_columns]
    factor_columns = [f"Factor {i}" for i in range(1, n_of_factors+1)]
    scaler = MinMaxScaler(feature_range=(-1,1))
    scaler.fit(cond_summarised_data)
    scaled_factors = scaler.transform(cond_summarised_data)
    # sns.clustermap(scaled_factors.T, cmap="viridis", xticklabels=cond_summarised_data.index.values, yticklabels=factor_columns)
    # plt.savefig(f'{PLOT_PATH}/{sample_type}_{n_of_factors}_factors_heatmap_sample.pdf');
    cond_df = pd.DataFrame(scaled_factors.T, index=factor_columns,   columns=cond_summarised_data.index.values)

    merged_df = pd.concat([ct_df, cond_df], axis=1) # pd.merge(ct_df, cell_type_summarised_data, right_index=True) 
    print()
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


adata_cc_merged = get_merged_adata()
print(adata_cc_merged)
adata_cc_merged = apply_nmf_on_merged_data(adata_cc_merged)

# TODO: Seperate plots
plt.rcParams['figure.figsize'] = [10, 10]
df_rmse = pd.read_csv(os.path.join(output_path, f"{analysis_name}_RMSE.csv"), index_col=0).plot.line()
plt.xticks(rotation=90)
plt.savefig(os.path.join(PLOT_PATH, f"{analysis_name}_RMSE.pdf"))
print(adata_cc_merged)
analyse_nmf_results(adata_cc_merged, 15,  42)
"""
for col in adata_cc_merged.var.columns:
    if col.lower().startswith("mt-") or col.lower().startswith("rp-"):
        adata_cc_merged.var[col] = adata_cc_merged.var[col].astype(str)
adata_cc_merged.write(os.path.join(OUT_DATA_PATH, f'{sample_type}_merged_nnmf.h5ad'))"""

"""adata_cc_merged = sc.read_h5ad(os.path.join(OUT_DATA_PATH, f'{sample_type}_merged_nnmf.h5ad'))
analyse_nmf_results(42)
adata_integ_clust = sc.read_h5ad(os.path.join(OUT_DATA_PATH, f'{sample_type}_integrated_clustered.h5ad'))"""

# WARNING: If you want to generate the aggregation across clusters use the script below. also check generateheatmap function.
# adata_cc_merged  = adata_cc_merged[adata_cc_merged.obs.condition.isin(['CD-AOM-DSS', 'CD-no-AOM-DSS', 'HFD-AOM-DSS', 'HFD-no-AOM-DSS'])]
# adata_cc_merged.obs[f"leiden_{l_param:.2f}"] = adata_integ_clust.obs[f"leiden_{l_param:.2f}"].values
# adata_cc_merged.write(os.path.join(OUT_DATA_PATH, f'{sample_type}_merged_nnmf.h5ad'))
# generate_heatmap(adata_cc_merged, n_of_factors=20)


#  python vis_nnmf_deconv_joint.py -i ../data/out_data/ -o  ../data/out_data -an visium_nnmf_deconv_joint
