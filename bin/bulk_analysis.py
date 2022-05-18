import os
import math
import utils
import numpy as np
import pandas as pd
import scanpy as sc
import qc_preprocess
import seaborn as sns
from anndata import AnnData
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from sklearn import preprocessing
from scipy.spatial import distance
from scipy.cluster import hierarchy
from utils import out_data_path, plot_path, data_path

sc.settings.verbosity = 3
np.random.seed(123)

meta = utils.get_meta_data()


def get_pseudobulk_sample(adata):
    """
    Given sample id, create pseudobulk sample
    Args:
        sample_id (int): Sample id in meta_data.csv file
    
    Returns:
        numpy array: pseudobulk array for the given sample
    """


    # Different normalization techniques will also be tried. 
    # sc.pp.normalize_total(adata, inplace=True)
    # sc.pp.log1p(adata)

    pseudo_bulk = list(np.squeeze(np.asarray(np.sum(adata.X, axis=0))))
    return pseudo_bulk



def get_pseoudobulk_samples_all():
    """
    This function creates pseudobulk profile for all samples
    
    Returns:
        list: list of pseudobulk samples
    """

    pseudobulk_lst = []

    for ind, row in meta.iterrows():
        sample_id = row["sample_id"]
        condition = row["condition"]

        if sample_id not in ["sample04"]:
            pseudo_bulk_samp = get_pseudobulk_sample(sample_id)
            pseudobulk_lst.append(pseudo_bulk_samp)
    
    return pseudobulk_lst

def get_pseoudobulk_samples():
    pseudobulk_lst = []
    merged_adata = qc_preprocess.get_concatenated_adata_files()

    for ind, row in meta.iterrows():
        sample_id = row["sample_id"]
        condition = row["condition"]

        if sample_id not in ["sample04"]:
            print(sample_id)
            adata_samp = merged_adata[merged_adata.obs_names.str.contains(sample_id),:]
            pseudo_bulk_samp = get_pseudobulk_sample(adata_samp)
            pseudobulk_lst.append(pseudo_bulk_samp)
    return pseudobulk_lst

def apply_MDS_pseudocount_plot_result():
    """
    This function creates pseudobulk profile for all samples
    
    """
    pseudobulk_lst = get_pseoudobulk_samples()
    ps_counts_adata = AnnData(np.array(pseudobulk_lst))
    
    X_norm = sc.pp.normalize_total(ps_counts_adata, target_sum=1e6, inplace=False)['X']
    embedding = MDS(n_components=2)
    X_transformed = embedding.fit_transform(X_norm)
    plt.figure(dpi=300)
    plt.plot(X_transformed[:, 0], X_transformed[:, 1], '.', markersize=20)
    lbl_list = list(meta["label"])
    lbl_list.remove('Control: Day 70')
    print(lbl_list)
    for x, y, text in zip(X_transformed[:, 0], X_transformed[:, 1], lbl_list):
        plt.text(x, y, text, fontsize=6)

    plt.yticks(fontsize=6)
    plt.xticks(fontsize=6) #, rotation=90)
    plt.xlabel("MDS Dimension - 1")
    plt.ylabel("MDS Dimension - 2")
    plt.savefig(os.path.join(plot_path, "preprocess", "bulk_mds_result.png"))
    



def calculate_dist_matrix_PS_divergence():
    pseudobulk_lst = get_pseoudobulk_samples()

    ps_counts_adata = AnnData(np.array(pseudobulk_lst))
    X_norm = sc.pp.normalize_total(ps_counts_adata, target_sum=1.0, inplace=False)['X']
    # X_norm = ps_counts_adata.X
    # X_norm = sc.pp.log1p(sc.pp.normalize_total(ps_counts_adata, target_sum=1.0, inplace=False))['X']

    print(X_norm.shape)
    dist_arr = np.zeros((5,5))
    # print(dist_arr)

    for samp_id in range(1,7):

        if samp_id!=4:
            samp_id = samp_id if samp_id<4 else samp_id-1
            row = meta[meta["sample_id"]==f"sample0{samp_id}"]
            condition = str(row["condition"].values[0])
            dist1 = X_norm[samp_id-1,:]

        for samp_id2 in range(1,7):
            if samp_id2!=4:
                samp_id2 = samp_id2 if samp_id2<4 else samp_id2-1
                row2 = meta[meta["sample_id"]==f"sample0{samp_id2}"]
                condition2 = str(row2["condition"].values[0])
                dist2 = X_norm[samp_id2-1,:]
                dist = math.sqrt(distance.jensenshannon(dist1, dist2))
                # print(samp_id, samp_id2, dist)
                dist_arr[samp_id-1,samp_id2-1] = dist

    return dist_arr

def apply_hierarchical_clustering_plot_dendogram():

    pseudobulk_lst = get_pseoudobulk_samples()


    dist_arr = calculate_dist_matrix_PS_divergence()
    
    lbl_list = list(meta["label"])
    lbl_list.remove('Control: Day 70')


    frame = pd.DataFrame(np.array(dist_arr), columns=lbl_list, index=lbl_list)
    Z = hierarchy.linkage(frame, 'complete')
    plt.figure(dpi=300)
    dendro = hierarchy.dendrogram(Z, leaf_rotation=45, labels=frame.index)
    plt.subplots_adjust(bottom=0.30)
    plt.yticks(fontsize=6)
    plt.xticks(fontsize=6) #, rotation=90)
    plt.title('Hierarchical Clustering of Samples')
    plt.ylabel("Jensen–Shannon Divergence", fontsize=6)
    plt.savefig(os.path.join(plot_path, "preprocess", "bulk_hierc_clust.png"))




"""
def normalize_sample(sample_num):
    sample_path = f"{path}/sample0{sample_num}/outs"
    print(sample_path)
    adata = sc.read_visium(path=sample_path, count_file="filtered_feature_bc_matrix.h5" )
    adata.var_names_make_unique()
    print("samp_num", adata.X.shape)
    x_arr = adata.X
    sum_of_rows = x_arr.sum(axis=1)
    normalized_array = adata.X / sum_of_rows[:, np.newaxis]

    return normalized_array

"""