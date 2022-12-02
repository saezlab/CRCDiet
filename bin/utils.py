import os
import pickle
import random
import pandas as pd
import scanpy as sc
import numpy as np
from tqdm import tqdm
from pathlib import Path
from IPython.display import Markdown, display
from collections import Counter
sc.settings.verbosity = 0


S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
SC_RAW_DATA_PATH = os.path.join(S_PATH, "../data", "sc_outputs")
VIS_RAW_DATA_PATH = os.path.join(S_PATH, "../data", "visium_outputs")
OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
PLOT_PATH =  os.path.join(S_PATH, "plots")



def read_pickle(fl_path):
    """Read a pickle file

    This function reads a pickle file and returns the object.

    Args:
        fl_path (str): pickle file to be read


    Returns:
        object: return the object stored in the pickle file.

    """

    p_file = open(fl_path,'rb')
    obj_r = pickle.load(p_file)
    p_file.close()

    return obj_r

def write_pickle(fl_path, obj_w):
    """Create a pickle file given object

    This function creates a pickle file to store the input object in the given path.

    Args:
        fl_path (str): pickle file to be created
        obj_w (object): object to be stored

    """
    p_file = open(fl_path, 'wb')    
    # dump information to that file
    pickle.dump(obj_w, p_file)
    p_file.close()


def get_meta_data(sc_or_vis):
    """
    Read  meta data file and return the dataframe

    Returns:
        Meta data as data frame 
    """

    meta = pd.read_csv(os.path.join(DATA_PATH, f'{sc_or_vis}_meta_data.csv'))
    return meta


def read_raw_sc_sample(sample_name):
    '''Read raw single-cell data'''

    adata = sc.read_10x_mtx(os.path.join(SC_RAW_DATA_PATH, sample_name, "outs", 'filtered_feature_bc_matrix'), cache=True)
    adata.var_names_make_unique()
    return adata

def read_raw_visium_sample(sample_name):
    '''Read raw visium data'''
    adata = sc.read_visium(os.path.join(VIS_RAW_DATA_PATH,sample_name , "outs"), count_file="filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    return adata


def printmd(string):
    display(Markdown(string))

def get_filtered_concat_data(sample_type):
    meta = get_meta_data(sample_type)
    # TODO: Use whole transcriptome instead of HVGs
    # Comment out the section below for running DEGs on HVGs
    # TODO: Refactor this script, it is ugly and inefficient
    adata_concat = []
    for ind, row in tqdm(meta.iterrows(), total=meta.shape[0]):
        
        
        sample_id = row["sample_id"]
        condition = row["condition"]
        
        print(f"Merging {sample_id}...")

        tmp = sc.read_h5ad(os.path.join(OUT_DATA_PATH,f"{sample_id}_filtered.h5ad"))
        # Fetch sample metadata
        m = meta[meta['sample_id'] == sample_id]
        # Add metadata to adata
        for col in m.columns:
            tmp.obs[col] = m[col].values[0]
        # Append
        adata_concat.append(tmp)
        del tmp
    
    # Merge objects and delete list
    adata_concat = adata_concat[0].concatenate(adata_concat[1:], join='outer')
    return adata_concat


def get_unfiltered_concat_data(sample_type):

    adata_concat = []
    meta = get_meta_data(sample_type)
    # TODO: Use whole transcriptome instead of HVGs
    # Comment out the section below for running DEGs on HVGs
    # TODO: Refactor this script, it is ugly and inefficient
    for _, row in meta.iterrows():
    
        sample_id = row["sample_id"]
        condition = row["condition"]
        print(f"Merging {sample_id}...")

        tmp = utils.read_raw_sc_sample(sample_id)
        # Fetch sample metadata
        m = meta[meta['sample_id'] == sample_id]
        # Add metadata to adata
        for col in m.columns:
            tmp.obs[col] = m[col].values[0]
        # Append
        adata_concat.append(tmp)
        del tmp

    # Merge objects and delete list
    adata_concat = adata_concat[0].concatenate(adata_concat[1:], join='outer')
    return adata_concat


def set_n_return_paths(name):
    """
    This function sets the paths
    """
    S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
    DATA_PATH = os.path.join(S_PATH, "../data")
    OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
    PLOT_PATH =  os.path.join(S_PATH, "../plots", name)

    Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
    Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = PLOT_PATH

    return S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH

def read_mtx_file(mtx_file_path, var_name_list= None, obs_name_list= None, transpose=False):
    adata = sc.read_mtx(mtx_file_path)
    if transpose:
        adata = adata.transpose()
    adata.var_names = var_name_list
    adata.obs_names = obs_name_list
    return adata


def calculate_proportions_from_list(lst_cell_type, cell_count_lst, c_type_list=None):
    
    sample_1_list = lst_cell_type[:cell_count_lst[0]]
    sample_2_list = lst_cell_type[cell_count_lst[0]:]
    c1 = Counter(sample_1_list) 
    c2 = Counter(sample_2_list) 
    proportion_s1 = []
    proportion_s2 = []



    for c_type in c_type_list:
        proportion_s1.append(100*(c1[c_type]/cell_count_lst[0]))
        proportion_s2.append(100*(c2[c_type]/cell_count_lst[1]))

    return proportion_s1, proportion_s2


def calculate_cell_type_proportions(sample_list, adata=None, obs_col = "cell_type_0.20", sample_type="sc"):
    input_path = "../data/out_data/sc_integrated_cluster_scannot.h5ad"
    sample_type="sc"
    adata = sc.read_h5ad(input_path)
    meta = get_meta_data(sample_type)
    condition = list(np.unique(meta['condition']))
    c_type_list = list(adata.obs[obs_col].cat.categories)
    cell_count_lst = []
    sample_list = sample_list.split(",")
    all_cells_cell_type_list = []
    samp_prop_dict = dict()
    samp_propor_arr = []
    for samp in sample_list:
        samp_prop_dict[samp] = dict()
        samp_propor_arr.append([])

        adata_tmp = adata[adata.obs["condition"]==samp,:]
        all_cells_cell_type_list.extend(list(adata_tmp.obs[obs_col]))
        cell_count_lst.append(adata_tmp.shape[0])
        for c_type in c_type_list:
            # print("c_type", c_type, adata_tmp[adata_tmp.obs[obs_col]==c_type].shape)
            proportion = 100*(adata_tmp[adata_tmp.obs[obs_col]==c_type].shape[0]/adata_tmp.shape[0])
            samp_prop_dict[samp][c_type] = proportion
            samp_propor_arr[-1].append(proportion)
    
    return c_type_list, samp_propor_arr[0], samp_propor_arr[1], samp_prop_dict, all_cells_cell_type_list, cell_count_lst

def random_populations(number_of_permutations):
    c_type_list, samp_propor_1, samp_propor_2, samp_prop_dict, all_cells_cell_type_list, cell_count_lst = calculate_cell_type_proportions("CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN", adata=None, obs_col = "cell_type_0.20", sample_type="sc")

    all_diffs_simulations = []
    cpy_all_cells_cell_type_list = all_cells_cell_type_list[:]
    for i in range(number_of_permutations):
        random.shuffle(cpy_all_cells_cell_type_list)
        rand_proportion_s1, rand_proportion_s2 = calculate_proportions_from_list(cpy_all_cells_cell_type_list, cell_count_lst, c_type_list)
        diff_proportion = np.array(rand_proportion_s1) - np.array(rand_proportion_s2) 
        all_diffs_simulations.append(diff_proportion)
    print(all_diffs_simulations)
random_populations(10)

