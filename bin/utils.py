import os
import pickle
import pandas as pd
import scanpy as sc
from tqdm import tqdm
from pathlib import Path
from IPython.display import Markdown, display

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


