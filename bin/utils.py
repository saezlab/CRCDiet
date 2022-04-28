import os
import pickle
import pandas as pd
import scanpy as sc
from IPython.display import Markdown, display

sc.settings.verbosity = 0


S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
DATA_PATH = os.path.join(S_PATH, "../data")
SC_RAW_DATA_PATH = os.path.join(S_PATH, "../data", "sc_outputs")
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
    adata = sc.read_10x_mtx(os.path.join(SC_RAW_DATA_PATH, sample_name, "outs", 'filtered_feature_bc_matrix'), cache=True)
    adata.var_names_make_unique()
    return adata

def read_raw_visium_sample(sample_name):
    adata = sc.read_10x_mtx(os.path.join(DATA_PATH, sample_name, 'filtered_feature_bc_matrix'), cache=True)
    adata.var_names_make_unique()
    return adata


def printmd(string):
    display(Markdown(string))