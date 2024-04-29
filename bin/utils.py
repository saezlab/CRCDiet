import os
import pickle
import random
import pandas as pd
import scanpy as sc
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
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
    adata=None
    if "Helminth" in sample_name:
        adata = sc.read_visium(os.path.join(VIS_RAW_DATA_PATH,sample_name ), count_file="h_poly_filtered_feature_bc_matrix.h5")
        adata = adata[:,adata.var_names.str.contains("refdata-gex-mm10")]
        adata.var.index = pd.Index(gen.split("refdata-gex-mm10__________")[-1] for gen in adata.var.index.values)
    else:
        adata = sc.read_visium(os.path.join(VIS_RAW_DATA_PATH,sample_name , "outs"), count_file="filtered_feature_bc_matrix.h5")

    adata.var_names_make_unique()
    return adata


def read_filtered_visium_sample(sample_id):
    sample = sc.read_h5ad(os.path.join(OUT_DATA_PATH,f"{sample_id}_filtered.h5ad"))
    return sample

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

        # if "no-AOM-DSS" not in sample_id:
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
        tmp = None
        if sample_type=="sc":
            tmp = read_raw_sc_sample(sample_id)
        elif sample_type=="visium":
            tmp = read_raw_visium_sample(sample_id)
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


def set_n_return_paths(name, plot=True):
    """
    This function sets the paths
    """
    S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
    DATA_PATH = os.path.join(S_PATH, "../data")
    OUT_DATA_PATH = os.path.join(DATA_PATH, "out_data")
    PLOT_PATH =  os.path.join(S_PATH, "../plots", name)

    Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
    if plot:
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
    # input_path = "../data/out_data/sc_integrated_cluster_scannot.h5ad"
    df_cond_cell_type = pd.read_csv("../data/out_data/sc_cond_cell_type.csv", index_col=0)
    # print(df_cond_cell_type)
    # adata = sc.read_h5ad(input_path)
    meta = get_meta_data(sample_type)
    # c_type_list = list(adata.obs[obs_col].cat.categories)
    c_type_list = list(set(df_cond_cell_type["cell_type_0.20"]))
    cell_count_lst = []
    all_cells_cell_type_list = []
    samp_prop_dict = dict()
    samp_cell_count_dict = dict()
    samp_propor_arr = []
    for samp in sample_list:
        samp_prop_dict[samp] = dict()
        samp_cell_count_dict[samp] = dict()
        samp_propor_arr.append([])

        # adata_tmp = adata[adata.obs["condition"]==samp,:]
        df_tmp =  df_cond_cell_type[df_cond_cell_type["condition"]==samp]
        # all_cells_cell_type_list.extend(list(adata_tmp.obs[obs_col]))
        all_cells_cell_type_list.extend(list(df_tmp[obs_col]))
        # cell_count_lst.append(adata_tmp.shape[0])
        cell_count_lst.append(len(df_tmp))
        
        for c_type in c_type_list:
            # c_type_count=adata_tmp[adata_tmp.obs[obs_col]==c_type].shape[0]
            c_type_count=len(df_tmp[df_tmp[obs_col]==c_type])
            # print(f"Sample:\t{samp}\tC_type:\t{c_type}\tCell Type Count:\t{c_type_count}\tTotal count:\t{adata_tmp.shape[0]}")
            samp_cell_count_dict[samp][c_type] = c_type_count
            # proportion = 100*(c_type_count/adata_tmp.shape[0])
            proportion = 100*(c_type_count/len(df_tmp))
            samp_prop_dict[samp][c_type] = proportion
            samp_propor_arr[-1].append(proportion)

        del df_tmp

    # print("sample proportions:", samp_propor_arr)
    all_ct_prop_diff_arr = np.array(samp_propor_arr[0])- np.array(samp_propor_arr[1])

    # del adata

    return c_type_list, all_ct_prop_diff_arr, samp_prop_dict, samp_cell_count_dict, all_cells_cell_type_list, cell_count_lst


# calculate_cell_type_proportions(["HFD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN", "CD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Immune", "LFD-AOM-DSS-Immune", "CD-AOM-DSS-Immune"])

def random_populations(str_sample_list, number_of_simulations):

    sample_list = str_sample_list.split(",")
    c_type_list, all_ct_prop_diff_arr, samp_prop_dict, samp_cell_count_dict, all_cells_cell_type_list, cell_count_lst = calculate_cell_type_proportions(sample_list,  adata=None, obs_col = "cell_type_0.20", sample_type="sc")
    
    all_diffs_simulations = []
    cpy_all_cells_cell_type_list = all_cells_cell_type_list[:]
    for i in range(number_of_simulations):
        random.shuffle(cpy_all_cells_cell_type_list)
        rand_proportion_s1, rand_proportion_s2 = calculate_proportions_from_list(cpy_all_cells_cell_type_list, cell_count_lst, c_type_list)
        diff_proportion = np.array(rand_proportion_s1) - np.array(rand_proportion_s2) 
        all_diffs_simulations.append(diff_proportion)
    all_diffs_simulations = np.array(all_diffs_simulations)
    

    
    
    

    """if "Epi_plus_DN" in sample_list[0]:
        for ct in [ "T cells", "Plasma cells", "Mast cells", "Neutrophils", "ILC2",  "Dendritic cells", "Myeloid cells", "B cells"]:
            # print(["Mast cells", "Dendritic cells", "Neutrophils", "ILC2"], sample_list[0])
            c_type_list.remove(ct)
        
    elif "Immune" in sample_list[0]:
        # print(["Tuft cells", "Goblet cells", "Prolif.", "Enteroendocrine", "Keratynocytes", "Prolif. + Mature enterocytes"], sample_list[0])
        for ct in ["Tuft cells", "Goblet cells", "Prolif.", "Enteroendocrine", "Keratynocytes", "Prolif. + Mature enterocytes"]:
            c_type_list.remove(ct)"""


    dict_cell_type_pval = dict()
    fig, axs = plt.subplots(len(c_type_list), 2, figsize=(15, 90))
    for ind, c_type in enumerate(c_type_list):
        if "Epi_plus_DN" in sample_list[0]:
            if c_type in [ "T cells", "Plasma cells", "Mast cells", "Neutrophils", "ILC2",  "Dendritic cells", "Myeloid cells", "B cells"]:
                # print(["Mast cells", "Dendritic cells", "Neutrophils", "ILC2"], sample_list[0])
                continue
        
        elif "Immune" in sample_list[0]:
            # print(["Tuft cells", "Goblet cells", "Prolif.", "Enteroendocrine", "Keratynocytes", "Prolif. + Mature enterocytes"], sample_list[0])
            if c_type in ["Tuft cells", "Goblet cells", "Prolif.", "Enteroendocrine", "Keratynocytes", "Prolif. + Mature enterocytes"]:
                continue
        samp_prop_list = []
        for samp in sample_list:
            samp_prop_list.append(samp_prop_dict[samp][c_type])
            # print(f"Condition-{samp}\tCell type-{c_type}: {samp_cell_count_dict[samp][c_type]}")
        
        real_diff = all_ct_prop_diff_arr[ind]
        ct_null_dist = all_diffs_simulations[:,ind]
        
        num_of_rare_or_rarer = 0.0
        if real_diff <0.0:
            num_of_rare_or_rarer = np.sum(ct_null_dist<real_diff)
        else:
            num_of_rare_or_rarer = np.sum(ct_null_dist>real_diff)

        p_val = float(num_of_rare_or_rarer)/float(number_of_simulations)
        dict_cell_type_pval[c_type] = p_val
        
        n, bins, patches = axs[ind][0].hist(ct_null_dist, bins=10, facecolor='g')
        axs[ind][0].axvline(real_diff, color='k', linestyle='dashed', linewidth=1)
        axs[ind][0].set_xlabel('Null distribution')
        axs[ind][0].set_ylabel('Frequency')
        axs[ind][0].set_title(f'Histogram of Null Distribution - Cell Type: {c_type}')
        axs[ind][0].grid(True)

        axs[ind][0].text(.79, 0.99, f"p =  {p_val:.2e}", ha='left', va='top', transform=axs[ind][0].transAxes)

        """plt.savefig(f"../plots/sc_cell_type_prop/hist_{str_sample_list}_{c_type}.pdf")
        plt.clf()"""

        axs[ind][1].bar(sample_list, samp_prop_list, color = (0.3,0.1,0.4,0.6))
        axs[ind][1].text(.79, 0.99, f"{p_val:.2e}", ha='left', va='top', transform=axs[ind][1].transAxes)
        # axs[ind][1].text(.05, 0.99, f"{samp_cell_count_dict[sample_list[0]][c_type]}, {cell_count_lst[0]}, {samp_cell_count_dict[sample_list[1]][c_type]}, {cell_count_lst[1]}", ha='right', va='top', transform=axs[ind][1].transAxes)
        axs[ind][1].set_title(f'Cell Type Proportion: {c_type}')
        axs[ind][1].set_xlabel('Condition')
        axs[ind][1].set_ylabel('Proportion (%)')

    
    fig.tight_layout()
    fig.savefig(f"../plots/sc_cell_type_prop/barplot_{str_sample_list}.pdf")
    fig.savefig(f"../plots/sc_cell_type_prop/barplot_{str_sample_list}.png")

    return c_type_list, samp_prop_dict, dict_cell_type_pval
        

def p_to_star(p_val):
    if p_val > 0.05:
        return "n.s."
    elif p_val<=0.05 and p_val>0.01:
        return "*"
    elif p_val<=0.01 and p_val>0.001:
        return "**"
    elif p_val<=0.001:
        return "***"

def calculate_pval(str_sample_list, number_of_simulations):

    sample_list = str_sample_list.split(",")
    c_type_list, all_ct_prop_diff_arr, samp_prop_dict, samp_cell_count_dict, all_cells_cell_type_list, cell_count_lst = calculate_cell_type_proportions(sample_list,  adata=None, obs_col = "cell_type_0.20", sample_type="sc")
    
    all_diffs_simulations = []
    cpy_all_cells_cell_type_list = all_cells_cell_type_list[:]
    for i in range(number_of_simulations):
        random.shuffle(cpy_all_cells_cell_type_list)
        rand_proportion_s1, rand_proportion_s2 = calculate_proportions_from_list(cpy_all_cells_cell_type_list, cell_count_lst, c_type_list)
        diff_proportion = np.array(rand_proportion_s1) - np.array(rand_proportion_s2) 
        all_diffs_simulations.append(diff_proportion)
    all_diffs_simulations = np.array(all_diffs_simulations)


    dict_cell_type_pval = dict()
    for ind, c_type in enumerate(c_type_list):
        if "Epi_plus_DN" in sample_list[0]:
            if c_type in [ "T cells", "Plasma cells", "Mast cells", "Neutrophils", "ILC2",  "Dendritic cells", "Myeloid cells", "B cells"]:
                # print(["Mast cells", "Dendritic cells", "Neutrophils", "ILC2"], sample_list[0])
                continue
        
        elif "Immune" in sample_list[0]:
            # print(["Tuft cells", "Goblet cells", "Prolif.", "Enteroendocrine", "Keratynocytes", "Prolif. + Mature enterocytes"], sample_list[0])
            if c_type in ["Tuft cells", "Goblet cells", "Prolif.", "Enteroendocrine", "Keratynocytes", "Prolif. + Mature enterocytes"]:
                continue
        samp_prop_list = []
        for samp in sample_list:
            samp_prop_list.append(samp_prop_dict[samp][c_type])
            # print(f"Condition-{samp}\tCell type-{c_type}: {samp_cell_count_dict[samp][c_type]}")
        
        real_diff = all_ct_prop_diff_arr[ind]
        ct_null_dist = all_diffs_simulations[:,ind]
        
        num_of_rare_or_rarer = 0.0
        if real_diff <0.0:
            num_of_rare_or_rarer = np.sum(ct_null_dist<real_diff)
        else:
            num_of_rare_or_rarer = np.sum(ct_null_dist>real_diff)

        p_val = float(num_of_rare_or_rarer)/float(number_of_simulations)
        dict_cell_type_pval[c_type] = p_val
        

    return c_type_list, samp_prop_dict, dict_cell_type_pval

def obs_to_csv(adata):
    a = pd.concat([df_cond, df_cell_type], axis=1)
    df_cond = adata.obs["condition"]
    df_cell_type = adata.obs["cell_type_0.20"]
    final_df = pd.concat([df_cond, df_cell_type], axis=1)


    final = pd.read_csv("../data/out_data/sc_cond_cell_type.csv", index_col=0)


def subset_ann_data(adata_pth, obs_col, group_lst):
    adata = sc.read_h5ad(adata_pth)
    group_lst = group_lst.split(",")
    return adata[adata.obs[obs_col].isin(group_lst),:]


def create_ena_read_submission_file(data_path="/Volumes/Transcend/", study="PRJEB64135"):
    lst_samples_path = [("ERS16007550", "ext_L24854_CD-AOM-DSS-Immune-single-cell_S7"), ("ERS16007552", "ext_L24854_HFD-AOM-DSS-Immune-single-cell_S9"), 
                        ("ERS16007551", "ext_L24854_LFD-AOM-DSS-Immune-single-cell_S8"), ("ERS15990591", "ext_L24854_CD-AOM-DSS-Epi_plus_DN-single-cell_S10"), 
                        ("ERS15990593", "ext_L24854_HFD-AOM-DSS-Epi_plus_DN-single-cell_S12"), ("ERS15990592", "ext_L24854_LFD-AOM-DSS-Epi_plus_DN-single-cell_S11"), 
                        ("ERS15990585", "ext_L24854_CD-no-AOM-DSS-colon-d81-visium_S1"), ("ERS15990586", "ext_L24854_LFD-no-AOM-DSS-colon-d81-visium_S5"), 
                        ("ERS15990587", "ext_L24854_HFD-no-AOM-DSS-colon-d81-visium_S6"), ("ERS15990588", "ext_L24854_CD-AOM-DSS-colon-d81-visium_S2"), 
                        ("ERS15990589", "ext_L24854_LFD-AOM-DSS-colon-d81-visium_S3"), ("ERS15990590", "ext_L24854_HFD-AOM-DSS-colon-d81-visium_S4")]
    for samp, folder in lst_samples_path:
        path = os.path.join(data_path, folder)
        for fl in os.listdir(path):
            if fl.endswith(".fastq.gz"):
                if "_I1_" in fl:
                    fl_2 = fl.replace("_I1_", "_I2_")
                    str_acc = [samp, study, "Illumina NovaSeq 6000", "", "TRANSCRIPTOMIC", "cDNA",	"RNA-Seq",	"PAIRED", fl, "", fl_2, ""]
                    print("\t".join(str_acc))
                elif "_R1_" in fl:
                    fl_2 = fl.replace("_R1_", "_R2_")
                    str_acc = [samp, study, "Illumina NovaSeq 6000", "", "TRANSCRIPTOMIC", "cDNA",	"RNA-Seq",	"PAIRED", fl, "", fl_2, ""]
                    print("\t".join(str_acc))
                else:
                    pass
                
# create_ena_read_submission_file()


"""sample_pair_list = ["CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN", "CD-AOM-DSS-Epi_plus_DN,HFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN",
                    "CD-AOM-DSS-Immune,LFD-AOM-DSS-Immune", "CD-AOM-DSS-Immune,HFD-AOM-DSS-Immune", "HFD-AOM-DSS-Immune,LFD-AOM-DSS-Immune"]
p_val_list= []

for samp_pair in tqdm(sample_pair_list):
    for i in tqdm(range(1000)):
        
        c_type_list12, samp_prop_dict12, dict_cell_type_pval12 = calculate_pval(samp_pair, 10000)
        for c_type in dict_cell_type_pval12.keys():
            print([samp_pair, c_type, dict_cell_type_pval12[c_type]])
            p_val_list.append([samp_pair, c_type, dict_cell_type_pval12[c_type]])"""
    


"""from plotting import plot_significance

c_type_list12, samp_prop_dict12, dict_cell_type_pval12 = random_populations("CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN", 10000)
c_type_list13, samp_prop_dict13, dict_cell_type_pval13 = random_populations("CD-AOM-DSS-Epi_plus_DN,HFD-AOM-DSS-Epi_plus_DN", 10000)
c_type_list23, samp_prop_dict23, dict_cell_type_pval23 = random_populations("HFD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN", 10000)"""
"""print(c_type_list12)
print("samp_prop_dict12", samp_prop_dict12)
print("samp_prop_dict23", samp_prop_dict23)
print("samp_prop_dict13", samp_prop_dict13)"""
"""for ind, c_type in enumerate(c_type_list12):
    # plot_significance("CD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN", samp_prop_dict12["CD-AOM-DSS-Epi_plus_DN"][c_type_list12[ind]], samp_prop_dict12["LFD-AOM-DSS-Epi_plus_DN"][c_type_list12[ind]], samp_prop_dict23["HFD-AOM-DSS-Epi_plus_DN"][c_type_list12[ind]], f"{dict_cell_type_pval12[c_type_list12[ind]]:.2e}", f"{dict_cell_type_pval23[c_type_list12[ind]]:.2e}", f"{dict_cell_type_pval13[c_type_list12[ind]]:.2e}")
    print(c_type)
    try:
        plot_significance("CD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN", samp_prop_dict12["CD-AOM-DSS-Epi_plus_DN"][c_type_list12[ind]], samp_prop_dict12["LFD-AOM-DSS-Epi_plus_DN"][c_type_list12[ind]], samp_prop_dict23["HFD-AOM-DSS-Epi_plus_DN"][c_type_list12[ind]], p_to_star(dict_cell_type_pval12[c_type_list12[ind]]), p_to_star(dict_cell_type_pval23[c_type_list12[ind]]), p_to_star(dict_cell_type_pval13[c_type_list12[ind]]), c_type, "Epi_plus_DN")
    except:
        pass


c_type_list12, samp_prop_dict12, dict_cell_type_pval12 = random_populations("CD-AOM-DSS-Immune,LFD-AOM-DSS-Immune", 10000)
c_type_list13, samp_prop_dict13, dict_cell_type_pval13 = random_populations("CD-AOM-DSS-Immune,HFD-AOM-DSS-Immune", 10000)
c_type_list23, samp_prop_dict23, dict_cell_type_pval23 = random_populations("HFD-AOM-DSS-Immune,LFD-AOM-DSS-Immune", 10000)
for ind, c_type in enumerate(c_type_list12):
    # plot_significance("CD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN", samp_prop_dict12["CD-AOM-DSS-Epi_plus_DN"][c_type_list12[ind]], samp_prop_dict12["LFD-AOM-DSS-Epi_plus_DN"][c_type_list12[ind]], samp_prop_dict23["HFD-AOM-DSS-Epi_plus_DN"][c_type_list12[ind]], f"{dict_cell_type_pval12[c_type_list12[ind]]:.2e}", f"{dict_cell_type_pval23[c_type_list12[ind]]:.2e}", f"{dict_cell_type_pval13[c_type_list12[ind]]:.2e}")
    print(c_type)
    try:
        plot_significance("CD-AOM-DSS-Immune", "LFD-AOM-DSS-Immune", "HFD-AOM-DSS-Immune", samp_prop_dict12["CD-AOM-DSS-Immune"][c_type_list12[ind]], samp_prop_dict12["LFD-AOM-DSS-Immune"][c_type_list12[ind]], samp_prop_dict23["HFD-AOM-DSS-Immune"][c_type_list12[ind]], p_to_star(dict_cell_type_pval12[c_type_list12[ind]]), p_to_star(dict_cell_type_pval23[c_type_list12[ind]]), p_to_star(dict_cell_type_pval13[c_type_list12[ind]]), c_type, "Immune")
    except:
        pass
"""


def get_subset_adata(adata, cell_types="epi", sample_type="sc"):

    epi_samples = ["CD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN"]
    immune_samples = ["CD-AOM-DSS-Immune", "LFD-AOM-DSS-Immune", "HFD-AOM-DSS-Immune"]

    samples = []
    if cell_types=="epi" or cell_types=="stroma":
        samples = epi_samples
    else:
        samples =immune_samples
    
    
    if cell_types=="epi" or cell_types=="stroma":
        samples = epi_samples
    else:
        samples = immune_samples
    
    # input_path = "../data/out_data/sc_integrated_cluster_scannot.h5ad"
    immune_cell_types = [ "T cells", "Plasma cells", "Mast cells", "Neutrophils", "ILC2",  "Dendritic cells", "Myeloid cells", "B cells"]
    epithelial_cell_types = ["Tuft cells", "Goblet cells", "Prolif.", "Enteroendocrine", "Keratynocytes", "Prolif. + Mature enterocytes"]
    stroma_cell_types = ["Stroma", "Myofibroblasts", "Endothelial cells"]
    sample_type="sc"
    # adata = sc.read_h5ad(input_path)
    if cell_types=="epi":
        adata = adata[adata.obs["condition"].isin(samples),:]
        adata=adata[adata.obs["cell_type_0.20"].isin(epithelial_cell_types),:]
        
    elif cell_types=="stroma":
        adata = adata[adata.obs["condition"].isin(samples),:]
        adata=adata[adata.obs["cell_type_0.20"].isin(stroma_cell_types),:]
        # print("stroma", adata)
        # for cat in stroma_cell_types:
            # print(cat)
            # print(adata[adata.obs["cell_type_0.20"]==cat,:  ])
    elif cell_types=="immune":
        adata = adata[adata.obs["condition"].isin(samples),:]
        adata=adata[adata.obs["cell_type_0.20"].isin(immune_cell_types),:]

    return adata


# Create proportion tables for each list
def create_proportion_table(data):
    categories = set(data)
    freq_table = {}
    for category in categories:
        freq_table[category] = data.count(category)/len(data)
        # print(category, data.count(category), len(data))
    return freq_table

def cal_significance_proportions(list1, list2, num_permutations=10000):
    freq_table1 = create_proportion_table(list1)
    freq_table2 = create_proportion_table(list2)
    # print(freq_table1)
    # Create a list of all unique categories
    unique_categories = set(freq_table1.keys()).union(set(freq_table2.keys()))

    # Initialize a dictionary to store p-values for each category
    p_values = {}


    # Perform the permutation test for each category individually
    for category in unique_categories:
        # print(category, freq_table1[category], freq_table2[category])
        observed_diff = freq_table1.get(category, 0) - freq_table2.get(category, 0)
        # print(observed_diff)

        all_data = list1 + list2
        num_obs = len(list1)
        num_perm = 0
        num_of_rare_or_rarer = 0.0
        list_of_means = []
        for _ in range(num_permutations):
            np.random.shuffle(all_data)
            # print(all_data)
            perm_list1 = all_data[:num_obs]
            perm_list2 = all_data[num_obs:]

            perm_diff = perm_list1.count(category)/len(perm_list1) - perm_list2.count(category)/len(perm_list2)
            list_of_means.append(perm_diff)
            """print("First proportion, second p", perm_list1.count(category)/len(perm_list1), perm_list2.count(category)/len(perm_list2))
            print("Cat, perm_diff, obs diff", category, perm_diff, observed_diff)
            print("observed_diff", observed_diff)
            print("Perm. diff", perm_diff)"""
            """if observed_diff <0.0 and perm_diff<observed_diff:

                num_of_rare_or_rarer += 1
                # print("First ", category,  observed_diff, perm_diff)
            elif observed_diff >=0.0 and perm_diff>observed_diff:
                num_of_rare_or_rarer +=1
                # print("Second ", category,  observed_diff, perm_diff)"""
        list_of_means = np.array(list_of_means)
        num_of_rare_or_rarer2 = 0.0
        
        if observed_diff <0.0:
            num_of_rare_or_rarer_l = np.sum(list_of_means<observed_diff)
            num_of_rare_or_rarer_r = np.sum(list_of_means>np.abs(observed_diff))
        else:
            num_of_rare_or_rarer_l = np.sum(list_of_means<-1*observed_diff)
            num_of_rare_or_rarer_r = np.sum(list_of_means>np.abs(observed_diff))
        
        num_of_rare_or_rarer = num_of_rare_or_rarer_l + num_of_rare_or_rarer_r
        # p_val = float(num_of_rare_or_rarer2)/float(num_permutations)
        # print("p_val", category, p_val)
        # plt.hist(list_of_means[100:])
        # plt.show()
        # print(num_of_rare_or_rarer)
        p_value = (num_of_rare_or_rarer + 1) / (num_permutations + 1)
        p_values[category] = p_value

    # Correct p-values using the Benjamini-Hochberg procedure
    p_values_sorted = {k: v for k, v in sorted(p_values.items(), key=lambda item: item[1])}
    q_values = {}
    m = len(p_values_sorted)
    for i, (category, p) in enumerate(p_values_sorted.items()):
        q_values[category] = (p * m) / (i + 1)

    # Set a significance level (e.g., 0.01)
    significance_level = 0.01

    # Report whether the change in proportion is statistically significant for each category
    for category, q in q_values.items():
        if q < significance_level:
            print(f"Category: {category} - q-value: {q:.4f} (Significant change)")
        else:
            print(f"Category: {category} - q-value: {q:.4f} (No significant change)")

    # Report whether the change in proportion is statistically significant for each category
    for category, p in p_values.items():
        if p < significance_level:
            print(f"Category: {category} - p-value: {p} (Significant change)")
        else:
            print(f"Category: {category} - p-value: {p} (No significant change)")



def get_proportion_significance(cell_types="epi", conditions="CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN"):

    adata = get_subset_adata(cell_types, sample_type="sc")
    lst_conditions = conditions.split(",")
    lst_cell_types_1 = list(adata[adata.obs["condition"]==lst_conditions[0],:].obs["cell_type_0.20"].values)
    lst_cell_types_2 = list(adata[adata.obs["condition"]==lst_conditions[1],:].obs["cell_type_0.20"].values)
    freq_table1 = create_proportion_table(lst_cell_types_1)
    freq_table2 = create_proportion_table(lst_cell_types_2)
    print(lst_conditions[0], freq_table1)
    print(lst_conditions[1], freq_table2)
    # cal_significance_proportions(lst_cell_types_1, lst_cell_types_2, num_permutations=10000)

# get_proportion_significance()
# get_proportion_significance(cell_types="immune", conditions="CD-AOM-DSS-Immune,LFD-AOM-DSS-Immune")
# get_proportion_significance(cell_types="immune", conditions="CD-AOM-DSS-Immune,HFD-AOM-DSS-Immune")"""


sample_pair_list = ["CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN", "CD-AOM-DSS-Epi_plus_DN,HFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN",
                    "CD-AOM-DSS-Immune,LFD-AOM-DSS-Immune", "CD-AOM-DSS-Immune,HFD-AOM-DSS-Immune", "HFD-AOM-DSS-Immune,LFD-AOM-DSS-Immune"]
