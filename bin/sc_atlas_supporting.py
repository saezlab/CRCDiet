import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import utils


def transfer_cell_type_annotation(analysis_name):
    S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name)
    plt.rcParams['figure.dpi']= 300
    plt.rcParams['figure.figsize']= (45, 30)

    adata = sc.read_h5ad("../data/out_data/atlas_integrated_clustered.h5ad")
    adata_n = sc.read_h5ad("../data/out_data/sc_integrated_cluster_scannot.h5ad")


    adata.obs["cell_type"] = "no-annot"

    crcdiet_samples = ["CD-AOM-DSS-Epi_plus_DN", "CD-AOM-DSS-Immune", "HFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Immune", "LFD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Immune"]
    adata.obs.loc[adata.obs.condition.isin(crcdiet_samples), "cell_type"] = adata_n.obs["cell_type_0.20"].values

    """sc.pl.umap(
        adata, color="cell_type",
        title= "Cell Type", legend_loc='on data', size=10, frameon=False, show=True, save=f"atlas_crcdiet_cell_type_harmony"
    )
    sc.pl.umap(
        adata, color="cell_type",
        title= "Cell Type", size=10, frameon=False, show=True, save=f"atlas_crcdiet_cell_type_on_data_harmony"
    )"""
    
    haber_samples = ["Haber_GSE92332_GSE92332_SalmonellaInfect", "Haber_GSE92332_GSE92332_SalmHelm", "Haber_GSE92332_GSE92332_atlas", "Haber_GSE92332_FAE"]

    annot_lst = []
    for val  in adata[adata.obs.condition.isin(haber_samples),:].obs_names:
        annot_lst.append(val.split("_")[-1].split("-")[0])

    adata.obs.loc[adata.obs.condition.isin(haber_samples), "cell_type"] = annot_lst

    

    sc.pl.umap(
        adata, color="cell_type",
        title= "Cell Type", legend_loc='on data', size=10, frameon=False, show=True, save=f"atlas_habercrc_all_cell_type_on_data_harmony"
    )
    sc.pl.umap(
        adata, color="cell_type",
        title= "Cell Type", size=10, frameon=False, show=True, save=f"atlas_habercrc_all_cell_type_harmony"
    )

# transfer_cell_type_annotation("atlas_cluster")

def create_subanndata():
    atlas_integrated = sc.read_h5ad("../data/out_data/atlas_integrated_clustered.h5ad")
    atlas_cell_type_annot = sc.read_h5ad("../data/out_data/sc_integrated_cluster_scannot.h5ad")
    atlas_integrated.obs["cell_type"] = atlas_cell_type_annot.obs["cell_type_level1"]
    anndata_b_cells = atlas_integrated[atlas_integrated.obs["cell_type"].isin(['B cells-1', 'B cells-2', 'B cells-3', 'B cells-5'])]
    anndata_bcell_neutrophils_iga =  atlas_integrated[atlas_integrated.obs["cell_type"].isin(['B cells-1', 'Neutrophils', 'IgA plasma cells-1']),:]