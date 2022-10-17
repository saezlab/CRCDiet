import pandas as pd
import scanpy as sc


def transfer_cell_type_annotation():
    adata = sc.read_h5ad("../data/out_data/atlas_integrated.h5ad")
    adata_n = sc.read_h5ad("../data/out_data/sc_integrated_cluster_scannot.h5ad")


    adata.obs["cell_type"] = "no-annot"

    crcdiet_samples = ["CD-AOM-DSS-Epi_plus_DN", "CD-AOM-DSS-Immune", "HFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Immune", "LFD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Immune"]
    adata.obs.loc[adata.obs.condition.isin(crcdiet_samples), "cell_type"] = adata_n.obs["cell_type_0.20"].values

    """haber_samples = ["Haber_GSE92332_GSE92332_SalmonellaInfect", "Haber_GSE92332_GSE92332_SalmHelm", "Haber_GSE92332_GSE92332_atlas", "Haber_GSE92332_FAE"]

    annot_lst = []
    for val  in adata[adata.obs.condition.isin(haber_samples),:].obs_names:
        annot_lst.append(val.split("_")[-1].split("-")[0])

    adata.obs.loc[adata.obs.condition.isin(haber_samples), "cell_type"] = annot_lst"""