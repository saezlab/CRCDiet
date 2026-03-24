import scanpy as sc
import decoupler as dc

# Only needed for processing
import numpy as np
import pandas as pd
from matplotlib.pyplot import rc_context
# Needed for some plotting
import matplotlib.pyplot as plt

# Plotting options, change to your liking
sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))

# decoupler version was 1.2.0 and changed 18/04/2024 to 1.6.0
ref = "CD"
contrast = "LFD"

adata = sc.read_h5ad("../data/out_data/sc_integrated_cluster_scannot.h5ad")
adata.X = adata.layers["log1p_transformed"]
diet_lst = [col.split("-")[0] for col in adata.obs["condition"].values]
adata.obs["diet"] = diet_lst



adata.uns['log1p']["base"] = None
for cat in adata.obs["major_cell_types"].cat.categories:
    adata_tmp = adata[adata.obs["major_cell_types"]==cat,:].copy()
    # ref = "LFD"
    # contrast = "HFD"
    sc.tl.rank_genes_groups(adata_tmp, "diet", groups=[contrast], reference=ref, method="wilcoxon", key_added=f"wilcoxon_{cat}_{ref}_{contrast}")
    sc.pl.rank_genes_groups(adata_tmp, groups=[contrast], n_genes=20,use_raw=False, key=f"wilcoxon_{cat}_{ref}_{contrast}", save=f"{cat}_{contrast}_vs_{ref}_DEGs.pdf")
    with rc_context({"figure.figsize": (9, 10)}):
        sc.pl.rank_genes_groups_violin(adata_tmp, n_genes=20, jitter=False, use_raw=False, key=f"wilcoxon_{cat}_{ref}_{contrast}", save=f"{cat}_{contrast}_vs_{ref}_DEGs_violin.pdf")