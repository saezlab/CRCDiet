import matplotlib as mpl
import pandas as pd
import scanpy as sc
import os

dot_path = "/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/csv"
vis_path = "/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data"
df_dot = pd.read_csv(os.path.join(dot_path, "DOT_HFD-AOM_sparsity_1_ratios_0.csv"), index_col=0)

df_dot = pd.read_csv(os.path.join(dot_path, "DOT_HFD-AOM_sparsity_0.5_ratios_0.5.csv"), index_col=0)
df_dot = pd.read_csv(os.path.join(dot_path, "DOT_HFD-AOM_sparsity_1_ratios_0.5.csv"), index_col=0)
df_dot = pd.read_csv(os.path.join(dot_path, "DOT_HFD-AOM_sparsity_0.5_ratios_0.csv"), index_col=0)
adata = sc.read_h5ad(os.path.join(vis_path, "ext_L24854_HFD-AOM-DSS-colon-d81-visium_filtered.h5ad"))


for col in df_dot.columns:
    adata.obs[col] = df_dot[col]

with mpl.rc_context({'axes.facecolor':  'black', 'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(adata, cmap='magma', color=df_dot.columns, ncols=4, size=1.3, img_key='hires', vmin=0, vmax='p99.2')

