from reprlib import aRepr
from pathlib import Path
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd

import argparse
import os


# Read command line and set args
parser = argparse.ArgumentParser(prog='qc', description='Run intergration by Harmony')
parser.add_argument('-i', '--input_path', help='Input path to merged object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-st', '--sample_type', help='Human, mouse or tumor', required=True)
args = vars(parser.parse_args())

input_path = args['input_path']
output_path = args['output_dir']
sample_type = args["sample_type"]
###############################

plot_path="../plots/integrate"
Path(plot_path).mkdir(parents=True, exist_ok=True)
sc.settings.figdir = plot_path


# Read merged object
adata = sc.read_h5ad(input_path)
# Run harmony
sce.pp.harmony_integrate(adata, 'batch', adjusted_basis='X_pca', max_iter_harmony=30)

# Run umap with updated connectivity
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(
    adata, color=["condition"], palette=sc.pl.palettes.default_20,  show=False, save=f"{sample_type}_harmony"
)

# Write to file
adata.write(os.path.join(output_path, f'{sample_type}_integrated.h5ad'))

#  python integrate.py -i ../data/out_data/mouse_merged.h5ad -o ../data/out_data -st mouse
#  python integrate.py -i ../data/out_data/human_merged.h5ad -o ../data/out_data -st human