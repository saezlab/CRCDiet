import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import argparse
import os


integrated_sample_path = "../data/out_data/integrated.h5ad"

def check_marker_vs_sample_gene_overlap():
    adata = sc.read_h5ad(integrated_sample_path)
    print(adata.var.index)

    df_markers = pd.read_csv("../data/Stroma_oldclusters_DEgenes_Venice_top100xlsx.csv", sep=";")
    marker_genes = []

    # adata_filtered_mouse.var.index = pd.Index(gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.index.values)
    # print([gen.split("mm10___")[1] for gen in adata_filtered_mouse.var.gene_ids.values])
    sample_gene_ids_set = set()
    for item in adata.var.index.values:
        sample_gene_ids_set.add(item.upper())
    print("Number of gene ids in the filtered data:", len(sample_gene_ids_set))

    all_marker_gene_set = set()
    for ind, row in df_markers.iterrows():

        for cell_type in df_markers.columns:
            marker_genes.append([row[cell_type], cell_type])
            all_marker_gene_set.add(row[cell_type])

    print("Number of gene ids in the marker gene set:", len(all_marker_gene_set))

    print("Overlap:", len(sample_gene_ids_set & all_marker_gene_set))

    marker_genes = pd.DataFrame(marker_genes, columns=['gene', 'cell_type']) 

check_marker_vs_sample_gene_overlap()