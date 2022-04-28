import utils
import qc_preprocess
import decoupler as dc
import pandas as pd

meta = utils.get_meta_data()

qc_preprocess.create_filtered_adata_files()



# qc_preprocess.analyze_pdx_samples("sample01")
"""marker_genes = dc.get_resource('PanglaoDB')
marker_db="PanglaoDB"
# Filter by canonical_marker and human
marker_genes = marker_genes[(marker_genes['mouse']=='True')&(marker_genes['canonical_marker']=='True')]

# Remove duplicated entries
marker_genes = marker_genes[~marker_genes.duplicated(['cell_type', 'genesymbol'])]
marker_genes.to_csv("../data/out_data/mouse_marker_genes_paglaodb.csv", index=False, index_label=False)"""