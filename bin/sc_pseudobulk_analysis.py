import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import warnings
import utils
import gseapy
# Only needed for processing
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

############################### BOOOORIING STUFF BELOW ############################### 
# Warning settings
warnings.simplefilter(action='ignore')
sc.settings.verbosity = 0
# Set figure params
sc.set_figure_params(scanpy=True, facecolor="white", dpi=80, dpi_save=300)
# Read command line and set args
parser = argparse.ArgumentParser(prog='nnmfbc', description='Pseudobulk analysis')
parser.add_argument('-i', '--input_path', help='Input path to integrated object', required=True)
parser.add_argument('-o', '--output_dir', help='Output directory where to store the object', required=True)
parser.add_argument('-an', '--analysis_name', help='Analysis name', required=True)
parser.add_argument('-ss', '--subset', help="Cell type or  subgroup (e.g. cluster) of interest", required=True)
parser.add_argument('-con', '--condition', help='Condition sample', required=True)
parser.add_argument('-cont', '--contrast', help="Samples to be compared, E.g., [condition, B, A] will measure the "+ 
                    "LFC of condition B compared to condition A. If None, the last variable from the design matrix "+
                    "is chosen as the variable of interest, and the reference level is picked alphabetically", required=True)
parser.add_argument('-mct', '--major_cell_type', help='Major cell type (immune, epithelial, stroma)', required=False)

args = vars(parser.parse_args())
input_path = args['input_path']
output_path = args['output_dir']
analysis_name = args['analysis_name'] # subset
subset = args['subset'] # subset
condition = args['condition'] 
contrast = args['contrast']
major_cell_type = args['major_cell_type']
# Get necesary paths and create folders if necessary
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths(analysis_name, plot=True)
############################### BOOOORIING STUFF ABOVE ###############################

np.seterr(all="ignore")
sample_type = "sc"   
# save the cells coming from major cell type as a seperate object
adata_integ_clust= sc.read_h5ad(input_path)
immune_cell_types = [ "T cells", "Plasma cells", "Mast cells", "Neutrophils", "ILC2",  "Dendritic cells", "Myeloid cells", "B cells"]
epithelial_cell_types = ["Tuft cells", "Goblet cells", "Prolif.", "Enteroendocrine", "Keratynocytes", "Prolif. + Mature enterocytes"]
stroma_cell_types = ["Stroma", "Myofibroblasts", "Endothelial cells"]
lst_contrast = contrast.split(",")
# immune_samples = lst_contrast# , "LFD-AOM-DSS-Immune"]


adata_integ_clust = adata_integ_clust[adata_integ_clust.obs["condition"].isin(lst_contrast)]

if major_cell_type=="immune":
    adata_integ_clust = adata_integ_clust[adata_integ_clust.obs['cell_type_0.20'].isin(immune_cell_types)]
elif major_cell_type=="stroma":
    adata_integ_clust = adata_integ_clust[adata_integ_clust.obs['cell_type_0.20'].isin(stroma_cell_types)]
elif major_cell_type=="epithelial":
    adata_integ_clust = adata_integ_clust[adata_integ_clust.obs['cell_type_0.20'].isin(stroma_cell_types)]
else:
    raise Exception("Sorry, this is not among sample type!")


print("###########")
print(lst_contrast)
print("###########")

str_comparison = f"condition_{lst_contrast[0]}_vs_{lst_contrast[1]}"
# adata_integ_clust = adata_integ_clust[adata_integ_clust.obs["cell_type_0.20"]=="B cells"]


# Get pseudo-bulk profile
pdata = dc.get_pseudobulk(adata_integ_clust,
                          sample_col=condition,
                          groups_col='cell_type_0.20',
                          layer='counts',
                          mode='sum',
                          min_cells=0,
                          min_counts=0
                         )
# print(pdata.obs["cell_type_0.20"])


# dc.plot_psbulk_samples(pdata, groupby=[condition, 'cell_type_0.20'], save=f"{PLOT_PATH}/pseudobulk_{str_comparison}.pdf", figsize=(24, 6))

# print(pdata)
# Get filtered pseudo-bulk profile
pdata = dc.get_pseudobulk(adata_integ_clust,
                          sample_col=condition,
                          groups_col='cell_type_0.20',
                          layer='counts',
                          mode='sum',
                          min_cells=10,
                          min_counts=1000
                         )

dc.plot_psbulk_samples(pdata, groupby=[condition, 'cell_type_0.20'], save=f"{PLOT_PATH}/pseudobulk_after_filtering_{str_comparison}.pdf", figsize=(24, 6))


# Select T cell profiles
xcells = pdata[pdata.obs['cell_type_0.20'] == subset].copy()

# dc.plot_filter_by_expr(tcells, group='condition', min_count=10, save="testpsd2.pdf", min_total_count=15)

# Obtain genes that pass the thresholds
genes = dc.filter_by_expr(xcells, group=condition, min_count=10, min_total_count=15)
# Filter by these genes
xcells = xcells[:, genes].copy()


dds = DeseqDataSet(
    adata=xcells,
    design_factors="condition",  # compare samples based on the "condition"
    # column ("B" vs "A")
    refit_cooks=True,
    n_cpus=8,
)


dds.fit_size_factors()
print(dds.obsm["size_factors"])

dds.fit_genewise_dispersions()
print(dds.varm["genewise_dispersions"])

dds.fit_dispersion_trend()
print(dds.uns["trend_coeffs"])
print(dds.varm["fitted_dispersions"])

dds.fit_dispersion_prior()
print(
    f"logres_prior={dds.uns['_squared_logres']}, sigma_prior={dds.uns['prior_disp_var']}"
)

dds.fit_MAP_dispersions()
print(dds.varm["MAP_dispersions"])
print(dds.varm["dispersions"])

dds.fit_LFC()
print(dds.varm["LFC"])

dds.calculate_cooks()
if dds.refit_cooks:
    # Replace outlier counts
    dds.refit()

stat_res = DeseqStats(dds, contrast= [condition]+ lst_contrast, alpha=0.05, cooks_filter=True, independent_filter=True)


stat_res.run_wald_test()
stat_res.p_values

if stat_res.independent_filter:
    stat_res._independent_filtering()
else:
    stat_res._p_value_adjustment()

stat_res.padj

stat_res.summary()

results_df = stat_res.results_df
# results_df.to_csv(f"{DATA_PATH}/analysis/{str_comparison}_before_shrink_deg.csv")
# dc.plot_volcano_df(results_df, x='log2FoldChange', y='padj', top=20, figsize=(14,10), dpi=300, save=f"{PLOT_PATH}/volcano_before_{str_comparison}.pdf")

stat_res.lfc_shrink(coeff=str_comparison)

results_df = stat_res.results_df


results_df.to_csv(f"{DATA_PATH}/analysis/{str_comparison}_deg.csv")

dc.plot_volcano_df(results_df, x='log2FoldChange', y='padj', top=20, figsize=(14,10), dpi=300, save=f"{PLOT_PATH}/volcano_after_{str_comparison}.pdf")
lFCs_thr, sign_thr = 0.5, 0.05
two_cond_upregul = (results_df['log2FoldChange'] >= lFCs_thr) & (results_df['padj'] <= sign_thr)

# sort by padj
results_df = results_df[two_cond_upregul].sort_values('padj', ascending=True)
results_df.to_csv(f"{PLOT_PATH}/{str_comparison}_upregulated_sorted_by_padj.csv")
lst_two_cond_upreg_genes= list(results_df[two_cond_upregul].sort_values('padj', ascending=True).index)
lst_two_cond_upreg_genes = [gene_id.upper() for gene_id in lst_two_cond_upreg_genes]

print(lst_two_cond_upreg_genes)

enr_res = gseapy.enrichr(gene_list=lst_two_cond_upreg_genes,
                     organism='Mouse',
                     gene_sets='GO_Biological_Process_2021',
                     cutoff = 0.05)

print(enr_res.res2d)
gseapy.barplot(enr_res.res2d, title='GO_Biological_Process_2021', ofname=f"{PLOT_PATH}/barplot_upregulated_gea_{str_comparison}.pdf")

"""
python sc_pseudobulk_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an sc_bcells_pseudobulk_deg_analysis -ss "B cells" -con condition -cont "LFD-AOM-DSS-Immune,CD-AOM-DSS-Immune" -mct "immune"
python sc_pseudobulk_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an sc_bcells_pseudobulk_deg_analysis -ss "B cells" -con condition -cont "HFD-AOM-DSS-Immune,CD-AOM-DSS-Immune" -mct "immune"
python sc_pseudobulk_analysis.py -i ../data/out_data/sc_integrated_cluster_scannot.h5ad -o ../data/out_data -an sc_bcells_pseudobulk_deg_analysis -ss "B cells" -con condition -cont "LFD-AOM-DSS-Immune,HFD-AOM-DSS-Immune" -mct "immune"
"""

