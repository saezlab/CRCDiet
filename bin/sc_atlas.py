import os
import utils
import scanpy as sc
import pandas as pd
from utils import OUT_DATA_PATH, PLOT_PATH, DATA_PATH, read_raw_sc_sample


S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_n_return_paths("sc_atlas_save")

lst_adata = []

toplam = 0





lst_Drokhlyansky_SCP1038 = ["msi.neur", "msi.glia", "mli.neur", "mli.glia", "msi", "mli"]
lst_Drokhlyansky_SCP1038 = ["mli"]

for samp_id in lst_Drokhlyansky_SCP1038:
    anndata = sc.read_mtx(os.path.join("../data", "atlas", "Drokhlyansky_SCP1038", f"gene_sorted-{samp_id}.matrix.mtx")).transpose()
    features = list(pd.read_csv(os.path.join("../data", "atlas", "Drokhlyansky_SCP1038", f"{samp_id}.genes.tsv"), sep="\t", header=None)[0])
    barcodes = list(pd.read_csv(os.path.join("../data", "atlas", "Drokhlyansky_SCP1038", f"{samp_id}.barcodes.tsv"), sep="\t", header=None)[0])
    anndata.var_names = features
    anndata.obs_names = barcodes
    anndata.var_names_make_unique()
    anndata.write(os.path.join(OUT_DATA_PATH, f"Drokhlyansky_SCP1038_{samp_id}.h5ad"))
    lst_adata.append(anndata)
    # anndata.var_names = [f_name.upper() for f_name in list(anndata.var_names)]
    toplam += anndata.shape[0]
    print(toplam)
    print(anndata)

"""
lst_Ayyaz_GSE123516 = ["GSM3308717_C04", "GSM3308718_C05", "GSM3308719_C06", "GSM3308720_C07"]

for samp_id in lst_Ayyaz_GSE123516:
    anndata = sc.read_mtx(os.path.join("../data", "atlas", "Ayyaz_GSE123516", f"{samp_id}.mtx.gz")).transpose()
    features = list(pd.read_csv(os.path.join("../data", "atlas", "Ayyaz_GSE123516", f"{samp_id}_genes.tsv.gz"), sep="\t", header=None)[1])
    barcodes = list(pd.read_csv(os.path.join("../data", "atlas", "Ayyaz_GSE123516", f"{samp_id}_barcodes.tsv.gz"), sep="\t", header=None)[0])
    anndata.var_names = features
    anndata.obs_names = barcodes
    anndata.var_names_make_unique()
    anndata.write(os.path.join(OUT_DATA_PATH, f"Ayyaz_GSE123516_{samp_id}.h5ad"))
    # anndata.var_names = [f_name.upper() for f_name in list(anndata.var_names)]

    lst_adata.append(anndata)
    toplam += anndata.shape[0]
    print(toplam) 
# ['B cells-1' 'Neutrophils' 'IgA plasma cells-1']



# a  = pd.read_csv("GSE92332_atlas_UMIcounts.txt")
# a  = pd.read_csv("GSE92332_SalmHelm_UMIcounts.txt", sep = "\t")
# Haber_GSE92332
lst_Haber_GSE92332 = ["FAE_UMIcounts.txt", "GSE92332_atlas_UMIcounts.txt", "GSE92332_SalmHelm_UMIcounts.txt", "GSE92332_SalmonellaInfect_UMIcounts.txt"]

for samp_id in lst_Haber_GSE92332:
    anndata = sc.AnnData(pd.read_csv(os.path.join("../data", "atlas", "Haber_GSE92332", samp_id), sep = "\t").transpose())
    # anndata.var_names = [f_name.upper() for f_name in list(anndata.var_names)]
    anndata.var_names_make_unique()
    anndata.write(os.path.join(OUT_DATA_PATH, f"Haber_GSE92332_{samp_id}.h5ad"))
    lst_adata.append(anndata)
    toplam += anndata.shape[0]
    print(anndata)





#lst_Kim_GSE116514 = ["GSM3242184_Sample_DS39_gene_exon_tagged.dge.txt.gz", "GSM3242185_Sample_DS19_gene_exon_tagged.dge.txt.gz"]

#for samp_id in lst_Kim_GSE116514:
#    anndata = sc.AnnData(pd.read_csv(os.path.join("../data", "atlas", "Kim_GSE116514", samp_id), sep = "\t").transpose())
#    toplam += anndata.shape[0]
#    print(toplam)
#    print(anndata)

lst_Morarach_GSE149524 = ["GSM4504451_P21_2", "GSM4504448_E15", "GSM4504449_E18", "GSM4504450_P21_1"]

for samp_id in lst_Morarach_GSE149524:
    anndata = sc.read_mtx(os.path.join("../data", "atlas", "Morarach_GSE149524", f"{samp_id}_matrix.mtx.gz")).transpose()
    features = list(pd.read_csv(os.path.join("../data", "atlas", "Morarach_GSE149524", f"{samp_id}_features.tsv.gz"), sep="\t", header=None)[1])
    barcodes = list(pd.read_csv(os.path.join("../data", "atlas", "Morarach_GSE149524", f"{samp_id}_barcodes.tsv.gz"), sep="\t", header=None)[0])
    anndata.var_names = features
    anndata.obs_names = barcodes
    anndata.var_names_make_unique()
    anndata.write(os.path.join(OUT_DATA_PATH, f"Morarach_GSE149524_{samp_id}.h5ad"))
    lst_adata.append(anndata)
    # anndata.var_names = [f_name.upper() for f_name in list(anndata.var_names)]
    toplam += anndata.shape[0]
    print(toplam)
    print(anndata)




lst_Niec_GSE190037 = ["GSM5712427_SISC_Lymphatics",  "GSM5712428_Colon_Lymphatics"]

for samp_id in lst_Niec_GSE190037:
    anndata = sc.read_mtx(os.path.join("../data", "atlas", "Niec_GSE190037", f"{samp_id}_matrix.mtx.gz")).transpose()
    features = list(pd.read_csv(os.path.join("../data", "atlas", "Niec_GSE190037", f"{samp_id}_features.tsv.gz"), sep="\t", header=None)[1])
    barcodes = list(pd.read_csv(os.path.join("../data", "atlas", "Niec_GSE190037", f"{samp_id}_barcodes.tsv.gz"), sep="\t", header=None)[0])
    anndata.var_names = features
    anndata.obs_names = barcodes
    anndata.var_names_make_unique()
    anndata.write(os.path.join(OUT_DATA_PATH, f"Niec_GSE190037_{samp_id}.h5ad"))
    lst_adata.append(anndata)
    # anndata.var_names = [f_name.upper() for f_name in list(anndata.var_names)]
    toplam += anndata.shape[0]
    print(toplam)
    print(anndata)



lst_Parigi_GSE163638 = ["GSM4983265_IEC-Stroma-control", "GSM4983266_IEC-Stroma-depleted"]

for samp_id in lst_Parigi_GSE163638:
    anndata = sc.read_mtx(os.path.join("../data", "atlas", "Parigi_GSE163638", f"{samp_id}-matrix.mtx.gz")).transpose()
    features = list(pd.read_csv(os.path.join("../data", "atlas", "Parigi_GSE163638", f"{samp_id}-features.tsv.gz"), sep="\t", header=None)[1])
    barcodes = list(pd.read_csv(os.path.join("../data", "atlas", "Parigi_GSE163638", f"{samp_id}-barcodes.tsv.gz"), sep="\t", header=None)[0])
    anndata.var_names = features
    anndata.obs_names = barcodes
    anndata.var_names_make_unique()
    anndata.write(os.path.join(OUT_DATA_PATH, f"Parigi_GSE163638_{samp_id}.h5ad"))
    lst_adata.append(anndata)
    # anndata.var_names = [f_name.upper() for f_name in list(anndata.var_names)]
    toplam += anndata.shape[0]
    print(toplam)
    print(anndata)




anndata = sc.read_mtx(os.path.join("../data", "atlas", "Xu_GSE124880", "GSE124880_PP_LP_mm10_count_matrix.mtx")).transpose()
features = list(pd.read_csv(os.path.join("../data", "atlas", "Xu_GSE124880", "GSE124880_PP_LP_mm10_count_gene.tsv"), sep="\t", header=None)[0])
barcodes = list(pd.read_csv(os.path.join("../data", "atlas", "Xu_GSE124880", "GSE124880_PP_LP_mm10_count_barcode.tsv"), sep="\t", header=None)[0])
anndata.var_names = features
anndata.obs_names = barcodes
anndata.var_names_make_unique()
anndata.write(os.path.join(OUT_DATA_PATH, f"Xu_GSE124880_GSE124880_PP_LP_mm10.h5ad"))
lst_adata.append(anndata)
# anndata.var_names = [f_name.upper() for f_name in list(anndata.var_names)]
toplam += anndata.shape[0]
print(toplam)
print(anndata)


"""