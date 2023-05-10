from cProfile import label
from turtle import color
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib import cm
import scanpy as sc
# import utils
import numpy as np
import matplotlib
import os
import utils
import matplotlib.image as mpimg

'''Plotting functions'''

cond_colors = {
    'Healthy' : '#440154FF',
    'HF-A' : '#21908CFF',
    'HF-CKD' : '#FDE725FF',
}

ctype_colors = {
    'adipocytes' : '#D51F26',
    'cardiomyocyte' : '#272E6A',
    'endothelial' : '#208A42',
    'fibroblast' : '#89288F',
    'lymphatic_endo' : '#F47D2B',
    'macrophages' : '#FEE500',
    'mast_cells' : '#8A9FD1',
    'neuronal' : '#C06CAB',
    'pericyte' : '#D8A767',
    'T-cells' : '#90D5E4',
    'vSMCs' : '#89C75F'
}

def plot_mt_vs_counts(data, ax, mt_thr=20, fontsize=11):    
    # Plot scatter
    ax.scatter(x=data.obs.total_counts, y=data.obs.pct_counts_mt, s=1, c='gray')
    ax.axhline(y=mt_thr, linestyle='--', color="black")
    ax.set_xlabel("Total counts", fontsize=fontsize)
    ax.set_ylabel("Fraction MT counts", fontsize=fontsize)


def plot_rp_vs_counts(data, ax, rp_thr=None, fontsize=11):    
    # Plot scatter
    ax.scatter(x=data.obs.total_counts, y=data.obs.pct_counts_rp, s=1, c='gray')
    if rp_thr:
        ax.axhline(y=rp_thr, linestyle='--', color="black")
    ax.set_xlabel("Total counts", fontsize=fontsize)
    ax.set_ylabel("Fraction RP counts", fontsize=fontsize)


def plot_ngenes_vs_counts(data, ax, gene_thr=6000, fontsize=11):
    # Plot scatter
    ax.scatter(x=data.obs.total_counts, y=data.obs.n_genes_by_counts, s=1, c='gray')
    ax.axhline(y=gene_thr, linestyle='--', color="black")
    ax.set_xlabel("Total counts", fontsize=fontsize)
    ax.set_ylabel("Number of genes expr", fontsize=fontsize)


def plot_doublet_scores(data, ax, doublet_thr=0.2, fontsize=11):
    # Plot histogram
    ax.hist(data.obs.doublet_score, bins=100, color='gray')
    ax.axvline(x=doublet_thr, linestyle='--', color="black")
    ax.set_xlabel('Droplet score distribution', fontsize=fontsize)
    

def plot_diss_scores(data, ax, diss_thr=0.5, fontsize=11):
    # Plot histogram
    ax.hist(data.diss_score, bins=100, color='gray')
    ax.axvline(x=diss_thr, linestyle='--', color="black")
    ax.set_xlabel('Dissociation score distribution', fontsize=fontsize)


def plot_ncell_diff(data, ax, labels, n_rem, fontsize=11):
    # Plot Barplot
    for label, n in zip(labels, n_rem):
        ax.bar(label, n)
    ax.set_title('Cells removed per filter', fontsize=fontsize)
    ax.tick_params(axis='x', rotation=45)
    
def plot_cell_type_proportion(cond_list, cond_name ="Immune", adata=None, obs_col = "cell_type_0.20", sample_type="sc"):
    input_path = "/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/sc_integrated_cluster_scannot.h5ad"
    sample_type="sc"
    adata = sc.read_h5ad(input_path)
    meta = utils.get_meta_data(sample_type)
    condition = list(np.unique(meta['condition']))
    c_type_list = list(adata.obs[obs_col].cat.categories)

    # CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN,HFD-AOM-DSS-Epi_plus_DN
    # CD-AOM-DSS-Immune,LFD-AOM-DSS-Immune,HFD-AOM-DSS-Immune
    cond_list = cond_list.split(",")
    # cond_list = ["CD-AOM-DSS-Immune", "LFD-AOM-DSS-Immune", "HFD-AOM-DSS-Immune"]
    # cond_list = ["CD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN"]
    cond_prop = dict()
    cond_arr = []
    for cond in cond_list:
        cond_arr.append([])
        # print(cond, cond_arr)
        cond_prop[cond] = []
        adata_tmp = adata[adata.obs["condition"]==cond,:]
        # print(adata_tmp.shape)
        sum = 0
        for c_type in c_type_list:
            print("c_type", c_type, adata_tmp[adata_tmp.obs[obs_col]==c_type].shape)
            cond_arr[-1].append(100*(adata_tmp[adata_tmp.obs[obs_col]==c_type].shape[0]/adata_tmp.shape[0]))
            #cond_prop[cond][c_type] = adata_tmp.obs["cell_type_0.20"].str.count(c_type).sum()/adata_tmp.shape[0]
            # cond_prop[cond][c_type] = adata_tmp.obs["cell_type_0.20"].str.count(c_type).sum()/adata_tmp.shape[0]
            #sum += adata_tmp.obs["cell_type_0.20"].str.count(c_type).sum()

    data = np.array(cond_arr).T

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # For loop for creating stacked bar chart
    cmap = matplotlib.cm.get_cmap('tab20')

    X = np.arange(data.shape[1])
    for i in range(data.shape[0]):
        ax1.bar(X, data[i],bottom = np.sum(data[:i], 
                    axis =0), width= 0.85, color = cmap.colors[i], label=c_type_list[i]  )

    ax1.set_xticks([0,1,2])
    ax1.set_xticklabels(cond_list) # , rotation=45)
    ax1.set_xlabel("Condition", fontweight='bold')
    ax1.set_ylabel("Proportion (%)", fontweight='bold')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.subplots_adjust(bottom=0.45)
    fig.tight_layout()
    plt.savefig(f"../plots/sc_cell_type_annot/major_cell_type_prop_{cond_name}_barplot.pdf")

    fig, axs = plt.subplots(5, 4, figsize=[15, 15])
    
    for ind, c_t in enumerate(c_type_list):
        fig_row, fig_col = int(ind/4), ind%4
        axs[fig_row][fig_col].plot(cond_list, data[ind,:], marker='o', mec = 'grey', mfc = 'grey', markersize=12, linewidth=5, color=cmap.colors[ind], label=c_t)
        axs[fig_row][fig_col].set_title(c_t)
        axs[fig_row][fig_col].set_xticklabels(["CD", "LFD", "HFD"])
        axs[fig_row][fig_col].tick_params(axis='y', which='major', labelsize=8)

    
    plt.subplots_adjust(bottom=0.45)
    fig.supxlabel('Condition', fontweight='bold')
    fig.supylabel('Proportion (%)', fontweight='bold')
    fig.tight_layout()
    
    fig.delaxes(axs[4][3])
    fig.delaxes(axs[4][2])
    fig.delaxes(axs[4][1])
    plt.savefig(f"../plots/sc_cell_type_annot/majorcell_type_prop_{cond_name}_line.pdf")



# plot_cell_type_proportion("CD-AOM-DSS-Epi_plus_DN,LFD-AOM-DSS-Epi_plus_DN,HFD-AOM-DSS-Epi_plus_DN", cond_name ="epithelial", adata=None, obs_col = "major_cell_types", sample_type="sc")
# plot_cell_type_proportion("CD-AOM-DSS-Immune,LFD-AOM-DSS-Immune,HFD-AOM-DSS-Immune", cond_name ="immune", adata=None, obs_col = "major_cell_types", sample_type="sc")


def plot_ngene_diff(adata, ax, fontsize=11):
    ax.set_title('Num genes filtered', fontsize=fontsize)
    ax.bar(x="Before", height=adata.uns['hvg']['ngene'])
    ax.bar(x="After", height=adata.shape[1])
    
    
def plot_hvg_nbatches(data, ax, fontsize=11):
    for nbatches in np.flip(np.unique(data.highly_variable_nbatches)):
        num = data[data.highly_variable_nbatches == nbatches].shape[0]
        ax.bar(str(nbatches), num, color="gray")
    ax.set_title('Num shared HVG by num samples',fontsize=fontsize)
    
    
def plot_sorted_rank(data, col, ax, fontsize=11):
    xdim = np.arange(len(data))
    ysort = np.flip(np.sort(data[col]))
    ax.set_title('Ranked {0}'.format(col), fontsize=fontsize)
    ax.plot(xdim, ysort, c='grey')

    
def stacked_barplot(data, feature_name, ax, cmap=cm.tab20):
    # cell type names
    type_names = data.var.index
    levels = pd.unique(data.obs[feature_name])
    n_levels = len(levels)
    feature_totals = np.zeros([n_levels, data.X.shape[1]])

    for level in range(n_levels):
        l_indices = np.where(data.obs[feature_name] == levels[level])
        feature_totals[level] = np.sum(data.X[l_indices], axis=0)
        
    y = feature_totals
    title=feature_name
    level_names=levels
    
    n_bars, n_types = y.shape
    r = np.array(range(n_bars))
    sample_sums = np.sum(y, axis=1)

    barwidth = 0.85
    cum_bars = np.zeros(n_bars)

    for n in range(n_types):
        bars = [i / j * 100 for i, j in zip([y[k][n] for k in range(n_bars)], sample_sums)]
        ax.bar(r, bars, bottom=cum_bars, color=cmap(n % cmap.N), width=barwidth, label=type_names[n], linewidth=0)
        cum_bars += bars
    ax.set_title(title)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    ax.set_xticks(r)
    ax.set_xticklabels(level_names, rotation=45)
    ax.set_ylabel("Proportion")

def volcano(name, lfc, pvals, ax, max_num=None,
                 lfc_thr=0.5, p_thr=0.05, s=10, fontsize=12):
    '''Volcano plot from a list of lfc and untrasnformed pvalues'''

    if max_num is None:
        max_num=np.max(np.abs(lfc))
    
    # Transform pvals
    pvals = -np.log10(pvals)
    
    # Mask significant genes
    msk = (pvals > -np.log10(p_thr)) & (np.abs(lfc) > lfc_thr)
    
    # Plot scatter
    ax.set_title(name)
    ax.scatter(lfc[~msk], pvals[~msk], c='gray', s=s)
    ax.scatter(lfc[msk], pvals[msk], c='red', s=s)
    ax.set_xlim(-max_num, max_num)
    ax.set_xlabel('LogFC', fontsize=fontsize)
    ax.set_ylabel('-log10(pvalue)', fontsize=fontsize)
    ax.set_box_aspect(1)
    
def dotplot(title, x, y, c, s, size_title, color_title, cmap='coolwarm', edgecolor=None, num=30, fontsize=9, figsize=(12,6)):
    # Define figure
    fig, ax = plt.subplots(1,1, dpi=150, figsize=figsize)
    ax.set_title(title, fontsize=fontsize+5)
    
    # Add grid and set it to background
    ax.grid(True)
    ax.set_axisbelow(True)
    
    # Dot plot
    max_num = np.max(np.abs(c))
    scatter = ax.scatter(
        x=x,
        y=y,
        c=c,
        s=s * num,
        cmap=cmap,
        vmax=max_num,
        vmin=-max_num,
        edgecolor=edgecolor
    )
    
    # Format dot plot ticks
    ax.tick_params(axis='x', rotation=90, labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.margins(y=0.05)

    # Plot pvalue dot sizes legend
    handles, labels = scatter.legend_elements("sizes", num=4)
    labels = ['$\\mathdefault{'+'{0}'.format(int(int(label.split('{')[1].split('}')[0])/num))+'}$' for label in labels]

    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1,1), frameon=False, title=size_title)
    
    # Add color bar
    cax = fig.add_axes([0.945, 0.25, 0.025, 0.35])
    cbar = fig.colorbar(scatter, cax=cax, orientation='vertical')
    cbar.ax.set_title(color_title)
    
    # Format figure
    fig.tight_layout()
    fig.set_facecolor('white')
    
    return fig

def plot_ora(name, df, ax, top=10, fontsize=11):
    df = df.sort_values('adj_pvalues', ascending=True).head(top)
    names = np.flip(df['descr'].tolist())
    pvals = np.flip(-np.log10(df['adj_pvalues']))
    ax.barh(names, pvals, color='gray')
    ax.axvline(x=-np.log10(0.05), c='black', ls='--')
    ax.set_xlabel('-log10(adj_pval)', fontsize=fontsize)
    ax.set_title(name, fontsize=fontsize)
    
def corr(name, x, y, ax, fontsize=11):
    from scipy import stats
    corr, _ = stats.spearmanr(x, y)
    ax.scatter(x, y, c='gray', s=5)
    ax.set_title('{0} | corr: {1}'.format(name, '{:.2f}'.format(corr)), fontsize=fontsize)
    
def violins(arr, meta, ax):
    data = arr.copy()
    data[data==0] = np.nan
    sns.violinplot(x=data.melt().variable, 
                   y=np.log2(data.melt().value), 
                   hue=meta.loc[data.melt().variable].condition.values)


# def plot_qc_after_filtering(sample_type):
def show_plot(plt_path, text=None):
    img = mpimg.imread(plt_path)
    plt.imshow(img)

def show_qc_filtering_plot(plot_fold_path, sample_type):
    meta = utils.get_meta_data(sample_type)

    for _, row in meta.iterrows():
        sample_id = row["sample_id"]

        print("QC metrics...")
        show_plot(os.path.join(plot_fold_path, f"{sample_type}_qc_preprocess", f"basic_stats_before_filtering_{sample_id}.pdf"))
        print("Plotting highest expressed genes after QC and filtering...")
        show_plot(os.path.join(plot_fold_path, f"{sample_type}_qc_preprocess", f"highest_expr_genesbasic_stats_after_filtering_{sample_id}.pdf"))

def plot_clusters(adata):
    plt.rcParams['figure.dpi']= 300
    plt.rcParams['figure.figsize']= (45, 30)
    for cat_n in adata.obs['leiden_0.40'].cat.categories:
        cat_n = int(cat_n)
        adata.obs['cluster_dummy'] = adata.obs['leiden_0.40'] == adata.obs['leiden_0.40'].cat.categories[cat_n]
        adata.obs["cluster_dummy"] = adata.obs['cluster_dummy'].astype(str).astype('category')
        sc.pl.umap(adata, color='cluster_dummy', size=10, title=f"Cluster {cat_n}", save=f"atlas_cluster_{cat_n}")

"""        adata.obs['cluster_dummy'] = adata.obs['leiden_0.40'] == adata.obs['leiden_0.40'].cat.categories[1]
        adata.obs["cluster_dummy"] = adata.obs['cluster_dummy'].astype(str).astype('category')
        sc.pl.umap(adata, color='cluster_dummy', size=10, title="Cluster 0", save=f"atlas_cluster_1")"""


# the functions for adding stat. sig. plots are coming from here
# https://stackoverflow.com/questions/11517986/indicating-the-statistically-significant-difference-in-bar-graph


# This plotting function comes from the answer in the following link
# https://stackoverflow.com/questions/11517986/indicating-the-statistically-significant-difference-in-bar-graph
def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    plt.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)


def test_significance():
    heights = [1.8, 2, 3]
    bars = np.arange(len(heights))

    plt.figure()
    plt.bar(bars, heights, align='center')
    plt.ylim(0, 5)
    barplot_annotate_brackets(0, 1, .1, bars, heights)
    barplot_annotate_brackets(1, 2, .001, bars, heights)
    barplot_annotate_brackets(0, 2, 'p < 0.0075', bars, heights, dh=.2)
    plt.savefig("test_sig.pdf")




# Custom function to draw the diff bars

def label_diff(i,j,text,X,Y, ax):
    x = (X[i]+X[j])/2
    print(Y[i], Y[j])

    y = max(Y[i], Y[j])
    
    print(y)
    dx = abs(X[i]-X[j])

    props = {'connectionstyle':'bar','arrowstyle':'-',\
                 'shrinkA':5,'shrinkB':5,'linewidth':2}

    ax.annotate(text, xy=(x,y), zorder=10, ha='center', fontsize=15)
    ax.annotate('', xy=(X[i],y), xytext=(X[j],y), arrowprops=props)
    



def plot_significance(first_label, second_label, third_label, first_prop, second_prop, third_prop, first_second_p_val, second_third_p_val, first_third_pval, c_type, fl_name):
    proportions   = (first_prop, second_prop, third_prop)
    # print(proportions, max(proportions))
    ind  = np.arange(3)    # the x locations for the groups
    width= 0.7
    labels = (first_label, second_label, third_label)

    # Pull the formatting out here
    bar_kwargs = {'width':width,'color':'y','linewidth':2,'zorder':5}
    err_kwargs = {'zorder':0,'fmt':None,'linewidth':2,'ecolor':'k'}  #for matplotlib >= v1.4 use 'fmt':'none' instead
    plt.rcParams['figure.dpi']= 300
    plt.rcParams['figure.figsize']= (15, 10)

    fig, ax = plt.subplots()
    ax.p1 = plt.bar(ind, proportions, **bar_kwargs)
    # ax.errs = plt.errorbar(ind, menMeans, yerr=menStd, **err_kwargs)
    # Call the function
    label_diff(0,1,f'{first_second_p_val}',ind,proportions, ax)
    label_diff(1,2,f'{second_third_p_val}',ind,proportions, ax)
    label_diff(0,2,f'{first_third_pval}',ind,proportions, ax)
    # ax.set_ylim(max(proportions)+ max(proportions)*0.20)

    ax.set_ylim([0.0, (max(proportions)+ max(proportions)*0.30)])
    ax.set_title(f'Cell Type Proportion: {c_type}')
    ax.set_xlabel('Condition', fontsize=20)
    ax.set_ylabel('Proportion (%)', fontsize=20)
    ax.spines[['right', 'top']].set_visible(False)



    
    plt.xticks(ind, labels, color='k')
    plt.savefig(f"../plots/sc_cell_type_prop/{fl_name}_{c_type}.pdf")
    plt.savefig(f"../plots/sc_cell_type_prop/{fl_name}_{c_type}.png")


def deg_venn_diagram(sign_thr=0.05, lFCs_thr=0.5, topn=50):
    from matplotlib_venn import venn3

    df_hfd_vs_cd = pd.read_csv("../data/analysis/condition_HFD-AOM-DSS-Immune_vs_CD-AOM-DSS-Immune_deg.csv")
    df_lfd_vs_cd = pd.read_csv("../data/analysis/condition_LFD-AOM-DSS-Immune_vs_CD-AOM-DSS-Immune_deg.csv")
    df_lfd_vs_hfd = pd.read_csv("../data/analysis/condition_LFD-AOM-DSS-Immune_vs_HFD-AOM-DSS-Immune_deg.csv")


    """df_hfd_vs_cd = pd.read_csv("../data/analysis/condition_HFD-AOM-DSS-Immune_vs_CD-AOM-DSS-Immune_before_shrink_deg.csv")
    df_lfd_vs_cd = pd.read_csv("../data/analysis/condition_LFD-AOM-DSS-Immune_vs_CD-AOM-DSS-Immune_before_shrink_deg.csv")
    df_lfd_vs_hfd = pd.read_csv("../data/analysis/condition_LFD-AOM-DSS-Immune_vs_HFD-AOM-DSS-Immune_before_shrink_deg.csv")"""

    df_hfd_vs_cd_mask_upreg = (df_hfd_vs_cd['log2FoldChange'] >= lFCs_thr) & (df_hfd_vs_cd['padj'] <= sign_thr)
    df_hfd_vs_cd_mask_downreg = (df_hfd_vs_cd['log2FoldChange'] <= -lFCs_thr) & (df_hfd_vs_cd['padj'] <= sign_thr)
    set_hfd_vs_cd_mask_upreg= set(df_hfd_vs_cd[df_hfd_vs_cd_mask_upreg].sort_values('padj', ascending=True).iloc[:topn]["Unnamed: 0"])
    set_hfd_vs_cd_mask_downreg= set(df_hfd_vs_cd[df_hfd_vs_cd_mask_downreg].sort_values('padj', ascending=True).iloc[:topn]["Unnamed: 0"])

    df_lfd_vs_cd_mask_upreg = (df_lfd_vs_cd['log2FoldChange'] >= lFCs_thr) & (df_lfd_vs_cd['padj'] <= sign_thr)
    df_lfd_vs_cd_mask_downreg = (df_lfd_vs_cd['log2FoldChange'] <= -lFCs_thr) & (df_lfd_vs_cd['padj'] <=sign_thr)
    set_lfd_vs_cd_mask_upreg= set(df_lfd_vs_cd[df_lfd_vs_cd_mask_upreg].sort_values('padj', ascending=True).iloc[:topn]["Unnamed: 0"])
    set_lfd_vs_cd_mask_downreg= set(df_lfd_vs_cd[df_lfd_vs_cd_mask_downreg].sort_values('padj', ascending=True).iloc[:topn]["Unnamed: 0"])

    df_lfd_vs_hfd_mask_upreg = (df_lfd_vs_hfd['log2FoldChange'] >= lFCs_thr) & (df_lfd_vs_hfd['padj'] <= sign_thr)
    df_lfd_vs_hfd_mask_downreg = (df_lfd_vs_hfd['log2FoldChange'] <= -lFCs_thr) & (df_lfd_vs_hfd['padj'] <= sign_thr)
    set_lfd_vs_hfd_mask_upreg= set(df_lfd_vs_hfd[df_lfd_vs_hfd_mask_upreg].sort_values('padj', ascending=True).iloc[:topn]["Unnamed: 0"])
    set_lfd_vs_hfd_mask_downreg= set(df_lfd_vs_hfd[df_lfd_vs_hfd_mask_downreg].sort_values('padj', ascending=True).iloc[:topn]["Unnamed: 0"])

    # venn3([set_hfd_vs_cd_mask_upreg.union(set_hfd_vs_cd_mask_downreg), set_lfd_vs_cd_mask_upreg.union(set_lfd_vs_cd_mask_downreg), set_lfd_vs_hfd_mask_upreg.union(set_lfd_vs_hfd_mask_downreg)], ('HFD vs. CD', 'LFD vs. CD', 'LFD vs. HFD'))
    venn3([set_hfd_vs_cd_mask_upreg, set_lfd_vs_cd_mask_upreg, set_lfd_vs_hfd_mask_upreg], ('HFD vs. CD', 'LFD vs. CD', 'LFD vs. HFD'))

 
    plt.show()
    venn3([set_hfd_vs_cd_mask_downreg, set_lfd_vs_cd_mask_downreg, set_lfd_vs_hfd_mask_downreg], ('HFD vs. CD', 'LFD vs. CD', 'LFD vs. HFD'))

    plt.show()

# deg_venn_diagram()

def plot_custom_corr_heatmap():
    # colocalization analysis performed based on https://www.nature.com/articles/s41467-021-21892-z#Sec8
    sample_type = "visium"
    meta = utils.get_meta_data(sample_type)
    ref_ct= "B cells-1"
    cell_types = ["IgA plasma cells-1", "Neutrophils"]
    str_cell_types = "_".join(cell_types)
    lst_corr = []
    lst_cond = []
    dict_colocalization = dict()
    for ind, row in meta.iterrows():
        lst_cell_types = []
        sample_id = row["sample_id"]
        condition = row["condition"]
        df_abundance = pd.read_csv(f"../data/out_data/cell2location_map/cell_type_abundances_{sample_id}_filtered_deconv_15_20.csv", index_col=0)
        for ct in df_abundance.columns:
                ct = ct.split("_")[-1]
                lst_cell_types.append(ct)
        df_abundance.columns = lst_cell_types
        
        corr = df_abundance.corr(method='pearson', min_periods=1)
        # print(corr)
        lst_cond.append(condition)
        lst_corr.append(list(corr[cell_types].loc[[ref_ct]].values[0]))

    pd_results = pd.DataFrame(lst_corr, columns=cell_types, index=lst_cond)
    print(pd_results)
    
    dict_cell_type =dict()

    for c_t in cell_types:
        dict_cell_type[c_t] = []


    for ind, row in pd_results.iterrows():
        diet = ind.split("-")[0]
        
        for c_t in cell_types:
            
            if "no-AOM-DSS" in ind:
                dict_cell_type[c_t].append(["no-AOM-DSS", diet, row[c_t]])
            else:
                dict_cell_type[c_t].append(["AOM-DSS", diet, row[c_t]])

    print(dict_cell_type)

    for c_t in cell_types:
        print(c_t)
        df = pd.DataFrame(dict_cell_type[c_t], columns=['Group', "Diet", "Score"])
        print(df)
        # plot with seaborn barplot
        sns.barplot(data=df, x='Diet', y='Score', hue='Group').set(title=f'{ref_ct} vs. {c_t}')
        plt.tight_layout()
        plt.rcParams['figure.dpi']= 300
        plt.rcParams['figure.figsize']= (240, 120)
        plt.savefig(f"../plots/vis_deconvolution/barplot_{ref_ct}_vs_{c_t}.pdf")
        plt.clf()
    
    pd_results = pd_results.transpose()
    # print(pd_results)
    # cmap = sns.diverging_palette(230, 20, as_cmap=True)
    sns.heatmap(pd_results, annot=False, xticklabels=True, yticklabels=True)
    plt.tight_layout()
    plt.rcParams['figure.dpi']= 300
    plt.rcParams['figure.figsize']= (240, 120)
    plt.savefig(f"../plots/vis_deconvolution/heatmap_{ref_ct}_vs_{str_cell_types}.pdf")
    plt.close()




# plot_custom_corr_heatmap()

def plot_dotplot_tumor_markers_vs_conditions():
    adata_integ_clust = sc.read_h5ad("/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/sc_epicells_integrated_clustered.h5ad")
    df_over_threshold = pd.read_csv(f"/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/sc_epicells_tumor_markers_combined_thr_5.csv")
    
    adata = utils.get_filtered_concat_data("sc_epicells")
    adata = adata[adata_integ_clust.obs_names,:]
    adata.obs["leiden_0.20"] = adata_integ_clust.obs["leiden_0.20"]
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    print(df_over_threshold["Over threshold"].values)
    adata.obs["tumor_markers"] = df_over_threshold["Over threshold"].values
    # adata = adata[adata.obs["tumor_markers"]==1,:]
    # filter out the cells missing in adata_integ_clust
    # 
    # print(adata[adata.obs["tumor_markers"].isin([1]),:])
    
    print(adata)
    adata.var.index = pd.Index(gen.upper() for gen in adata.var.index.values)
    # markers = ["WIF1", "AXIN2", "NKD1", "NOTUM", "MMP7", "PROX1", "SOX4"]
    markers = ["WIF1", "NKD1", "NOTUM", "MMP7", "PROX1"]
    # sc.pl.dotplot(adata, markers, groupby='condition', standard_scale='group', dendrogram=False, show=False)
    # sc.pl.dotplot(adata, markers, groupby='condition', standard_scale='var', dendrogram=False, show=False)
    sc.pl.dotplot(adata, markers, groupby='condition', dendrogram=True, show=False, save="epicells_tumor_markers_no-AXIN1-no-SOX4.pdf")
    sc.pl.dotplot(adata, markers, groupby='condition', standard_scale='var', dendrogram=True, show=False, save="epicells_tumor_markers_ss_var-no-AXIN1-no-SOX4.pdf")
    sc.pl.dotplot(adata, markers, groupby='condition', standard_scale='group', dendrogram=True, show=False, save="epicells_tumor_markers_ss_group-no-AXIN1-no-SOX4.pdf")
    """sc.pl.dotplot(adata, markers, groupby='condition', dendrogram=True, show=False, save="epicells_tumor_markers.pdf")
    sc.pl.dotplot(adata, markers, groupby='condition', standard_scale='var', dendrogram=True, show=False, save="epicells_tumor_markers_ss_var.pdf")
    sc.pl.dotplot(adata, markers, groupby='condition', standard_scale='group', dendrogram=True, show=False, save="epicells_tumor_markers_ss_group.pdf")"""
    # plt.tight_layout()
    # plt.savefig("../plots/sc_epi_cells_aom_noaom_visualize_markers/dot_plot.pdf")
 
# plot_dotplot_tumor_markers_vs_conditions()


def plot_dotplot_tumor_markers_vs_individual_condition(condition, group_by):
    adata_integ_clust = sc.read_h5ad("/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/sc_epicells_aom_noaom_integrated_clustered.h5ad")
    adata = utils.get_filtered_concat_data("sc_epicells_aom_noaom")
    adata_integ_clust = adata_integ_clust[adata_integ_clust.obs["condition"]==condition,:]
    adata = adata[adata_integ_clust.obs_names,:]
    adata.obs["leiden_0.20"] = adata_integ_clust.obs["leiden_0.20"]
    adata_integ_clust.obs["condition"][condition] = condition.split("_")[0]
    condition = condition.split("_")[0]

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    adata_integ_clust.var.index = pd.Index(gen.upper() for gen in adata_integ_clust.var.index.values)
    adata.var.index = pd.Index(gen.upper() for gen in adata.var.index.values)
    # markers = ["EPCAM", "CDH1", "MUC3", "CDHR5", "WNT6", "WNT10A", "FZD10", "DKK3", "WIF1", "NKD1", "AXIN2", "NOTUM", "SOX4", "PROX1", "MUC2", "REG4", "CCL9", "MMP7",  "IFITM3"]
    markers = [ "EPCAM", "CDH1", "CDHR5", "MUC3", "MUC2", "REG4", "AXIN2", "SOX4", "IFITM3", "WNT6", "WNT10A", "FZD10", "DKK3", "WIF1", "NKD1", "NOTUM", "PROX1", "MMP7"]
    # markers = ["EPCAM", "MUC3", "CDHR5", "WNT6", "WNT10A", "DKK3", "WIF1", "NKD1", "NOTUM", "SOX4", "PROX1", "MUC2", "REG4", "CCL9", "MMP7",  "IFITM3"]

    sc.pl.dotplot(adata, markers, groupby=group_by, dendrogram=False, vmax=10, show=False, save=f"epicells_{condition}.pdf")
    # sc.pl.dotplot(adata, markers, groupby=group_by, dendrogram=False, standard_scale='var', show=False, save=f"epicells_{condition}.pdf")
    
    # sc.pl.dotplot(adata, markers, groupby=group_by, standard_scale='var', dendrogram=True, show=False, save=f"epicells_{condition}_ss_var.pdf")
    # sc.pl.dotplot(adata, markers, groupby=group_by, standard_scale='group', dendrogram=True, show=False, save=f"epicells_{condition}_ss_group.pdf")


"""plot_dotplot_tumor_markers_vs_individual_condition("CD-AOM-DSS-Epi_plus_DN", "leiden_0.20")
plot_dotplot_tumor_markers_vs_individual_condition("Control-no-AOM-DSS-Immune", "leiden_0.20")
plot_dotplot_tumor_markers_vs_individual_condition("HFD-AOM-DSS-Epi_plus_DN", "leiden_0.20")
plot_dotplot_tumor_markers_vs_individual_condition("LFD-AOM-DSS-Epi_plus_DN", "leiden_0.20")"""


def plot_dotplot_cancer_vs_control(group_by):
    adata_integ_clust = sc.read_h5ad("/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/sc_epicells_integrated_clustered.h5ad")
    adata = utils.get_filtered_concat_data("sc_epicells_aom_noaom")
    adata = adata[adata_integ_clust.obs_names,:]
    adata.obs["leiden_0.20"] = adata_integ_clust.obs["leiden_0.20"]
    adata.obs["crc"] = "no-AOM-DSS"
    adata.obs.loc[adata[adata.obs["condition"].isin(["CD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN"]),:].obs_names, 'crc'] = "AOM-DSS"
    
    """for cond in ["CD-AOM-DSS-Epi_plus_DN", "Control-no-AOM-DSS-Immune", "HFD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN"]:
        
        if "Immune" in cond:
            adata[adata.obs["condition"]==cond,: ].obs["crc"] = "no-AOM-DSS"
        else:
            adata[adata.obs["condition"]==cond,: ].obs["crc"] = "AOM-DSS"
        # adata.obs["condition"]["Control-no-AOM-DSS-Immune"] = "No-AOM-DSS"
        # adata.obs["condition"]["HFD-AOM-DSS-Epi_plus_DN"] = "AOM-DSS"
        # adata.obs["condition"]["LFD-AOM-DSS-Epi_plus_DN"] = "AOM-DSS"
    """

    print(adata.obs["crc"])
    print(adata)

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    adata_integ_clust.var.index = pd.Index(gen.upper() for gen in adata_integ_clust.var.index.values)
    adata.var.index = pd.Index(gen.upper() for gen in adata.var.index.values)
    # markers = ["EPCAM", "CDH1", "MUC3", "CDHR5", "WNT6", "WNT10A", "FZD10", "DKK3", "WIF1", "NKD1", "AXIN2", "NOTUM", "SOX4", "PROX1", "MUC2", "REG4", "CCL9", "MMP7",  "IFITM3"]
    # markers = ["WNT6", "WNT10A", "FZD10", "DKK3", "WIF1", "NKD1", "NOTUM", "PROX1", "MMP7"] 
    markers = [ "EPCAM", "CDH1", "CDHR5", "MUC3", "MUC2", "REG4", "AXIN2", "SOX4", "IFITM3", "WNT6", "WNT10A", "FZD10", "DKK3", "WIF1", "NKD1", "NOTUM", "PROX1", "MMP7"]
    # markers = ["EPCAM", "MUC3", "CDHR5", "WNT6", "WNT10A", "DKK3", "WIF1", "NKD1", "NOTUM", "SOX4", "PROX1", "MUC2", "REG4", "CCL9", "MMP7",  "IFITM3"]

    sc.pl.dotplot(adata, markers, groupby=group_by, dendrogram=False, show=False, save=f"epicells_aom-dss_no-aom-dss.pdf")

    markers = ["WNT6", "WNT10A", "FZD10", "DKK3", "WIF1", "NKD1", "NOTUM", "PROX1", "MMP7"] 
    sc.pl.dotplot(adata, markers, groupby=group_by, dendrogram=False, show=False, save=f"epicells_aom-dss_no-aom-dss_excluding_HEGs.pdf")    
    
    top10_16_clust = ["NKD1", "ASPSCR1", "IFITM3", "NOTUM", "PROX1", "CYBA", "WIF1", "FAM89A", "MT3", "IFITM1"]
    sc.pl.dotplot(adata, top10_16_clust, groupby=group_by, dendrogram=False, show=False, save=f"epicells_aom-dss_no-aom-dss_clust_16_top10_DEGs.pdf")  
    # sc.pl.dotplot(adata, markers, groupby=group_by, standard_scale='var', dendrogram=True, show=False, save=f"epicells_{condition}_ss_var.pdf")
    # sc.pl.dotplot(adata, markers, groupby=group_by, standard_scale='group', dendrogram=True, show=False, save=f"epicells_{condition}_ss_group.pdf")

# plot_dotplot_cancer_vs_control("crc")


def plot_barplot_cell_proportion_in_cluster(adata_path, sample_type, obs_cluster):
    # colocalization analysis performed based on https://www.nature.com/articles/s41467-021-21892-z#Sec8
    adata = sc.read_h5ad(adata_path)

    meta = utils.get_meta_data(sample_type)
    condition = list(np.unique(meta['condition']))
    condition = ["Control-no-AOM-DSS-Immune", "CD-AOM-DSS-Epi_plus_DN", "LFD-AOM-DSS-Epi_plus_DN", "HFD-AOM-DSS-Epi_plus_DN"]
    prop_lst = []


    for clust in adata.obs[obs_cluster].cat.categories:
        adata_tmp = adata[adata.obs[obs_cluster]==clust,:]
        n_of_cells = adata_tmp.shape[0]
        # prop_dict[clust] = []
        
        for cond in condition:
            adata_tmp2 = adata_tmp[adata_tmp.obs["condition"]==cond,:]
            proportion = adata_tmp2.shape[0]/n_of_cells

            print(clust, cond, n_of_cells, adata_tmp2.shape)    
            cond_lbl = cond.split("_")[0]
            if "Control" in cond:
                cond_lbl = "no-AOM-DSS"
            
            prop_lst.append([cond_lbl, clust, proportion])
    




    df = pd.DataFrame(prop_lst, columns=['Group', "Cluster", "Proportion"])
    sns.barplot(data=df, x='Cluster', y='Proportion', hue='Group').set(title='Proportions of Cells in Clusters')
    plt.tight_layout()
    plt.rcParams['figure.dpi']= 300
    plt.rcParams['figure.figsize']= (240, 120)
    plt.savefig(f"../plots/../plots/sc_epicells_aom_noaom_cluster/barplot_proportion_cluster.pdf") 
    plt.clf()

# plot_barplot_cell_proportion_in_cluster("/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/sc_epicells_integrated_clustered.h5ad", "sc_epicells_aom_noaom", "leiden_0.20")
# plot_barplot_cell_proportion_in_cluster("/Users/ahmet/Google Drive/Projects/saezlab/CRCDiet/data/out_data/sc_epicells_aom_noaom_integrated_clustered.h5ad", "sc_epicells_aom_noaom", "leiden_0.20")

def plot_barplot_markers_vs_groups(adata, markers, group_by):

    expression_cutoff = 0.0
    markers = ['CD3D', 'CD79A','CST3']
    dp = sc.pl.DotPlot(adata, markers, groupby=group_by)
    obs_bool = dp.obs_tidy > expression_cutoff
    dot_color_df = (dp.obs_tidy.mask(~obs_bool).groupby(level=0).mean().fillna(0))