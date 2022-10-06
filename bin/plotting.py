from cProfile import label
from turtle import color
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib import cm
import scanpy as sc
import utils
import numpy as np
import matplotlib

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
    
def proj(x, y, ax):
    for label in np.unique(y):
        m = y == label
        ax.scatter(x[m,0], x[m,1], label=label)
        ax.legend()
        
        
def hextofloats(h):
    '''Takes a hex rgb string (e.g. #ffffff) and returns an RGB tuple (float, float, float).'''
    return tuple(int(h[i:i + 2], 16) / 255. for i in (1, 3, 5))
        
    
def varshift_condition(df, name):
    # Define fig
    fig, ax = plt.subplots(1,1, figsize=(3,3), dpi=150, facecolor='white')
    
    # Order by median
    order = df.groupby('condition').median().sort_values('dist', ascending=False).index

    # Plot boxplot
    sns.boxplot(data=df, x='condition', y='dist', ax=ax, order=order, palette=cond_colors, linewidth=0.5, fliersize=0)
    ax.axhline(1, ls='--', c='black')
    ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
    ax.set_title('Variation shifts', fontsize=11)
    ax.set_ylabel('inter / intra dist')
    ax.set_xlabel('')

    # Save
    fig.savefig('../plots/{0}'.format(name), bbox_inches='tight')
    

def varshift_cell_type(df, table, name, order, hue_order):
    # Define fig
    fig, ax = plt.subplots(1,1, figsize=(9,3), tight_layout=True, dpi=150, facecolor='white')
    
    # Order by median
    #order = df.groupby('cell_type').median().sort_values('dist', ascending=False).index

    # Plot boxplot
    sns.boxplot(data=df, x='cell_type', y='dist', hue='condition', ax=ax, palette=cond_colors, linewidth=0.5, 
                fliersize=0, hue_order=hue_order, order=order)
    ax.axhline(1, ls='--', c='black')
    ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
    ax.set_title('Cell type variation shifts', fontsize=11)
    ax.set_xlabel('')
    ax.set_ylabel('inter dist / intra dist')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    
    n = len(cond_colors)
    cell_types = np.repeat(order, n)
    for i,patch in enumerate(ax.artists):
        hue = hue_order[i % n]
        ctype = cell_types[i]
        adj_pval = table.loc[hue].loc[ctype]['adj_pval']
        if adj_pval < 0.05:
            a = 1
        else:
            a = 0.33
        r, g, b = hextofloats(cond_colors[hue])
        patch.set_facecolor((r, g, b, a))
    

    # Save
    fig.savefig('../plots/{0}'.format(name), bbox_inches='tight')