# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/06/19
content:    Explore merged datasets
'''
import os
import sys
import numpy as np
import scipy as sp
import pandas as pd
import loompy
import argparse

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
import singlet


data_fdn = '../data/all_viruses/'


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--normalize', action='store_true')
    args = pa.parse_args()

    if args.normalize:
        fn = data_fdn+'correlations_vRNA_gene_expression_normalized_2+viral_reads.tsv'
    else:
        fn = data_fdn+'correlations_vRNA_gene_expression_notnormalized_2+viral_reads.tsv'
    if not os.path.isfile(fn):
        print('Load merged data from loom file')
        ds = singlet.Dataset(dataset='all')

        print('Calculate correlations')
        corrs = []
        viruses = ds.samplesheet['virus'].unique()
        for virus in viruses:
            print(virus)
            dsi = ds.query_samples_by_metadata('(virus == @virus) & (n_viral_reads >= 2)', local_dict=locals())
            print(dsi.n_samples)

            if args.normalize:
                dsi.samplesheet['coverage'] = dsi.counts.sum(axis=0)
                normg = 1e6 * dsi.samplesheet['coverage']
                normv = 1e6 * (dsi.samplesheet['coverage'] + dsi.samplesheet['n_viral_reads'])
                dsi.counts.loc[:, :] = dsi.counts / normg
                dsi.samplesheet['viral_reads_per_million'] = dsi.samplesheet['n_viral_reads'] / normv

                corr = dsi.correlation.correlate_features_phenotypes('viral_reads_per_million').fillna(0)
            else:
                corr = dsi.correlation.correlate_features_phenotypes('n_viral_reads').fillna(0)

            corr.name = virus
            corrs.append(corr)
        corrs = pd.DataFrame(corrs).T
        corrs.to_csv(fn, sep='\t', index=True)
    else:
        print('Load correlations from file')
        corrs = pd.read_csv(fn, sep='\t', index_col=0)

    print('Rename viruses')
    corrs.rename(
            columns={'dengue': 'DENV', 'zika': 'ZIKV', 'veev': 'VEEV', 'wnv': 'WNV', 'flu': 'influenza'},
        inplace=True)

    print('Keep only VEEV, DENV, ZIKV, FLU, WNV')
    corrs = corrs[['DENV', 'ZIKV', 'VEEV', 'WNV', 'influenza']]

    print('Calculate ranks of the correlation coefficients')
    corranks = corrs.copy()
    ci = corrs.index
    for col in corranks.columns:
        c = corranks[[col]].sort_values(by=col, ascending=False)
        c['rank'] = np.arange(c.shape[0]) + 1
        c = c.loc[ci]
        corranks[col] = c['rank']

    print('Dimensionality reduction')
    ntop = 200
    ind = (corranks <= ntop).any(axis=1) | (corranks >= corranks.shape[0] - ntop).any(axis=1)
    data = corrs.values[ind]

    # Rescale between -1 and 1 within each virus
    data = (data - data.min(axis=0)) / (data.max(axis=0) - data.min(axis=0)) * 2 - 1

    index = corrs.index[ind]

    print('PCA 2 dimensional')
    from sklearn.decomposition.pca import PCA
    model = PCA(n_components=2)
    res = model.fit_transform(data)
    pcs = pd.DataFrame(data=res, index=index, columns=['pc1', 'pc2'])

    if False:
        print('PCA 4 dimensional')
        from sklearn.decomposition.pca import PCA
        model = PCA(n_components=4)
        pcsh = model.fit_transform(data)

        print('t-SNE')
        from sklearn.manifold.t_sne import TSNE
        perplexity = 20
        model = TSNE(perplexity=perplexity)
        res = model.fit_transform(pcsh)
        ts = pd.DataFrame(data=res, index=index, columns=['tsne1', 'tsne2'])

        print('Knn graph')
        from scipy.spatial.distance import pdist, squareform
        kk = 15
        thr = 0.3
        #dist = squareform(pdist(ts.values, 'euclidean'))
        dist = squareform(pdist(ts.values, 'euclidean'))
        edges = set()
        simi = dist.max()-dist
        for ir, row in enumerate(simi):
            neis = np.argpartition(row, -kk)[-kk:]
            neis = neis[simi[ir, neis] >= thr]
            for nei in neis:
                if nei == ir:
                    continue
                edges.add(frozenset([ir, nei]))
        edges = [tuple(e) for e in edges]

        print('Leiden clusters')
        import igraph as ig
        import leidenalg
        g = ig.Graph(edges)
        partition = leidenalg.CPMVertexPartition(
                g,
                resolution_parameter=0.001,
                )
        leidenalg.Optimiser().optimise_partition(partition)
        comms = partition.membership
        clusters = pd.Series(comms, index=index)

        df = ts.copy()
        df['cluster'] = clusters
        df.to_csv('../data/all_viruses/gene_clusters_compare5.tsv', sep='\t', index=True)

    else:
        df = pd.read_csv('../data/all_viruses/gene_clusters_compare5.tsv', sep='\t', index_col=0)
        ts = df[['tsne1', 'tsne2']]
        clusters = df['cluster']

    n_clusters = len(clusters.unique())

    sys.exit()

    print('Viral phylogeny')
    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage, leaves_list
    bins = [-1, -0.3, -0.1, 0.1, 0.3, 1]
    datap = np.zeros_like(data.T)
    for b in bins[:-1]:
        datap[data.T >= b] += 1
    dv = pdist(datap, metric='euclidean')
    zv = linkage(dv, 'average', optimal_ordering=True)
    indv = leaves_list(zv)

    print('Plot dimensionality reduction')
    import matplotlib.gridspec as gridspec
    vecs = {'PCA': pcs, 'tsne': ts}

    vs = vecs['tsne']

    fig = plt.figure(figsize=(12, 4))
    axs = []
    gs = gridspec.GridSpec(2, 6)
    axs.append(fig.add_subplot(gs[0, 0]))
    axs.append(fig.add_subplot(gs[0, 1], sharex=axs[0], sharey=axs[0]))
    axs.append(fig.add_subplot(gs[0, 2], sharex=axs[0], sharey=axs[0]))
    axs.append(fig.add_subplot(gs[1, 0], sharex=axs[0], sharey=axs[0]))
    axs.append(fig.add_subplot(gs[1, 1], sharex=axs[0], sharey=axs[0]))
    axs.append(fig.add_subplot(gs[1, 2], sharex=axs[0], sharey=axs[0]))
    axs.append(fig.add_subplot(gs[:, 3]))
    axs.append(fig.add_subplot(gs[:, 4:]))
    viruses = corrs.columns
    plotnames = ['DENV', 'ZIKV', 'WNV', 'VEEV', 'influenza'] + ['Gene cluster']
    for plotname, ax in zip(plotnames, axs):
        kwargs = {}
        if plotname in viruses:
            c = corrs[plotname].loc[ind].values
            c = (c - c.min()) / (c.max() - c.min())
            cmap = sns.diverging_palette(250, 15, s=75, l=40, center="dark", as_cmap=True)
            kwargs['cmap'] = cmap
        elif plotname == 'Gene cluster':
            palette = sns.color_palette('Set1', n_colors=n_clusters)
            c = [palette[x] for x in clusters.values]
            for cn in np.arange(n_clusters):
                xm = vs.iloc[:, 0].loc[clusters == cn].mean()
                ym = vs.iloc[:, 1].loc[clusters == cn].mean()
                ax.text(xm, ym, str(cn), ha='center', va='center')

        ax.scatter(vs.iloc[:, 0], vs.iloc[:, 1], s=30, c=c, alpha=0.15, **kwargs)
        ax.set_title(plotname)
        ax.set_axis_off()
    fig.text(0.25, 0.04, 'dimension 1', va='bottom', ha='center')
    fig.text(0.01, 0.52, 'dimension 2', va='center', ha='left', rotation=90)

    # Heatmap with averages within clusters
    heatmap = pd.DataFrame(
            data=np.zeros((n_clusters, len(viruses))),
            index=pd.Index(np.arange(n_clusters), name='cluster'),
            columns=corrs.columns,
            )
    for v in viruses:
        for cn in np.arange(n_clusters):
            val = corrs.loc[ind, v].loc[clusters == cn].mean()
            heatmap.loc[cn, v] = val
    vmax = np.abs(heatmap.values).max()

    # Normalize within each virus
    heatmap_norm = (heatmap - heatmap.min(axis=0)) / (heatmap.max(axis=0) - heatmap.min(axis=0)) * 2 - 1
    vmax = 1

    # Double hierarchical clustering
    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage, leaves_list
    d0 = pdist(heatmap_norm.values, metric='euclidean')
    z0 = linkage(d0, 'average', optimal_ordering=True)
    ind0 = leaves_list(z0)
    #d1 = pdist(heatmap_norm.values.T, metric='euclidean')
    #z1 = linkage(d1, 'average', optimal_ordering=True)
    #ind1 = leaves_list(z1)
    z1, ind1 = zv, indv
    heatmap_hh = heatmap_norm.iloc[ind0].T.iloc[ind1].T

    # Plot dendrogram
    from scipy.cluster.hierarchy import dendrogram
    dendrogram(z1, orientation='left', ax=axs[-2])
    axs[-2].set_axis_off()
    axs[-2].set_ylim(*(axs[-2].get_ylim())[::-1])

    ax = axs[-1]
    cmap = sns.diverging_palette(250, 15, s=75, l=40, center="dark", as_cmap=True)
    # Choose order to be consistent with dengrogram
    sns.heatmap(
            heatmap_hh.T.iloc[::-1], ax=ax, vmin=-vmax, vmax=vmax, cmap=cmap,
            xticklabels=True)
    for tk in ax.get_yticklabels():
        tk.set_rotation(0)
    ax.set_title('Rescaled correlation between\ngene expression and vRNA')
    ax.set_xlabel('Gene cluster')
    ax.set_ylim(0, 5)
    ax.set_xlim(*(ax.get_xlim()[::-1]))

    cax = fig.get_axes()[-1]

    fig.tight_layout(rect=(0.03, 0.012, 0.98, 1), w_pad=0.3, h_pad=0.5)
    #fig.savefig('../figures/compare_5_viruses.png', dpi=600)
    #fig.savefig('../figures/compare_5_viruses.pdf')
    #fig.savefig('../figures/compare_5_viruses.svg')

    if False:
        print('Gene Ontology analysis')
        gene_clusters = {key: val.index.tolist() for key, val in clusters.groupby(clusters)}
        for i, gclu in gene_clusters.items():
            with open('../data/all_viruses/gene_cluster_{:}.csv'.format(i), 'wt') as f:
                f.write('\n'.join(gclu))
        print('A bunch of clicks online...')

    if False:
        go_res = []
        for i, gclu in gene_clusters.items():
            go_resi = pd.read_csv(
                    '../data/all/GO_analysis_cluster_{:}.txt'.format(i),
                    sep='\t',
                    skiprows=11,
                    )
            cols = []
            for icol, col in enumerate(go_resi.columns):
                if icol == 2:
                    cols.append('observed')
                    continue
                if col.startswith('go_cluster_{:}_input.csv'.format(i)):
                    col = col[len('go_cluster_{:}_input.csv '.format(i)):]
                if col.startswith('(') and col.endswith(')'):
                    col = col.strip('()')
                cols.append(col)
            go_resi.columns = cols
            go_resi['cluster'] = i
            ii = go_resi.index
            go_resi.sort_values('FDR', inplace=True)
            go_resi['rank'] = np.arange(len(go_resi)) + 1
            go_resi = go_resi.loc[ii]
            go_res.append(go_resi)
        go_res = pd.concat(go_res)
        go_res_filt = go_res.loc[go_res['rank'] < 30]
        fdr = go_res_filt[
                ['GO biological process complete', 'FDR', 'cluster']
                ].set_index(
                    ['GO biological process complete', 'cluster'],
                    )['FDR'].unstack().fillna(1)
        nlfdr = -np.log10(fdr)

        #fig, ax = plt.subplots(figsize=(8, 25))
        #sns.heatmap(
        #        np.log2(nlfdr + 0.1),
        #        vmin=0,
        #        vmax=np.log2(nlfdr.values.max() + 0.1),
        #        cmap='plasma',
        #        ax=ax,
        #        xticklabels=True,
        #        yticklabels=True)
        #for tk in ax.get_xticklabels():
        #    tk.set_rotation(90)
        #fig.tight_layout(rect=(0, 0.1, 1, 1))

        data = np.log2(nlfdr + 0.1).loc[:, [ii for ii in ind0 if ii in nlfdr.columns]]
        data = data.loc[data.max(axis=1) > 1]
        d00 = pdist(data.values, metric='correlation')
        z00 = linkage(d00, 'average', optimal_ordering=True)
        ind00 = leaves_list(z00)
        g = sns.clustermap(
                data,
                vmin=0,
                vmax=np.log2(nlfdr.values.max() + 0.1),
                cmap='inferno_r',
                xticklabels=True,
                yticklabels=True,
                col_cluster=False,
                row_linkage=z00,
                )
        for tk in g.ax_heatmap.get_yticklabels():
            tk.set_fontsize(8)
        g.fig.set_size_inches(6, 20)
        g.fig.subplots_adjust(bottom=0.03, top=0.98, right=0.65, left=0.02)

    #TODO: check that SLITs/ROBOs are not specifically expressed by U87
    fn = '../data/all_viruses/metascape_result_cluster6.xlsx'
    col_id = 'Input ID'
    col_slro = 'R-HSA-9010553 Regulation of expression of SL'
    dfex = pd.read_excel(fn, 'Annotation')
    genes_slro = dfex[[col_id, col_slro]].set_index(col_id).iloc[:, 0]
    genes_slro = genes_slro.index[genes_slro == 1]

    ds = singlet.Dataset(dataset='all')
    ds.query_features_by_name(genes_slro, inplace=True)

    dsu87 = ds.query_samples_by_metadata('(virus == "veev") & (n_viral_reads < 1)')
    dshuh7 = ds.query_samples_by_metadata('(virus == "dengue") & (n_viral_reads < 1)')
    dsflu = ds.query_samples_by_metadata('(virus == "flu") & (n_viral_reads < 1)')
    stu87 = dsu87.counts.get_statistics(metrics=('mean', 'std'))
    sthuh7 = dshuh7.counts.get_statistics(metrics=('mean', 'std'))
    stflu = dsflu.counts.get_statistics(metrics=('mean', 'std'))
    rel_exp = sthuh7['mean'] / stu87['mean']
    print(rel_exp.quantile([0, 0.25, 0.5, 0.75, 1]))
    rel_exp = stflu['mean'] / stu87['mean']
    print(rel_exp.quantile([0, 0.25, 0.5, 0.75, 1]))



    # Yes, these are genes that have higher baseline expression in U87



    plt.ion()
    plt.show()
