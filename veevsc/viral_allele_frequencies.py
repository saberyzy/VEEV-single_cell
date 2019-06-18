# vim: fdm=indent
'''
author:     Fabio Zanini
date:       09/05/19
content:    Check viral allele frequencies
'''
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
import singlet


if __name__ == '__main__':

    print('Loading mRNA counts...')
    counts = singlet.CountsTable(pd.read_csv(
            '../data/dataset/VEEVcounts.tsv',
            sep='\t',
            index_col=0,
            ).astype(np.float32),
            )
    counts._spikeins = list(counts.index[-103:-6])
    counts._otherfeatures = list(counts.index[-6:])
    print('Counts loaded.')

    print('Load cell metadata')
    samplesheet = singlet.SampleSheet(pd.read_csv(
            '../data/dataset/VEEVcellmetadata.tsv',
            sep='\t',
            index_col=0),
            )
    print('Done')

    ds = singlet.Dataset(
            counts_table=counts,
            samplesheet=samplesheet,
            )

    print('Load viral allele frequencies')
    fn = '../data/dataset/viral_allele_frequencies_100000_reads.cdf'
    data = xr.open_dataset(fn)

    la, n, l = data['allele_frequencies'].shape

    cov = data['coverage']
    covm = np.log10(np.maximum(9e-1, cov.mean(axis=1).data))

    print('Focus on cells with average coverage above 100')
    ind = (covm > 2).nonzero()[0]
    afs = data['allele_frequencies'][:, ind]
    dsh = ds.query_samples_by_name(afs['sample'].data)
    covh = cov[ind]

    # Use N if no coverage
    afs[5] = 1e-2
    maj = afs.argmax(axis=0)

    # Find non-reference major alleles
    maj_snv = data.coords['nucleotide'].data[maj] != data['reference'].data
    maj_snv &= data.coords['nucleotide'].data[maj] != 'N'
    maj_snv = pd.DataFrame(
            maj_snv,
            index=data.sample.data[ind],
            columns=data.position,
            )

    print('Plot clustermaps')
    maj_snv_plot = maj_snv.loc[:, maj_snv.sum(axis=0) != 0]
    pos_poly = maj_snv.columns[maj_snv.sum(axis=0) != 0]
    alphal = list(data.nucleotide.data)
    refi = np.array([alphal.index(x) for x in data.reference], int)
    afs_poly = afs[:, :, pos_poly]
    refi_poly = refi[pos_poly]
    covh_poly = covh[:, pos_poly]

    afs_nonref_poly = np.zeros((afs_poly.shape[1], afs_poly.shape[2]))
    for j, refij in enumerate(refi_poly):
        afs_nonref_poly[:, j] = afs_poly[refij, :, j]
    afs_nonref_poly = 1 - afs_nonref_poly
    # Default missing coverage to majority allele in the population
    afs_default = np.zeros((afs_poly.shape[1], afs_poly.shape[2]))
    afs_default[:] = afs_nonref_poly.mean(axis=0) > 0.5
    #ind_noncov = (afs_poly.sum(axis=0) < 1).data
    ind_noncov = (covh_poly < 10).data
    afs_nonref_poly[ind_noncov] = afs_default[ind_noncov]

    afs_nonref_poly = pd.DataFrame(
            afs_nonref_poly,
            index=data.sample.data[ind],
            columns=pos_poly,
            )

    # Exclude some cells with spotty coverage (manual annotation)
    #exclude_spotty = [
    #    '1001703011_P15', '1001703009_P24', '1001703011_L22', '1001703007_I20',
    #    '1001703011_O18', '1001703009_O17', '1001703011_F11', '1001703011_F9',
    #    '1001703011_H19', '1001703011_J14', '1001703009_E24', '1001703007_N20']
    #afs_nonref_poly.drop(exclude_spotty, axis=0, inplace=True)

    # Cell colors by total amount of reads
    v = dsh.samplesheet['log_virus_reads_per_million'].fillna(-1).values
    v = (v - v.min()) / (v.max() - v.min())
    cmap = sns.dark_palette('Blue', as_cmap=True)
    colors_viral = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)

    v = dsh.samplesheet['time_index'].values
    v = (v - v.min()) / (v.max() - v.min())
    cmap = sns.dark_palette('Green', as_cmap=True)
    colors_time = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)

    v = dsh.samplesheet['MOI'].values
    v[v>0] += 1
    v = (v - v.min()) / (v.max() - v.min())
    cmap = sns.dark_palette('Orange', as_cmap=True)
    colors_moi = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)

    colors = pd.concat([colors_viral, colors_time, colors_moi], axis=1)
    colors.columns = ['vRNA', 'time', 'MOI']

    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage, leaves_list
    D = pdist(afs_nonref_poly.values, 'euclidean')
    Z = linkage(D, 'average', optimal_ordering=True)

    # Sort cells by vRNA
    ind_sort = dsh.samplesheet['log_virus_reads_per_million'].fillna(-1).values.argsort()

    # Allele frequencies
    g = sns.clustermap(
        afs_nonref_poly.iloc[ind_sort],
        figsize=(18, 8),
        xticklabels=True,
        yticklabels=True,
        row_colors=colors.iloc[ind_sort],
        row_linkage=Z,
        col_cluster=False,
        row_cluster=False,
        )
    plt.subplots_adjust(0.01, 0.08, 0.92, 0.95)
    ax = g.ax_col_dendrogram
    ax.set_xlim(0, afs_nonref_poly.shape[1])
    i = afs_nonref_poly.columns.tolist().index(9226)
    ax.arrow(i + 0.5, 0.3, 0, -0.25, ec='red', fc='red', head_width=2, head_length=0.15, length_includes_head=True)
    ax.text(i + 0.5, 0.35, '9227 A>G, Ala>Thr, E3 protein', ha='center', va='bottom', fontsize=10)
    ax1 = g.fig.get_axes()[-1]
    ax1.set_ylabel('Frequency of\nalternate allele', rotation=0, labelpad=15, va='center', ha='left')
    g.ax_heatmap.set_yticklabels(['' for i in range(afs_nonref_poly.shape[0])])
    g.ax_heatmap.set_position([0.07, 0.10, 0.91, 0.76])
    g.ax_col_dendrogram.set_position([0.07, 0.86, 0.91, 0.1])
    g.ax_row_colors.set_position([0.01, 0.10, 0.05, 0.76])
    g.fig.get_axes()[-1].set_position([0.01, 0.88, 0.05, 0.1])
    g.ax_heatmap.set_ylabel('cells')

    # Coverage
    g = sns.clustermap(
        np.log10(0.1 + covh_poly.to_pandas().iloc[ind_sort]),
        figsize=(18, 8),
        xticklabels=True,
        yticklabels=True,
        row_colors=colors.iloc[ind_sort],
        row_linkage=Z,
        col_cluster=False,
        row_cluster=False,
        )
    plt.subplots_adjust(0.01, 0.08, 0.92, 0.95)
    ax = g.ax_col_dendrogram
    ax.set_xlim(0, afs_nonref_poly.shape[1])
    i = afs_nonref_poly.columns.tolist().index(9226)
    ax.arrow(i + 0.5, 0.3, 0, -0.25, ec='red', fc='red', head_width=2, head_length=0.15, length_includes_head=True)
    ax.text(i + 0.5, 0.35, '9227 A>G, Ala>Thr, E3 protein', ha='center', va='bottom', fontsize=10)
    ax1 = g.fig.get_axes()[-1]
    ax1.set_ylabel('Coverage', rotation=0, labelpad=15, va='center', ha='left')
    g.ax_heatmap.set_yticklabels(['' for i in range(afs_nonref_poly.shape[0])])
    g.ax_heatmap.set_position([0.07, 0.10, 0.91, 0.76])
    g.ax_col_dendrogram.set_position([0.07, 0.86, 0.91, 0.1])
    g.ax_row_colors.set_position([0.01, 0.10, 0.05, 0.76])
    g.fig.get_axes()[-1].set_position([0.01, 0.88, 0.05, 0.1])
    g.ax_heatmap.set_ylabel('cells')


    if False:
        # Summary
        # position 9226 is polymorphic within our population AND across phylogeny
        # and is around 50 bases after the end of GFP. It is first codon position
        # and makes for a nonsynonymous mutation Ala -> Thr at codon 292 (1-index)
        # - 556 including GFP - of the structural polyprotein which starts at 7561
        # in this reference. So it is inside E3, a protein part of envelope
        # other polymorphic positions: 400, 1612, 1615, 1618, 2594, 8301, 10551, 11691

        # pos 400 is usually a T, but we have mostly G here and a C in the reference.
        # It is also a synonymous mutation in a Leucin at third position, i.e. all
        # 4 nucleotides are fine.

        from seqanpy import align_overlap
        from Bio import SeqIO
        gene_hypo_aa = list(SeqIO.parse('../data/genbank/veev_genes_aa.fasta', 'fasta'))[0]

        # Start of first hypothetical CDS
        start = 44

        gene_struct_aa = list(SeqIO.parse('../data/genbank/veev_genes_aa.fasta', 'fasta'))[2]
        start = 7561

        # Split capsid, E3, E2, E1, fusion
        capsid_coo_genbank = [7562, 8386 + 1] # 275 residues
        e3_coo_genbank = [8387, 8563 + 1] # 59 residues
        capsid_coo = [x - capsid_coo_genbank[0] for x in capsid_coo_genbank]
        e3_coo = [x - capsid_coo_genbank[0] for x in e3_coo_genbank]


    print('Correlate gene expression with specific mutation frequencies')
    ind = covm > 1.5
    afs = data['allele_frequencies'][:, ind]
    dsh = ds.query_samples_by_name(afs['sample'].data)
    dsh.counts.normalize('counts_per_million', inplace=True)

    def plot_mutation_scatter(dsh, afs, mutation):

        refa = str(data['reference'][mutation[1]-1].values)
        label = mutation[0] + '_'+str(mutation[1])
        dsh.samplesheet[label] = afs.sel(position=mutation[1]-1, nucleotide=mutation[0]).values
        corr = dsh.correlation.correlate_features_phenotypes(label, method='spearman').fillna(0)

        # running average
        bins = np.linspace(0, 1, 21)
        wl = 0.1
        def sliding_window(x, y, bins, wl):
            bins = bins[bins + wl <= 1]
            xw = bins + 0.5 * wl
            yw = []
            for wli in bins:
                yw.append(y[(x >= wli) & (x <= wli + wl)].mean())
            return xw, yw

        fig = plt.figure(figsize=(13, 4))
        gs = fig.add_gridspec(2, 6)
        axs = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(4)]
        genes = corr.nlargest(4).index.tolist() + corr.nsmallest(4).index.tolist()
        x = dsh.samplesheet[label].values
        for gene, ax in zip(genes, axs):
            y = np.log10(0.1 + dsh.counts.loc[gene].values)
            ax.scatter(x, y, color='steelblue', alpha=0.2)

            xw, yw = sliding_window(x, y, bins, wl)
            ax.plot(xw, yw, lw=1.3, color='darkred')

            ax.set_ylabel('{:} [cpm]'.format(gene))
            ax.grid(True)
            ax.text(0.95, 0.95, '$\\rho = {:.2f}$'.format(corr.loc[gene]),
                    transform=ax.transAxes, ha='right', va='top')
            if ax == axs[-2]:
                ax.set_xlabel('$\\nu [ {:} \\rightarrow {:} ] \\quad at \\quad {:}$'.format(
                    refa, mutation[0], mutation[1]))
            ax.set_yticks([0, 2, 4, 6])
            if ax in [axs[0], axs[len(axs) // 2]]:
                ax.set_yticklabels(['$1$', '$10^2$', '$10^4$', '$10^6$'])
            else:
                ax.set_yticklabels(['' for x in ax.get_yticks()])

        from scipy.stats import spearmanr
        ax = fig.add_subplot(gs[:, 4:])
        y = dsh.samplesheet['log_virus_reads_per_million'].values
        rho = spearmanr(x, y)[0]
        ax.scatter(x, y, color='steelblue', alpha=0.2)
        ax.set_xlabel('$\\nu [ {:} \\rightarrow {:} ] \\quad at \\quad {:}$'.format(
            refa, mutation[0], mutation[1]))
        ax.set_ylabel('vRNA [cpm]')
        ax.grid(True)
        ax.text(0.95, 0.05, '$\\rho = {:.2f}$'.format(rho),
                transform=ax.transAxes, ha='right', va='bottom')
        ax.set_yticks([0, 2, 4, 6])
        ax.set_yticklabels(['$1$', '$10^2$', '$10^4$', '$10^6$'])

        # Plot running average
        xw, yw = sliding_window(x, y, bins, wl)
        ax.plot(xw, yw, lw=2, color='darkred')

        fig.tight_layout()

        return fig, corr

    mutation = ('A', 9227)
    fig, corr = plot_mutation_scatter(dsh, afs, mutation)

    plt.ion()
    plt.show()

    # Differential expression for mutations with small populations
    # TODO

    # Look at specific freak cells
    strange_cells = [
            '1001703009_B19',
            '1001703011_E15',
            '1001703011_G16',
            ]
    print('Coverage for strange cells')
    cov = data['coverage']
    covm = np.log10(np.maximum(9e-1, cov.mean(axis=1).data))
    fig, ax = plt.subplots(
        1, 1, figsize=(18, 4),
        )
    for cn in strange_cells:
        covi = cov.sel(sample=cn)
        x = np.arange(l)
        y = np.log10(np.maximum(9e-1, covi.data))
        ax.plot(x, y, lw=1, alpha=0.5, label=cn)#, color=colors[j])
    ax.grid(True)
    ax.set_xlabel('Position in VEEV-GFP [bases]')
    ax.set_ylabel('Coverage [# reads]')
    ax.set_yticks([0, 1, 2, 3, 4, 5])
    ax.set_yticklabels(['$0$', '$1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
    ax.legend()
    fig.tight_layout()

    print('Check coverage of the beginning VS end of the genome')
    def scatter_differential_viral_coverage(ind1, ind2, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(3.4, 3))
        else:
            fig = None

        snames = afs.sample.data

        # Normalize by total coverage in uniquely mapped reads
        norm = 1e-6 * ds.samplesheet.loc[snames, 'coverage'].values

        cov1 = cov.sel(sample=snames)[:, ind1[0]: ind1[1]]
        cov2 = cov.sel(sample=snames)[:, ind2[0]: ind2[1]]

        x = np.log10(0.1 + cov1.data).mean(axis=1) - np.log10(0.1 + norm)
        y = np.log10(0.1 + cov2.data).mean(axis=1) - np.log10(0.1 + norm)
        ax.scatter(x, y, color='darkred', alpha=0.2, zorder=3)
        xl = np.array([0, 5])
        ax.plot(xl, xl, lw=2, color='k', zorder=2, alpha=0.8)
        ax.plot(xl + 0.25, xl, lw=2, color='k', zorder=2, alpha=0.3)
        ax.plot(xl + 0.45, xl, lw=2, color='k', zorder=2, alpha=0.3)
        sns.kdeplot(x, y, ax=ax, alpha=0.5, zorder=4, n_levels=30)

        bins = np.linspace(0, 4, 21)
        wl = 0.2
        def sliding_window_y(x, y, bins, wl):
            bins = bins[bins + wl <= bins[-1]]
            yw = bins + 0.5 * wl
            xw = []
            for wli in bins:
                xw.append(x[(y >= wli) & (y <= wli + wl)].mean())
            xw = np.array(xw)
            ind = ~np.isnan(xw)
            xw = xw[ind]
            yw = yw[ind]
            return xw, yw

        # average in sliding window
        xw, yw = sliding_window_y(x, y, bins, wl)
        ind = yw >= 1.7
        yw = yw[ind]
        xw = xw[ind]

        from scipy.stats import gaussian_kde
        X, Y = 3. * np.mgrid[0:100, 0:100] / 100
        positions = np.vstack([X.ravel(), Y.ravel()])
        kernel = gaussian_kde([x, y])
        xmax, ymax = positions[:, kernel(positions).argmax()]

        # smoothen
        from scipy.interpolate import interp1d
        from scipy.signal import savgol_filter
        yy = np.linspace(yw.min(), yw.max(), 250)
        xx = interp1d(yw, xw, kind='linear')(yy)
        window_size, poly_order = 151, 3
        xx_sg = savgol_filter(xx, window_size, poly_order)

        # only plot smooth interpolation until the max of y, else artifacts
        ind = yy <= ymax
        yy = yy[ind]
        xx_sg = xx_sg[ind]
        ax.plot(xx_sg, yy, lw=2.5, color='lightgreen', alpha=0.7, zorder=5)
        ax.arrow(xx_sg[-2], yy[-2], xx_sg[-1] - xx_sg[-2], yy[-1] - yy[-2],
                 zorder=10, width=0.07, edgecolor='none', facecolor='lightgreen')

        # plot vertical bars
        ax.plot([2.3] * 2, [-1, 4], color='grey', alpha=0.4)
        ax.plot([2.9] * 2, [-1, 4], color='grey', alpha=0.4)

        ax.set_xlabel('Mean coverage bases {:}-{:}'.format(ind1[0]+1, ind1[-1]))
        ax.set_ylabel('Mean coverage bases {:}-{:}'.format(ind2[0]+1, ind2[-1]))
        ax.set_xticks([0, 1, 2, 3, 4, 5])
        ax.set_xticklabels(['$0$', '$1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
        ax.set_yticks([0, 1, 2, 3, 4, 5])
        ax.set_yticklabels(['$0$', '$1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
        ax.set_xlim(-0.3, 3.9)
        ax.set_ylim(-0.3, 3.9)
        ax.grid(True)
        if fig is not None:
            fig.tight_layout()

        return fig

    fig, axs = plt.subplots(1, 3, figsize=(9, 3))
    scatter_differential_viral_coverage([0, 1700], [7500, l], ax=axs[0])
    scatter_differential_viral_coverage([7500, 9000], [9000, l], ax=axs[1])
    scatter_differential_viral_coverage([0, 850], [850, 1700], ax=axs[2])
    fig.tight_layout()

    # Play some more
    if False:
        fig, axs = plt.subplots(1, 5, figsize=(15, 3))
        scatter_differential_viral_coverage([0, 1700], [7500, l], ax=axs[0])
        scatter_differential_viral_coverage([7500, 9000], [9000, l], ax=axs[1])
        scatter_differential_viral_coverage([0, 850], [850, 1700], ax=axs[2])
        scatter_differential_viral_coverage([0, 1700], [1700, 4000], ax=axs[3])
        scatter_differential_viral_coverage([1700, 2900], [2900, 4000], ax=axs[4])
        fig.tight_layout()

    print('Export data on 5p and 3p coverage')
    ind1, ind2 = [0, 1700], [7500, l]
    snames = afs.sample.data
    norm = 1e-6 * ds.samplesheet.loc[snames, 'coverage'].values
    cov1 = cov.sel(sample=snames)[:, ind1[0]: ind1[1]]
    cov2 = cov.sel(sample=snames)[:, ind2[0]: ind2[1]]
    x = np.log10(0.1 + cov1.data).mean(axis=1) - np.log10(0.1 + norm)
    y = np.log10(0.1 + cov2.data).mean(axis=1) - np.log10(0.1 + norm)
    df = pd.DataFrame(
        data=(10**np.vstack([x, y]).T) - 0.1,
        index=snames,
        columns=['vRNA_5p', 'vRNA_3p'],
        )
    df.index.name = 'cell'
    #df.to_csv('../data/dataset/vRNA_5p_3p.tsv', sep='\t', index=True)

    print('Correlate host gene expression with both')
    dsh.samplesheet['cov_5p'] = df['vRNA_5p']
    dsh.samplesheet['cov_3p'] = df['vRNA_3p']
    corr = dsh.correlation.correlate_features_phenotypes(
            phenotypes=['cov_5p', 'cov_3p']).fillna(0)

    fig, ax = plt.subplots(figsize=(10, 10))

    ind1 = (corr.max(axis=1) >= 0.2) | (corr.min(axis=1) <= -0.37)
    ax.scatter(corr.loc[ind1, 'cov_5p'], corr.loc[ind1, 'cov_3p'], alpha=0.25)
    ax.scatter(corr.loc[~ind1, 'cov_5p'], corr.loc[~ind1, 'cov_3p'], alpha=0.05)
    texts = []
    for i in ind1.values.nonzero()[0]:
        t = ax.text(
                corr['cov_5p'].iloc[i],
                corr['cov_3p'].iloc[i],
                corr.index[i],
                ha='center', va='center',
                fontsize=6,
                )
        texts.append(t)
    from adjustText import adjust_text
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='k'))

    ax.plot([-1, 1], [-1, 1], lw=1.5, color='k', alpha=0.5)

    #sns.kdeplot(
    #    corr.loc[:, 'cov_5p'],
    #    corr.loc[:, 'cov_3p'],
    #    ax=ax,
    #    n_levels=100)

    vmax = np.abs(corr.values).max() * 1.1
    ax.set_xlim(-vmax, vmax)
    ax.set_ylim(-vmax, vmax)
    ax.grid(True)
    ax.set_xlabel("correlation with 5' vRNA [bases 1-1700]")
    ax.set_ylabel("correlation with 3' vRNA [bases 7500-{:}]".format(l))

    fig.tight_layout()

    print('Same plot, compare with ZY')
    cellnames = pd.read_csv(
            '../data/Zhiyuan/list_cells_corrplots.tsv',
            sep='\t', index_col=0, squeeze=True, skiprows=1, header=None).values
    dshZY = dsh.query_samples_by_name(cellnames)
    corrZY = dshZY.correlation.correlate_features_phenotypes(
            phenotypes=['cov_5p', 'cov_3p']).fillna(0)

    fig, ax = plt.subplots(figsize=(10, 10))
    genes_write = {
            'orange': ['PFN2', 'TMED2', 'DPYSL2', 'CAPZA1', 'BROX', 'AXL', 'PUF60_1', 'EEF1A1P5'],
            'green': [
                'TSPAN6', 'SNAPC1', 'AHR', 'CHIC2', 'PRKAB1', 'ARRDC3', 'F3', 'DUSP1',
                'INTS12', 'BCL10', 'ARNT', 'PLK2', 'WTAP', 'FLCN', 'RCAN1', 'GEM',
                'ABCA1', 'TSC22D2',
                ],
            'purple': corr.loc[(corr['cov_3p'] > corr['cov_5p'] + 0.31)].index.tolist() + corr.loc[(corr['cov_5p'] > -0.05) & (corr['cov_3p'] > corr['cov_5p'] + 0.25)].index.tolist(),
            }
    genes_write['purple'] = [g for g in genes_write['purple'] if (g not in genes_write['orange']) and (g not in genes_write['green'])]
    ind1 = (corrZY.max(axis=1) >= 0.2) | (corrZY.min(axis=1) <= -0.37)
    ax.scatter(corrZY.loc[ind1, 'cov_5p'], corrZY.loc[ind1, 'cov_3p'], alpha=0.25)
    ax.scatter(corrZY.loc[~ind1, 'cov_5p'], corrZY.loc[~ind1, 'cov_3p'], alpha=0.05)
    texts = []
    for color, genes_writei in genes_write.items():
        for gene in genes_writei:
            tx, ty = corrZY['cov_5p'].loc[gene], corrZY['cov_3p'].loc[gene]
            t = ax.text(
                    tx, ty,
                    gene,
                    ha='center', va='center',
                    fontsize=10,
                    )
            ax.scatter([tx], [ty], s=100, color=color, alpha=0.7)
            texts.append(t)
    from adjustText import adjust_text
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='k'))

    ax.plot([-1, 1], [-1, 1], lw=1.5, color='k', alpha=0.5)

    #sns.kdeplot(
    #    corr.loc[:, 'cov_5p'],
    #    corr.loc[:, 'cov_3p'],
    #    ax=ax,
    #    n_levels=100)

    vmax = 0.65
    ax.set_xlim(-vmax, vmax)
    ax.set_ylim(-vmax, vmax)
    ax.set_xlabel("correlation with 5' vRNA [bases 1-1700]")
    ax.set_ylabel("correlation with 3' vRNA [bases 7500-{:}]".format(l))
    ax.set_xticks([-0.6, -0.3, 0, 0.3, 0.6])
    ax.set_xticks([-0.45, -0.15, 0.15, 0.45], minor=True)
    ax.set_yticks([-0.6, -0.3, 0, 0.3, 0.6])
    ax.set_yticks([-0.45, -0.15, 0.15, 0.45], minor=True)
    ax.grid(True, which='both')

    fig.tight_layout()

    plt.ion()
    plt.show()

