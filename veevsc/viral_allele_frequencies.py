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
    ind = covm > 2
    afs = data['allele_frequencies'][:, ind]
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
    afs_nonref_poly = np.zeros((afs_poly.shape[1], afs_poly.shape[2]))
    for j, refij in enumerate(refi_poly):
        afs_nonref_poly[:, j] = afs_poly[refij, :, j]
    afs_nonref_poly = 1 - afs_nonref_poly
    afs_nonref_poly = pd.DataFrame(
            afs_nonref_poly,
            index=data.sample.data[ind],
            columns=pos_poly,
            )

    # Cell colors by total amount of reads
    dsh = ds.query_samples_by_name(afs_nonref_poly.index)
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

    # Order
    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage
    from polo import optimal_leaf_ordering
    D = pdist(afs_nonref_poly.values, 'euclidean')
    Z = linkage(D, 'average')
    optimal_Z = optimal_leaf_ordering(Z, D)

    #g = sns.clustermap(maj_snv_plot, figsize=(21, 8), xticklabels=True, yticklabels=True); plt.subplots_adjust(0.05, 0.05, 0.95, 0.95)
    g = sns.clustermap(
        afs_nonref_poly,
        figsize=(21, 8),
        xticklabels=True,
        yticklabels=True,
        row_colors=colors,
        row_linkage=optimal_Z,
        )
    plt.subplots_adjust(0.03, 0.08, 0.92, 0.95)

    # replot without clustering positions
    g = sns.clustermap(
        afs_nonref_poly,
        figsize=(21, 8),
        xticklabels=True,
        yticklabels=True,
        row_colors=colors,
        row_linkage=optimal_Z,
        col_cluster=False,
        )
    plt.subplots_adjust(0.03, 0.08, 0.92, 0.95)


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

    # Differential expression for mutations with small populations


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
        cov1 = cov.sel(sample=snames)[:, ind1[0]: ind1[1]]
        cov2 = cov.sel(sample=snames)[:, ind2[0]: ind2[1]]
        x = np.log10(0.1 + cov1.data).mean(axis=1)
        y = np.log10(0.1 + cov2.data).mean(axis=1)
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
        ind = yw >= 1.2
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
        ax.set_xlim(-0.3, 3.5)
        ax.set_ylim(-0.3, 3.5)
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
    fig, axs = plt.subplots(1, 5, figsize=(15, 3))
    scatter_differential_viral_coverage([0, 1700], [7500, l], ax=axs[0])
    scatter_differential_viral_coverage([7500, 9000], [9000, l], ax=axs[1])
    scatter_differential_viral_coverage([0, 850], [850, 1700], ax=axs[2])
    scatter_differential_viral_coverage([0, 1700], [1700, 4000], ax=axs[3])
    scatter_differential_viral_coverage([1700, 2900], [2900, 4000], ax=axs[4])
    fig.tight_layout()

    plt.ion()
    plt.show()


