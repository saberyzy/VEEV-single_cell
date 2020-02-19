# vim: fdm=indent
'''
author:     Fabio Zanini
date:       24/11/19
content:    Correlate viral SNPs with host gene expression
'''
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, leaves_list

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

    # FIXME: somehow MOI for column 24 in the plate is set to 0, reset it to 1
    dsh.samplesheet['plate'] = [x[8:10] for x in dsh.samplesheet.index]
    dsh.samplesheet['well'] = [x[11:] for x in dsh.samplesheet.index]
    dsh.samplesheet['well_col'] = [x[12:] for x in dsh.samplesheet.index]
    dsh.samplesheet.loc[dsh.samplesheet['well_col'] == '24', 'MOI'] = 1

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

    print('Prepare clustermaps')
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

    print('Plot allele frequencies')
    if False:
        # Cell colors
        v = dsh.samplesheet['log_virus_reads_per_million'].fillna(-1).values
        v = (v - v.min()) / (v.max() - v.min())
        cmap = sns.dark_palette('Blue', as_cmap=True)
        colors_viral = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)

        v = dsh.samplesheet['time_index'].values
        v = (v - v.min()) / (v.max() - v.min())
        cmap = sns.dark_palette('Green', as_cmap=True)
        colors_time = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)

        v = dsh.samplesheet['MOI'].values
        v[v > 0] += 1
        v = (v - v.min()) / (v.max() - v.min())
        cmap = sns.dark_palette('Orange', as_cmap=True)
        colors_moi = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)

        colors = pd.concat([colors_viral, colors_time, colors_moi], axis=1)
        colors.columns = ['vRNA', 'time', 'MOI']

        ind_poss = afs_nonref_poly.std(axis=0).nlargest(10).sort_index().index
        datap = afs_nonref_poly[ind_poss]
        datap.columns = pd.Index(
            [x + 1 for x in datap.columns],
            name='Position in VEEV-GFP genome [bases, 1-indexed]',
            )

        D = pdist(datap.values, 'euclidean')
        Z = linkage(D, 'average', optimal_ordering=True)

        # Sort cells by vRNA
        ind_sort = dsh.samplesheet['log_virus_reads_per_million'].fillna(-1).values.argsort()

        g = sns.clustermap(
            datap.iloc[ind_sort],
            figsize=(4, 8),
            xticklabels=True,
            yticklabels=True,
            row_colors=colors.iloc[ind_sort],
            row_linkage=Z,
            col_cluster=False,
            row_cluster=False,
            )
        plt.subplots_adjust(0.01, 0.08, 0.92, 0.95)
        ax = g.ax_col_dendrogram
        ax.set_xlim(0, datap.shape[1])
        i = datap.columns.tolist().index(9227)
        ax.arrow(i + 0.5, 0.3, 0, -0.25, ec='red', fc='red', head_width=0.5, head_length=0.15, length_includes_head=True)
        ax.text(i + 0.5, 0.35, '9227 A>G, Ala>Thr, E3 protein', ha='center', va='bottom', fontsize=10)
        ax1 = g.fig.get_axes()[-1]
        ax1.set_ylabel('Frequency of\nalternate allele', rotation=0, labelpad=15, va='center', ha='left')
        g.ax_heatmap.set_yticklabels(['' for i in range(datap.shape[0])])
        g.ax_heatmap.set_position([0.07, 0.10, 0.91, 0.76])
        g.ax_col_dendrogram.set_position([0.07, 0.86, 0.91, 0.1])
        g.ax_row_colors.set_position([0.01, 0.10, 0.05, 0.76])
        g.fig.get_axes()[-1].set_position([0.01, 0.88, 0.05, 0.1])
        g.ax_heatmap.set_ylabel('cells')

    if False:
        print('Plot coverage')
        # Cell colors
        v = dsh.samplesheet['log_virus_reads_per_million'].fillna(-1).values
        v = (v - v.min()) / (v.max() - v.min())
        cmap = sns.dark_palette('Blue', as_cmap=True)
        colors_viral = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)

        v = dsh.samplesheet['time_index'].values
        v = (v - v.min()) / (v.max() - v.min())
        cmap = sns.dark_palette('Green', as_cmap=True)
        colors_time = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)

        v = dsh.samplesheet['MOI'].values
        v[v > 0] += 1
        v = (v - v.min()) / (v.max() - v.min())
        cmap = sns.dark_palette('Orange', as_cmap=True)
        colors_moi = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)

        colors = pd.concat([colors_viral, colors_time, colors_moi], axis=1)
        colors.columns = ['vRNA', 'time', 'MOI']

        ind_poss = afs_nonref_poly.std(axis=0).nlargest(10).sort_index().index
        datap = np.log10(0.1 + covh_poly.to_pandas()[ind_poss])
        datap.columns = pd.Index(
            [x + 1 for x in datap.columns],
            name='Position in VEEV-GFP genome [bases, 1-indexed]',
            )

        # Sort cells by vRNA
        ind_sort = dsh.samplesheet['log_virus_reads_per_million'].fillna(-1).values.argsort()

        g = sns.clustermap(
            datap.iloc[ind_sort],
            figsize=(4, 8),
            xticklabels=True,
            yticklabels=True,
            row_colors=colors.iloc[ind_sort],
            row_linkage=Z,
            col_cluster=False,
            row_cluster=False,
            )
        plt.subplots_adjust(0.01, 0.08, 0.92, 0.95)
        ax = g.ax_col_dendrogram
        ax.set_xlim(0, datap.shape[1])
        i = datap.columns.tolist().index(9227)
        ax.arrow(i + 0.5, 0.3, 0, -0.25, ec='red', fc='red', head_width=0.5, head_length=0.15, length_includes_head=True)
        ax.text(i + 0.5, 0.35, '9227 A>G, Ala>Thr, E3 protein', ha='center', va='bottom', fontsize=10)
        ax1 = g.fig.get_axes()[-1]
        ax1.set_ylabel('Coverage', rotation=0, labelpad=15, va='center', ha='left')
        g.ax_heatmap.set_yticklabels(['' for i in range(datap.shape[0])])
        g.ax_heatmap.set_position([0.07, 0.10, 0.91, 0.76])
        g.ax_col_dendrogram.set_position([0.07, 0.86, 0.91, 0.1])
        g.ax_row_colors.set_position([0.01, 0.10, 0.05, 0.76])
        g.ax_heatmap.set_ylabel('cells')
        cax = g.fig.get_axes()[-1]
        cax.set_position([0.01, 0.88, 0.05, 0.1])

    if True:
        print('Plot allele freqs and coverage')
        fig = plt.figure(figsize=(7, 6))
        gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 4, 4], height_ratios=[1, 3])
        axs = []
        axs.append(fig.add_subplot(gs[0, 3]))
        axs.append(fig.add_subplot(gs[0, 4]))
        axs.append(fig.add_subplot(gs[1, 0]))
        axs.append(fig.add_subplot(gs[1, 1]))
        axs.append(fig.add_subplot(gs[1, 2]))
        axs.append(fig.add_subplot(gs[1, 3]))
        axs.append(fig.add_subplot(gs[1, 4]))

        # Plot data
        ind_poss = afs_nonref_poly.std(axis=0).nlargest(10).sort_index().index
        ind_poss = [i for i in ind_poss if np.log10(0.1 + covh_poly.sel(position=i)).mean() >= 2]
        datap = afs_nonref_poly[ind_poss]
        datap.columns = pd.Index(
            [x + 1 for x in datap.columns],
            name='Position in VEEV-GFP genome [bases, 1-indexed]',
            )
        datac = np.log10(0.1 + covh_poly.to_pandas()[ind_poss])
        datac.columns = pd.Index(
            [x + 1 for x in datac.columns],
            name='Position in VEEV-GFP genome [bases, 1-indexed]',
            )
        # Sort cells by vRNA
        ind_sort = dsh.samplesheet['log_virus_reads_per_million'].fillna(-1).values.argsort()

        # Cell colors
        v = dsh.samplesheet['log_virus_reads_per_million'].fillna(-1).values
        vabs = v.copy()
        v = (v - v.min()) / (v.max() - v.min())
        cmap = sns.dark_palette('Blue', as_cmap=True)
        colors_viral = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)
        ax = axs[2]
        ax.imshow(np.array([v]).T[ind_sort], cmap=cmap, vmin=0, vmax=1)
        ax.set_aspect(0.05)
        ax.set_xticks([0])
        ax.set_xticklabels(['vRNA'], rotation=90)
        yticks = [0, 0]
        yticklabels = ['$10^4$', '$10^5$']
        ytickvalues = [4, 5]
        vabs_sorted = vabs[ind_sort]
        for ii, vi in enumerate(vabs_sorted[:-1]):
            for it, tick in enumerate(ytickvalues):
                if (vi < tick) and (vabs_sorted[ii + 1] >= tick):
                    yticks[it] = ii
                    break
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.set_ylim(datap.shape[0] - 0.5, -0.5)

        v = dsh.samplesheet['time_index'].values
        vabs = v.copy()
        v = (v - v.min()) / (v.max() - v.min())
        cmap = sns.dark_palette('Green', as_cmap=True)
        colors_time = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)
        ax = axs[3]
        ax.imshow(np.array([v]).T[ind_sort], cmap=cmap, vmin=0, vmax=1)
        ax.set_aspect(0.05)
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_xticks([0])
        ax.set_xticklabels(['Time'], rotation=90)
        vu = np.unique(v)
        vau = np.unique(vabs)
        cu = cmap(vu)
        handles, labels = [], []
        d = {}
        for ti, tim in dsh.samplesheet[['time_index', 'time [h]']].values:
            d[ti] = tim
        from matplotlib.patches import Rectangle
        for vi, vai, ci in zip(vu, vau, cu):
            label = '{:.0f}'.format(d[vai])
            h = Rectangle((0, 0), 0, 0, color=ci)
            handles.append(h)
            labels.append(label)
        leg = ax.legend(
            handles, labels,
            title='Time [hpi]\n',
            loc='upper left',
            bbox_to_anchor=(-2, 1.49),
            bbox_transform=ax.transAxes,
            )
        leg.get_frame().set_linewidth(0.0)

        v = dsh.samplesheet['MOI'].values.copy()
        vabs = v.copy()
        v[v > 0] += 1
        v = (v - v.min()) / (v.max() - v.min())
        cmap = sns.dark_palette('Red', as_cmap=True)
        colors_moi = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)
        ax = axs[4]
        ax.imshow(np.array([v]).T[ind_sort], cmap=cmap, vmin=0, vmax=1)
        ax.set_aspect(0.05)
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_xticks([0])
        ax.set_xticklabels(['MOI'], rotation=90)
        vu = np.unique(v)
        vau = np.unique(vabs)
        cu = cmap(vu)
        handles, labels = [], []
        d = {0: '0', 0.1: '0.1', 1: '1'}
        from matplotlib.patches import Rectangle
        for vi, vai, ci in zip(vu, vau, cu):
            label = d[np.round(vai, 1)]
            h = Rectangle((0, 0), 0, 0, color=ci)
            handles.append(h)
            labels.append(label)
        leg = ax.legend(
            handles, labels,
            title='MOI\n',
            loc='upper left',
            bbox_to_anchor=(-0.75, 1.49),
            bbox_transform=ax.transAxes,
            )
        leg.get_frame().set_linewidth(0.0)

        colors = pd.concat([colors_viral, colors_time, colors_moi], axis=1)
        colors.columns = ['vRNA', 'time', 'MOI']

        # Allele freqs
        from matplotlib.ticker import Formatter
        class Fmter(Formatter):
            def format_ticks(self, values):
                dic = dict(zip(
                    [0, 0.25, 0.5, 0.75, 1],
                    ['0%', '25%', '50%', '75%', '100%']),
                    )
                return [dic[v] for v in values]
        ax = axs[5]
        sns.heatmap(
            data=datap.iloc[ind_sort],
            ax=ax,
            cmap='plasma',
            xticklabels=True,
            cbar=True,
            cbar_ax=axs[0],
            vmin=0,
            vmax=1,
            cbar_kws={
                'ticks': [0, 0.25, 0.5, 0.75, 1],
                'format': Fmter(),
                'orientation': 'horizontal',
                }
            )
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_ylabel('')
        ax.set_xlabel('')
        axs[0].set_aspect(1./6)
        axs[0].set_xlabel('')
        axs[0].xaxis.set_ticks_position('top')
        axs[0].set_title('Frequency of\nmutant allele', pad=10, fontsize=10)
        axs[0].set_xlim(-0.2, 1.2)

        # Coverage
        from matplotlib.ticker import Formatter
        class Fmter(Formatter):
            def format_ticks(self, values):
                dic = dict(zip(
                    [-1, 0, 1, 2, 3, 4, 5],
                    ['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$']),
                    )
                return [dic[v] for v in values]
        ax = axs[6]
        sns.heatmap(
            data=datac.iloc[ind_sort],
            ax=ax,
            cmap='plasma',
            xticklabels=True,
            cbar=True,
            cbar_ax=axs[1],
            cbar_kws={
                'ticks': [-1, 0, 1, 2, 3, 4, 5],
                'format': Fmter(),
                'orientation': 'horizontal',
                }
            )
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_ylabel('')
        ax.set_xlabel('Position in VEEV-GFP genome [bases]                                             ', labelpad=10)
        axs[1].set_aspect(1./6)
        axs[1].set_xlabel('')
        axs[1].xaxis.set_ticks_position('top')
        axs[1].set_title('Coverage\n[n reads]', pad=10, fontsize=10)
        axs[1].set_xlim(-1.5, 4.5)

        #fig.savefig('../figures/snp_and_coverage_VEEV.png', dpi=600)
        #fig.savefig('../figures/snp_and_coverage_VEEV.pdf')
        #fig.savefig('../figures/snp_and_coverage_VEEV.svg')

    if False:
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

            fig, axs = plt.subplots(1, 4, figsize=(6.2, 1.8))
            genes = corr.nsmallest(3).index.tolist()
            x = dsh.samplesheet[label].values
            for gene, ax in zip(genes, axs):
                y = np.log10(0.1 + dsh.counts.loc[gene].values)
                ax.scatter(x, y, color='steelblue', alpha=0.2)

                xw, yw = sliding_window(x, y, bins, wl)
                ax.plot(xw, yw, lw=1.3, color='darkred')

                ax.set_title(gene)
                ax.grid(True)
                ax.text(0.95, 0.95, '$\\rho = {:.2f}$'.format(corr.loc[gene]),
                        transform=ax.transAxes, ha='right', va='top')
                ax.set_xlabel('$\\nu [ {:} \\rightarrow {:} ] \\quad at \\quad {:}$'.format(
                    refa, mutation[0], mutation[1]))
                ax.set_yticks([0, 2, 4, 6])
                if ax == axs[0]:
                    ax.set_yticklabels(['$1$', '$10^2$', '$10^4$', '$10^6$'])
                    ax.set_ylabel('cpm')
                else:
                    ax.set_yticklabels(['' for x in ax.get_yticks()])
                    ax.set_ylabel('')

            from scipy.stats import spearmanr
            ax = axs[-1]
            y = dsh.samplesheet['log_virus_reads_per_million'].values
            rho = spearmanr(x, y)[0]
            ax.scatter(x, y, color='steelblue', alpha=0.2)
            ax.set_xlabel('$\\nu [ {:} \\rightarrow {:} ] \\quad at \\quad {:}$'.format(
                refa, mutation[0], mutation[1]))
            ax.set_title('vRNA')
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
        #fig.savefig('../figures/snp_and_host_gene_expression_VEEV.png', dpi=600)
        #fig.savefig('../figures/snp_and_host_gene_expression_VEEV.pdf')
        #fig.savefig('../figures/snp_and_host_gene_expression_VEEV.svg')

    if True:
        print('Plot allele freqs and coverage, smaller')
        fig = plt.figure(figsize=(7, 4))
        gs = fig.add_gridspec(4, 3, height_ratios=[3, 4.2, 1, 1], width_ratios=[1, 1.8, 1.5])
        axs = []
        axs.append(fig.add_subplot(gs[0, 1]))
        axs.append(fig.add_subplot(gs[0, 2]))
        axs.append(fig.add_subplot(gs[1, :]))
        axs.append(fig.add_subplot(gs[2, :]))
        axs.append(fig.add_subplot(gs[3, :]))

        # Plot data
        ind_poss = afs_nonref_poly.std(axis=0).nlargest(1).sort_index().index
        ind_poss = [i for i in ind_poss if np.log10(0.1 + covh_poly.sel(position=i)).mean() >= 2]
        datap = afs_nonref_poly[ind_poss]
        datap.columns = pd.Index(
            [x + 1 for x in datap.columns],
            name='Position in VEEV-GFP genome [bases, 1-indexed]',
            )
        datac = np.log10(0.1 + covh_poly.to_pandas()[ind_poss])
        datac.columns = pd.Index(
            [x + 1 for x in datac.columns],
            name='Position in VEEV-GFP genome [bases, 1-indexed]',
            )
        # Sort cells by vRNA
        ind_sort = dsh.samplesheet['log_virus_reads_per_million'].fillna(-1).values.argsort()

        # Cell colors
        v = dsh.samplesheet['log_virus_reads_per_million'].fillna(-1).values
        vabs = v.copy()
        v = (v - v.min()) / (v.max() - v.min())
        cmap = sns.dark_palette('Blue', as_cmap=True)
        colors_viral = pd.Series(data=[tuple(x) for x in cmap(v)], index=dsh.samplenames)
        ax = axs[2]
        ax.imshow(np.array([v]).T[ind_sort].T, cmap=cmap, vmin=0, vmax=1)
        ax.set_aspect(15)
        ax.set_yticks([0])
        ax.set_yticklabels(['vRNA\n[cpm]'], rotation=0)
        xticklabels = ['$10^4$', '$10^5$', '$2 x 10^5$', '$3 x 10^5$']
        xtickvalues = [4, 5, np.log10(2e5), np.log10(3e5)]
        xticks = [0 for x in xticklabels]
        vabs_sorted = vabs[ind_sort]
        for ii, vi in enumerate(vabs_sorted[:-1]):
            for it, tick in enumerate(xtickvalues):
                if (vi < tick) and (vabs_sorted[ii + 1] >= tick):
                    xticks[it] = ii
                    break
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_xlim(-0.5, datap.shape[0] - 0.5)
        ax.xaxis.set_ticks_position('top')

        # Allele freqs
        from matplotlib.ticker import Formatter

        class Fmter(Formatter):
            def format_ticks(self, values):
                dic = dict(zip(
                    [0, 0.25, 0.5, 0.75, 1],
                    ['0%', '25%', '50%', '75%', '100%']),
                    )
                return [dic[v] for v in values]
        ax = axs[3]
        sns.heatmap(
            data=datap.iloc[ind_sort].T,
            ax=ax,
            cmap='plasma',
            xticklabels=True,
            cbar=True,
            cbar_ax=axs[0],
            vmin=0,
            vmax=1,
            cbar_kws={
                'ticks': [0, 0.25, 0.5, 0.75, 1],
                'format': Fmter(),
                'orientation': 'horizontal',
                }
            )
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_ylabel('')
        ax.set_yticklabels(['Allele\nfrequency'], rotation=0, ha='right')
        axs[0].set_aspect(1./6)
        axs[0].set_xlabel('')
        axs[0].xaxis.set_ticks_position('top')
        for tk in axs[0].get_xticklabels():
            tk.set_rotation(-30)
        axs[0].set_title('Frequency of\nmutant allele', pad=10, fontsize=10)
        axs[0].set_xlim(-0.2, 1.2)

        # Coverage
        from matplotlib.ticker import Formatter
        class Fmter(Formatter):
            def format_ticks(self, values):
                dic = dict(zip(
                    [-1, 0, 1, 2, 3, 4, 5],
                    ['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$']),
                    )
                return [dic[v] for v in values]
        ax = axs[4]
        sns.heatmap(
            data=datac.iloc[ind_sort].T,
            ax=ax,
            vmin=0, vmax=3.2,
            cmap='plasma',
            xticklabels=True,
            cbar=True,
            cbar_ax=axs[1],
            cbar_kws={
                'ticks': [0, 1, 2, 3],
                'format': Fmter(),
                'orientation': 'horizontal',
                }
            )
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_yticklabels(['Coverage'], rotation=0, ha='right')
        #ax.set_ylabel('Position in VEEV-GFP genome [bases]                                             ', labelpad=10)
        axs[1].set_aspect(1./6)
        axs[1].set_xlabel('')
        axs[1].xaxis.set_ticks_position('top')
        axs[1].set_title('Coverage\n[n reads]', pad=10, fontsize=10)
        axs[1].set_xlim(-1.5, 4.5)

        fig.savefig('../figures/snp_and_coverage_VEEV_only9227_raw.svg')

    plt.ion()
    plt.show()


