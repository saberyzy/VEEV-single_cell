# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/05/19
content:    Classic correlations, but do a better job at fitting curves
'''
import os
import sys
import numpy as np
import pandas as pd

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

    ds.query_samples_by_metadata('coverage >= 100000', inplace=True)

    ds.counts.normalize('counts_per_million', inplace=True)

    print('Filter mito/ribo/genes')
    feas = ds.featurenames
    good_genes = pd.Series(data=np.ones(ds.n_features, bool), index=ds.featurenames)
    good_genes[feas.str.startswith('MT-')] = False
    good_genes[feas.str.startswith('RPL')] = False
    good_genes[feas.str.startswith('RPS')] = False
    good_genes = good_genes.index[good_genes]
    ds.query_features_by_name(good_genes, inplace=True)

    print('Get average expression')
    ge_mean = ds.counts.get_statistics(metrics=('mean',)).iloc[:, 0]

    print('Exclude controls')
    dsnc = ds.query_samples_by_metadata('(MOI > 0) & (VEEV_counts >= 2)')

    print('Correlate')
    corr = dsnc.correlation.correlate_features_phenotypes('virus_reads_per_million').fillna(0)

    print('Choose a few good genes to start')
    genes = ['ARRDC3', 'SERF2', 'CXCL3', 'CXCL2', 'DDIT3', 'HILPDA', 'NIPA2',
             'DEGS1', 'SBNO1', 'BLCAP', 'EIF1AD']

    fig, axs = plt.subplots(2, 5, figsize=(15, 5.5), sharex=True, sharey=True)
    axs = axs.ravel()
    x = np.log10(0.1 + dsnc.samplesheet['virus_reads_per_million'].values)
    for gene, ax in zip(genes, axs):
        y = np.log10(0.1 + dsnc.counts.loc[gene].values)
        ax.scatter(x, y, alpha=0.1)
        ax.grid(True)
        ax.set_xticks([-1, 0, 1, 2, 3, 4, 5])
        ax.set_xticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
        ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
        ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
        ax.text(5.6, 4.8,
                gene,
                ha='right', va='top',
                bbox=dict(facecolor=(1, 1, 1, 0.6), edgecolor='none', pad=2.0))

        # Get averages by window
        if False:
            bins = np.linspace(-0.5, 6, 50)
            overlap = 5
            left = bins[:-overlap]
            right = bins[overlap:]
            ym = []
            for l, r in zip(left, right):
                i = (x >= l) & (x < r)
                ym.append(y[i].mean())
            xm = 0.5 * (left + right)
            ax.plot(xm, ym, lw=2, color='darkred', alpha=0.5)

        # Interpolation
        from scipy.interpolate import interp1d
        from scipy.signal import savgol_filter
        xx = np.linspace(x.min(), x.max(), 100)
        itp = interp1d(x, y, kind='linear')
        window_size, poly_order = 51, 2
        yy_sg = savgol_filter(itp(xx), window_size, poly_order)
        ax.plot(xx, yy_sg, lw=2, color='seagreen', alpha=0.8, zorder=9)

        # Fill quartiles
        def plot_percentile_shadow(low=25, high=75, color='seagreen', alpha=0.2, zorder=2):
            bins = np.linspace(-0.5, 6, 50)
            overlap = 5
            left = bins[:-overlap]
            right = bins[overlap:]
            y25 = np.empty(len(left))
            y75 = np.empty(len(left))
            for ii, (l, r) in enumerate(zip(left, right)):
                i = (x >= l) & (x < r)
                y25[ii], y75[ii] = np.percentile(y[i], [low, high])
            xm = 0.5 * (left + right)
            # Interpolate quartiles
            xxq = np.linspace(xm.min(), xm.max(), 100)
            itp25 = interp1d(xm, y25, kind='linear')
            itp75 = interp1d(xm, y75, kind='linear')
            yy_sg25 = savgol_filter(itp25(xxq), window_size, poly_order)
            yy_sg75 = savgol_filter(itp75(xxq), window_size, poly_order)
            ax.fill_between(
                    xxq, yy_sg25, yy_sg75,
                    color=color, alpha=alpha, zorder=zorder)
        plot_percentile_shadow(low=25, high=75, color='seagreen', alpha=0.2, zorder=2)
        plot_percentile_shadow(low=10, high=90, color='seagreen', alpha=0.07, zorder=1)

        # Fit threshold-linear
        from scipy.optimize import curve_fit
        def fun(x, b, i, s):
            y = i + s * x
            t = (b - i) / s
            y[x <= t] = b
            return y
        kwargs = {}
        (b, i, s), pcov = curve_fit(
            fun,
            xdata=xx,
            ydata=yy_sg,
            **kwargs)
        t = (b - i) / s
        yy_fit = fun(xx, b, i, s)
        ax.plot(xx, yy_fit, lw=2, alpha=0.7, color='tomato')

        if t >= 0:
            ax.plot([t] * 2, [-1.1, 6], lw=1.5, color='dimgray', alpha=0.4)
        ax.text(0.15, 4.8,
                '$t = {:.1f}, \\Delta = {:.1f}$'.format(t, yy_fit.max() - yy_fit.min()),
                ha='left', va='top',
                bbox=dict(facecolor='white', edgecolor='dimgrey', pad=2.0))

        ax.set_xlim(0, 5.8)
        ax.set_ylim(-1.1, 5)

    fig.text(0.52, 0.01, 'vRNA [cpm]', ha='center', va='bottom')
    fig.text(0.01, 0.52, 'gene expression [cpm]', ha='left', va='center', rotation=90)
    fig.tight_layout(rect=(0.015, 0.015, 1, 1))


    if False:
        print('Cross-check the list with controls')
        dsc = ds.query_samples_by_metadata('MOI == 0')
        corr0 = dsc.correlation.correlate_features_phenotypes('time [h]').fillna(0)
        corr0.name = 'control with time'

        df = pd.concat([corr, corr0], axis=1)
        df.rename(columns={'correlation': 'vRNA'}, inplace=True)

