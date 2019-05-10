# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/19
content:    Explore the correlation between vRNA abundance and gene expression
            like in the eLife paper. Zhiyuan finds a negative median?
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

    print('Loading data...')
    counts = singlet.CountsTable(pd.read_csv(
            '../data/dataset/VEEVcounts.tsv',
            sep='\t',
            index_col=0,
            ).astype(np.float32),
            )
    print('Data loaded.')


    ds = singlet.Dataset(
            counts_table=counts,
            )

    ds.counts._spikeins = list(ds.counts.index[-103:-6])
    ds.counts._otherfeatures = list(ds.counts.index[-6:])
    ds.samplesheet['coverage'] = ds.counts.iloc[:-103].sum(axis=0)
    ds.samplesheet['VEEV_counts'] = ds.counts.loc['VEEVGFP']
    ds.samplesheet['virus_reads_per_million'] = 1e6 * ds.samplesheet['VEEV_counts'] / (ds.samplesheet['VEEV_counts'] + ds.counts.iloc[:-103].sum(axis=0))
    ds.samplesheet['log_virus_reads_per_million'] = np.log10(0.1 + ds.samplesheet['virus_reads_per_million'])

    # Metadata
    ds.samplesheet['plate'] = np.array([x[:10] for x in ds.samplesheet.index])
    ds.samplesheet['well'] = np.array([x[11:] for x in ds.samplesheet.index])
    ds.samplesheet['plate_col'] = np.array([x[12:] for x in ds.samplesheet.index], int)
    ds.samplesheet['time_index'] = (np.array([x[8:10] for x in ds.samplesheet.index], int) - 1) // 2
    ds.samplesheet['MOI'] = 0
    ds.samplesheet.loc[(ds.samplesheet['time_index'] == 0) & (ds.samplesheet['plate_col'] >= 6) & (ds.samplesheet['plate_col'] < 12), 'MOI'] = 0.1
    ds.samplesheet.loc[(ds.samplesheet['time_index'] == 0) & (ds.samplesheet['plate_col'] >= 12) & (ds.samplesheet['plate_col'] < 18), 'MOI'] = 1
    ds.samplesheet.loc[(ds.samplesheet['time_index'] == 0) & (ds.samplesheet['plate_col'] >= 18) & (ds.samplesheet['plate_col'] < 24), 'MOI'] = 3
    ds.samplesheet.loc[(ds.samplesheet['time_index'] != 0) & (ds.samplesheet['plate_col'] >= 8) & (ds.samplesheet['plate_col'] < 16), 'MOI'] = 0.1
    ds.samplesheet.loc[(ds.samplesheet['time_index'] != 0) & (ds.samplesheet['plate_col'] >= 16) & (ds.samplesheet['plate_col'] < 24), 'MOI'] = 1
    times = np.array([0.5, 1.5, 4, 6, 12, 24])
    ds.samplesheet['time [h]'] = times[ds.samplesheet['time_index']]

    sys.exit()


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

    print('Exclude controls')
    dsnc = ds.query_samples_by_metadata('(MOI > 0) & (VEEV_counts >= 2)')

    print('Correlate')
    corr = dsnc.correlation.correlate_features_phenotypes('virus_reads_per_million').fillna(0)

    print('Plot histogram and extremities')
    fig, axs = plt.subplots(1, 3, figsize=(10, 4))
    axs = axs.ravel()

    ax = axs[0]
    bins = [-1, -0.9, -0.8, -0.7, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3,
            -0.25, -0.2, -0.15, -0.09, -0.03]
    bins = np.array(bins + [-b for b in bins[::-1]])
    binsc = 0.5 * (bins[1:] + bins[:-1])
    h = np.histogram(corr.values, bins=bins, density=True)[0]
    ymin = 0.1 * h[h>0].min()
    ax.bar(
            bins[:-1], h + ymin, bins[1:] - bins[:-1],
            color='steelblue',
            )
    ax.set_xlabel('Spearman rho')
    ax.set_ylabel('Probability density')
    ax.set_ylim(ymin, 1.2 * h.max())
    ax.set_yscale('log')
    ax.grid(True)

    gtop = corr.idxmax()
    gbot = corr.idxmin()
    x = np.log10(0.1 + dsnc.samplesheet['virus_reads_per_million'].values)
    for (ax, gene) in zip(axs[1:], [gtop, gbot]):
        y = np.log10(0.1 + dsnc.counts.loc[gene].values)
        ax.scatter(x, y, alpha=0.2, color='steelblue')
        ax.set_xlim(-1.1, 6)
        ax.set_ylim(-1.1, 6)
        ax.grid(True)
        ax.set_xlabel('vRNA [cpm]')
        ax.set_ylabel('{:} expression [cpm]'.format(gene))

    fig.tight_layout()


    print('Plot a few more examples')
    gtops = corr.nlargest(30).index.tolist()
    gbots = corr.nsmallest(30).index.tolist()
    gplot = gtops + gbots[::-1]

    # Get bootstrapped values
    dsi = dsnc.query_features_by_name(gplot)
    corrb = []
    for ib in range(100):
        dsib = dsi.bootstrap()
        corri = dsib.correlation.correlate_features_phenotypes('virus_reads_per_million').fillna(0)
        corri.name = 'bs_{:}'.format(ib+1)
        corrb.append(corri)
    corrb = pd.concat(corrb, axis=1)
    corravg = corrb.mean(axis=1)
    corrstd = corrb.std(axis=1)

    # Threshold-linear fits
    dsnc.counts.log(inplace=True)
    fits = dsnc.fit.fit_single(xs=['log_virus_reads_per_million'], ys=gtops + gbots[::-1], model='threshold-linear')
    dsnc.counts.unlog(inplace=True)
    def fun(x, b, i, s):
        y = i + s * x
        t = (b - i) / s
        y[x <= t] = b
        return y

    fig, axs = plt.subplots(6, 10, figsize=(19.2, 9.2), sharex=True, sharey=True)
    axs = axs.ravel()
    x = np.log10(0.1 + dsnc.samplesheet['virus_reads_per_million'].values)
    for (ax, gene) in zip(axs, gplot):
        y = np.log10(0.1 + dsnc.counts.loc[gene].values)
        ax.scatter(x, y, alpha=0.2, color=('steelblue' if gene in gtops else 'darkred'))

        ax.set_xlim(-1.1, 6)
        ax.set_ylim(-1.1, 6)
        ax.grid(True)

        b, i, s, _ = fits.sel(y=gene).data[0]
        t = (b - i) / s
        xfit = np.linspace(-1, 6)
        yfit = fun(xfit, b, i, s)
        ax.plot(xfit, yfit, lw=1.5, color='k', alpha=0.7, ls='--')

        ax.text(0.02, 0.98, '{:}={:.2f}({:.0f})'.format(gene, corr.at[gene], corrstd.at[gene] * 100),
                ha='left', va='top',
                transform=ax.transAxes,
                bbox=dict(facecolor='none', edgecolor='black', pad=1.1))


    fig.text(0.52, 0.01, 'vRNA [cpm]')
    fig.text(0.01, 0.52, 'gene expression [cpm]', rotation=90)

    fig.tight_layout(rect=(0.01, 0.01, 1, 1), w_pad=0.1, h_pad=0.1)


    print('Plot a few Zhiyuan likes')
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))
    gene = 'BMP2K'
    x = np.log10(0.1 + dsnc.samplesheet['virus_reads_per_million'].values)
    y = np.log10(0.1 + dsnc.counts.loc[gene].values)
    ax.scatter(x, y, alpha=0.2, color='steelblue')
    ax.text(0.02, 0.98, '{:} = {:.2f}'.format(gene, corr.at[gene]),
            ha='left', va='top',
            transform=ax.transAxes,
            bbox=dict(facecolor='none', edgecolor='black', pad=1.1))

    ax.set_xlim(-1.1, 6)
    ax.set_ylim(-1.1, 6)
    ax.grid(True)
    ax.set_xlabel('vRNA [cpm]')
    ax.set_ylabel('{:} expression [cpm]'.format(gene))

    fig.tight_layout()

    plt.ion()
    plt.show()
