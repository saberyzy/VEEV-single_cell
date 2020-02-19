# vim: fdm=indent
'''
author:     Fabio Zanini
date:       10/11/19
content:    Plot 5' vs 3' VEEV genome coverage
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

    print('Loading host data...')
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

    print('Load viral allele frequencies')
    fn = '../data/dataset/viral_allele_frequencies_100000_reads.cdf'
    data = xr.open_dataset(fn)

    print('Restrict to data with host counts (1 sample missing)')
    data = data.sel({'sample': ds.samplenames})

    la, n, l = data['allele_frequencies'].shape

    if False:
        print('Coverage')
        cov = data['coverage']
        covm = np.log10(np.maximum(9e-1, cov.mean(axis=1).data))
        fig, ax = plt.subplots(
            1, figsize=(10, 3),
            )
        x = np.arange(l)
        ym = np.zeros(len(x))
        for j in range(n):
            covi = cov[j]
            y = np.log10(np.maximum(9e-1, covi.data))
            ym += y
            ax.plot(x, y, lw=1, alpha=0.05)#, label=covi.sample.to_dict()['data'])#, color=colors[j])
        ym /= (ym.mean() / 2.)
        ax.plot(x, ym, lw=1, alpha=1, color='red', label='Average [relative coverage]')

        ax.legend(loc='upper center')
        ax.grid(True)
        ax.set_xlabel('Position in VEEV-GFP [bases]')
        ax.set_ylabel('Coverage [# reads]')
        ax.set_yticks([0, 1, 2, 3, 4, 5])
        ax.set_yticklabels(['$0$', '$1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
        fig.tight_layout()
        fig.savefig('../figures/coverage_veev_with_average.pdf')
        fig.savefig('../figures/coverage_veev_with_average.svg')
        fig.savefig('../figures/coverage_veev_with_average.png', dpi=600)

    L = data['position'].shape[0]
    cov_host = ds.samplesheet['coverage']
    cov_virus = ds.samplesheet['VEEV_counts']
    cov_virus_host = cov_host + cov_virus
    cov_virus_aa = data['coverage'].sum(dim='position')

    bases = [
        [[0, 1700], [11000, L]],
        [[0, 350], [350, 700]],
        [[11000, 11700], [11700, L]],
        ]
    for bs1, bs2 in bases:
        cov_5p = data['coverage'].sel({'position': np.arange(bs1[0], bs1[1])}).mean(dim='position')
        cov_3p = data['coverage'].sel({'position': np.arange(bs2[0], bs2[1])}).mean(dim='position')

        if False:
            print('Compare viral coverage from allele counts and host dataset')
            fig, ax = plt.subplots()
            ax.scatter(
                cov_virus_aa,
                cov_virus,
                )
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(1, 1e7)
            ax.set_ylim(1, 1e7)
            ax.set_xlabel('Cov from AC')
            ax.set_xlabel('Cov from host')
            ax.grid(True)
            fig.tight_layout()

        from scipy.stats import gaussian_kde
        # We capped the allele count mapping to 100,000 reads
        def norm(xx):
            return 1.0 * xx / (0.1 + cov_virus_aa) * cov_virus / (0.1 + cov_virus + cov_host)
        x = np.log10(1e-5 + norm(cov_5p))
        y = np.log10(1e-5 + norm(cov_3p))
        dset = np.vstack([x, y])
        dset[:, np.isnan(dset).any(axis=0)] = -3
        dset = dset[:, (dset >= -4.8).any(axis=0)]
        kernel = gaussian_kde(dset)
        X, Y = np.mgrid[0:100:100j, 0:100:100j]
        xmax = dset.max()
        xmin = dset.min()
        xra = xmax - xmin
        X = (xra * X / 100.) + xmin
        Y = (xra * Y / 100.) + xmin
        positions = np.vstack([X.ravel(), Y.ravel()])
        Z = np.reshape(kernel(positions).T, X.shape)

        # Fit conditional probability until peak
        idx = Z.argmax()
        #xpeak = X.ravel()[idx]
        xbase = -4.8
        xpeak = -4.0
        tw = xpeak - xbase
        bw = 0.04
        nbins = int(tw / (2 * bw))
        bleft = np.linspace(xbase, xpeak - bw, nbins)
        bright = bleft + 0.5 * bw
        bcenter = 0.5 * (bleft + bright)
        h = []
        for i in range(nbins):
            #ind = (dset[0] >= bleft[i]) & (dset[0] <= bright[i])
            #hi = dset[1, ind].mean()
            ind = (X > bleft[i]) & (X <= bright[i])
            idx = Z.copy()
            idx[~ind] = 0
            idx = idx.argmax()
            hi = Y.ravel()[idx]
            h.append(hi)
        h = np.asarray(h)

        fig, ax = plt.subplots(figsize=(4, 3.5))
        ax.contourf(X, Y, Z, zorder=1, alpha=0.4, levels=30)
        ax.scatter(x, y, s=10, zorder=2, alpha=0.5)
        ax.plot(bcenter[:-1], h[:-1], lw=2, color='darkred')
        ax.arrow(
            bcenter[-2], h[-2], bcenter[-1] - bcenter[-2], h[-1] - h[-2],
            head_width=0.05,
            head_length=0.08,
            width=0.01,
            length_includes_head=True,
            edgecolor='none',
            facecolor='darkred',
            zorder=3,
            )
        xext = np.linspace(bcenter[-1], xpeak + 0.2 * (xpeak - xbase), 100)
        m = (h[-1] - h[-2]) / (bcenter[-1] - bcenter[-2])
        yext = h[-1] + m * (xext - xext[0])
        ax.plot(xext, yext, lw=2, ls='--', color='darkred', alpha=0.5)
        ax.plot([-5, 2 * xpeak - xbase], [-5, 2 * xpeak - xbase], lw=2, color='grey')
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(xmin, xmax)
        ax.set_xticks(np.log10(np.array([1e-4])))
        ax.set_yticks(np.log10(np.array([1e-4])))
        ax.set_xticks(np.log10(np.array([2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 2e-4])), minor=True)
        ax.set_yticks(np.log10(np.array([2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 2e-4])), minor=True)
        ax.set_xticklabels(['$10^{-4}$'])
        ax.set_yticklabels(['$10^{-4}$'])
        ax.grid(True)
        ax.set_xlabel("vRNA 5' coverage [normalized]")
        ax.set_ylabel("vRNA 3' coverage [normalized]")
        fig.tight_layout()


    fig, axs = plt.subplots(1, 3, figsize=(9, 3))
    bases = [
        [[0, 1700], [11000, L]],
        [[0, 900], [900, 1700]],
        [[11000, 11700], [11700, L]],
        ]
    for ib, (ax, (bs1, bs2)) in enumerate(zip(axs, bases)):
        cov_5p = data['coverage'].sel({'position': np.arange(bs1[0], bs1[1])}).mean(dim='position')
        cov_3p = data['coverage'].sel({'position': np.arange(bs2[0], bs2[1])}).mean(dim='position')

        r = ((0.1 + cov_3p) / (0.1 + cov_5p)).data
        cov_ercc = ds.counts.loc[ds.counts.index.str.startswith('ERCC-')].sum(axis=0)
        m = cov_virus / cov_ercc

        from scipy.stats import spearmanr
        ind = m >= 0.17
        x = np.log10(1e-2 + m[ind])
        y = r[ind]
        rho = spearmanr(x, y)[0]

        ax.scatter(x, y, s=10, zorder=2, alpha=0.25, label='$\\rho = {:.2f}$'.format(rho), color='darkred')
        sns.kdeplot(x, y, shade=True, zorder=1, alpha=0.3, cmap='viridis', shade_lowest=False, ax=ax)
        ax.set_ylim(0, [4, 1.5, 4][ib])
        ax.set_xlim(np.log10(0.17), 3)
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['$1$', '$10$', '$10^2$', '$10^3$'])
        ax.grid(True)
        ax.set_ylabel("3'/5' read ratio")
        ax.set_xlabel('vRNA/ERCC')
        ax.legend(loc='best')
        fig.tight_layout()

        fig.savefig('../figures/scatter_kde_35ratio.png')
        fig.savefig('../figures/scatter_kde_35ratio.pdf')
        fig.savefig('../figures/scatter_kde_35ratio.svg')

    plt.ion()
    plt.show()


    plt.ion()
    plt.show()
