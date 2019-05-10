# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/05/19
content:    Explore the viral coverage
'''
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == '__main__':

    print('Load viral allele frequencies')
    fn = '../data/dataset/viral_allele_frequencies_100000_reads.cdf'
    data = xr.open_dataset(fn)

    la, n, l = data['allele_frequencies'].shape

    print('Coverage')
    cov = data['coverage']
    covm = np.log10(np.maximum(9e-1, cov.mean(axis=1).data))
    fig, axs = plt.subplots(
        1, 2, figsize=(22, 4),
        gridspec_kw={'width_ratios': [10, 1]},
        sharey=True,
        )
    ax = axs[0]
    for j in range(n):
        covi = cov[j]
        x = np.arange(l)
        y = np.log10(np.maximum(9e-1, covi.data))
        ax.plot(x, y, lw=1, alpha=0.05, label=covi.sample.to_dict()['data'])#, color=colors[j])
    ax.grid(True)
    ax.set_xlabel('Position in VEEV-GFP [bases]')
    ax.set_ylabel('Coverage [# reads]')
    ax.set_yticks([0, 1, 2, 3, 4, 5])
    ax.set_yticklabels(['$0$', '$1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])

    # Smoothed A-richness
    ax2 = ax.twinx()
    for kl in (11, 21, 51, 101, 301):
        kernel = np.ones(kl)
        shift = (len(kernel) - 1) // 2
        ref_at = (data.reference.data == 'A')# | (data.reference.data == 'T')
        ref_at_smooth = 1.0 * np.convolve(ref_at, kernel, mode='valid') / kernel.sum()
        x = np.arange(l - len(kernel) + 1) + shift
        y = ref_at_smooth
        ax2.plot(x, y, lw=1, alpha=0.3, color='k', label='A-richness [w={:}]'.format(len(kernel)))
    ax2.set_ylim(bottom=0)

    ax = axs[1]
    bins = np.array([-0.1, 0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7, 4])
    h = np.histogram(covm, bins=bins, density=False)[0]
    ax.barh(bins[:-1], width=h, height=bins[1:] - bins[:-1], align='edge', left=0.9)
    ax.set_xlim(left=0.8)
    ax.set_xscale('log')
    ax.set_xticks([1, 10, 100, 1000])
    ax.grid(True)
    ax.set_xlabel('cells with mean coverage')
    fig.tight_layout()

    print('Correlate coverage and AT richness')
    from scipy.stats import spearmanr
    for kl in (5, 7, 11, 21, 51, 101, 151, 201, 251, 301, 351, 401, 451, 501):
        kernel = np.ones(kl)
        shift = (len(kernel) - 1) // 2
        ref_at = (data.reference.data == 'A') | (data.reference.data == 'T')
        ref_at_smooth = 1.0 * np.convolve(ref_at, kernel, mode='valid') / kernel.sum()
        x = np.arange(l - len(kernel) + 1) + shift
        y = ref_at_smooth
        covs_smooth = np.convolve(cov.sum(axis=0), kernel, mode='valid') / kernel.sum()
        rho = spearmanr(covs_smooth[5000:], ref_at_smooth[5000:])
        print('Kernel size: {:}, rho = {:.2f}, P = {:.2e}'.format(kl, rho[0], rho[1]))

    a = []
    ref_at = (data.reference.data == 'A')# | (data.reference.data == 'T')
    for kl in (7, 11, 21, 51, 101, 151, 201, 251, 301, 351, 401, 451, 501):
        kernel = np.ones(kl)
        shift = (len(kernel) - 1) // 2
        ref_at_smooth = 1.0 * np.convolve(ref_at, kernel, mode='valid') / kernel.sum()
        x = np.arange(l - len(kernel) + 1) + shift
        y = ref_at_smooth
        covs_smooth = np.convolve(cov.sum(axis=0), kernel, mode='valid') / kernel.sum()
        for sh in (-500, -100, -20, -10, 1, 10, 20, 50, 100, 200, 500):
            if sh > 0:
                rho = spearmanr(covs_smooth[4000:-sh], ref_at_smooth[4000 + sh:])
            else:
                rho = spearmanr(covs_smooth[4000 - sh:], ref_at_smooth[4000:sh])
            a.append({'kl': kl, 'sh': sh, 'rho': rho[0]})
            #print('Kernel size: {:}, rho shifted = {:.2f}, P = {:.2e}'.format(kl, rho[0], rho[1]))
    a = pd.DataFrame(a).set_index(['kl', 'sh']).unstack()
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(a, ax=ax, cmap=sns.diverging_palette(240, 10, as_cmap=True), vmin=-0.6, vmax=0.6)
    fig.tight_layout()

    a = []
    for kl in (11, 51, 151, 251, 351, 451):
        for sh in (-3000, -2500, -2000, -1500, -1000, -500, 1, 100, 200, 500, 1000, 1200, 1500, 1750, 2000, 2250, 2500, 2750, 3000):
            kernel = np.ones(kl)
            shift = (len(kernel) - 1) // 2
            ref_a = data.reference.data == 'A'
            ref_a_smooth = 1.0 * np.convolve(ref_a, kernel, mode='valid') / kernel.sum()
            x = np.arange(l - len(kernel) + 1) + shift
            y = ref_at_smooth
            covs_smooth = np.convolve(cov.sum(axis=0), kernel, mode='valid') / kernel.sum()
            if sh > 0:
                rho = spearmanr(covs_smooth[4000:-sh], ref_a_smooth[4000 + sh:])
            else:
                rho = spearmanr(covs_smooth[4000 - sh:], ref_a_smooth[4000:sh])
            a.append({'kl': kl, 'sh': sh, 'rho': rho[0]})
            print('Kernel size: {:}, rho shifted = {:.2f}, P = {:.2e}'.format(kl, rho[0], rho[1]))
    a = pd.DataFrame(a).set_index(['kl', 'sh']).unstack()
    fig, ax = plt.subplots(figsize=(6, 3))
    sns.heatmap(a, ax=ax, cmap=sns.diverging_palette(240, 10, as_cmap=True), vmin=-0.6, vmax=0.6)
    fig.tight_layout()

    print('Correlate coverage in a block with A richness in another block of size 20, downstream')
    kl_a = 20
    kernel_a = np.ones(kl_a)
    ref_a = data.reference.data == 'A'
    ref_a_smooth = 1.0 * np.convolve(ref_a, kernel_a, mode='valid') / kernel_a.sum()
    # beginning of A-rich block
    x_a_smooth = np.arange(len(ref_a_smooth))

    pad = 5

    kl = 200
    kernel = np.ones(kl)
    covs_smooth = np.convolve(cov.sum(axis=0), kernel, mode='valid') / kernel.sum()
    # end of coverage block
    x_smooth = np.arange(len(covs_smooth)) + kl

    rho = spearmanr(covs_smooth[:-(kl_a + pad)], ref_a_smooth[kl + pad:])
    print(rho)

    print('Consider averaged, normalized coverage')
    covh = cov[covm > 2]
    covhnor = covh / covh.sum(axis=1)
    y = covhnor.mean(axis=0)
    x = y.position
    ys = np.convolve(y, np.ones(10)) / 10.
    xs = np.arange(len(ys))
    fig, ax = plt.subplots(figsize=(24, 3))
    ax.plot(xs, ys[ys>0].min() * 0.5 + ys, label='<norm coverage>', lw=2, alpha=0.7)
    ax.set_yscale('log')
    ax.set_ylabel('Coverage (blue)')
    ax.set_title('A richness in viral genome versus coverage')
    ax.set_xlabel('Position in VEEV-GFP genome')

    ax2 = ax.twinx()
    kl_a = 20
    kernel_a = np.ones(kl_a)
    ref_a = data.reference.data == 'A'
    ref_a_smooth = 1.0 * np.convolve(ref_a, kernel_a, mode='valid') / kernel_a.sum()
    x_a_smooth = np.arange(len(ref_a_smooth))
    ax2.plot(x_a_smooth, ref_a_smooth, color='k', label='A richness', lw=2, alpha=0.05)
    ind = ref_a_smooth >= 0.5
    ax2.scatter(x_a_smooth[ind], ref_a_smooth[ind], color='k', label='A richness', lw=2, alpha=0.7)
    ax2.set_ylabel('A richness, w=20 (black + dots)')

    fig.tight_layout()


    plt.ion()
    plt.show()
