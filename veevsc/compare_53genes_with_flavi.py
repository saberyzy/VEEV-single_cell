# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/11/19
content:    Compare VEEV with flavis
'''
import os
import sys
import numpy as np
import pandas as pd
import loompy

import ternary
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
import singlet


if __name__ == '__main__':

    print('Load VEEV + DENV + ZIKA')
    ds0 = singlet.Dataset(dataset={
        'path': '../data/dengue_Zika/all3.loom',
        'index_samples': 'CellID',
        'index_features': 'GeneName',
        'bit_precision': 32,
        'format': 'loom',
        })
    ds = ds0.query_samples_by_metadata('virus_reads_per_million >= 2')

    print('Find genes that start at similar levels and diverge')
    dsp = ds.split('virus')

    print('Load correlations for VEEV 5/3 ratio')
    fn_ratio = '../data/Zhiyuan/VEEV_correlation_53ratio.tsv'
    corr_53 = pd.read_csv(fn_ratio, sep='\t', index_col=0)

    print('Calculate correlations with DENV/ZIKV')
    corrp = {
        vir: dsp[vir].correlation.correlate_features_phenotypes('virus_reads_per_million').fillna(0)
        for vir in ['DENV', 'ZIKV']
        }

    print('Find overlapping genes')
    genes_ovl = np.array(list(set(corr_53.index.values) & set(corrp['DENV'].index.values)))

    print('Make dataframe of correlations')
    df = corr_53.loc[genes_ovl, ['spearman_ratio_3p_5p']]
    df.rename(columns={'spearman_ratio_3p_5p': 'VEEV_35'}, inplace=True)
    df['DENV'] = corrp['DENV'].loc[genes_ovl]
    df['ZIKV'] = corrp['ZIKV'].loc[genes_ovl]

    print('Calculate ranks')
    L = len(df)
    for key in ['VEEV_35', 'DENV', 'ZIKV']:
        df.sort_values(key, inplace=True, ascending=False)
        df['rank_{:}'.format(key)] = np.arange(L) + 1
    for key in ['VEEV_35', 'DENV', 'ZIKV']:
        df.sort_values(key, inplace=True, ascending=False)
        df['enrich_{:}'.format(key)] = (200 * (-np.linspace(0, 1, L) + 0.5)).astype(int)
    for key in ['VEEV_35', 'DENV', 'ZIKV']:
        df.sort_values(key, inplace=True, ascending=False)
        df['percentile_{:}'.format(key)] = 100 * (1.0 - np.linspace(0, 1, L))
    df = df.loc[genes_ovl]

    print('Plot percentiles')
    ngenes = [30, 200]
    fig1, axs1 = plt.subplots(1, 2, figsize=(6, 3), sharex=True, sharey=True)
    fig2, axs2 = plt.subplots(1, 2, figsize=(6, 3), sharex=True, sharey=True)
    axs = [axs1, axs2]
    for ir, ax_row in enumerate(axs):
        genes = df.nlargest(ngenes[ir], 'VEEV_35').index
        genes_random = df.index[np.random.randint(L, size=ngenes[ir])]
        for ax, genesi in zip(ax_row, [genes, genes_random]):
            if ax == ax_row[0]:
                s = 108 - (10*np.linspace(0, 1, len(genesi)))**2
            else:
                s = 30

            x, y = df.loc[genesi, ['percentile_DENV', 'percentile_ZIKV']].values.T

            if ir == 0:
                ax.scatter(x, y, s=s, color='grey', alpha=0.6)
            else:
                sns.kdeplot(x, y, bw=5, cmap='Greys', ax=ax)

            if (ir == 0) and (ax == ax_row[0]):
                for ig, gene in enumerate(genesi[:10]):
                    xi = x[ig]
                    yi = y[ig]
                    ax.scatter([xi], [yi], s=5, color='black', alpha=0.8)
                    ax.text(xi, yi, gene, ha='center', va='center')

            ax.set_xlim(-1, 101)
            ax.set_ylim(-1, 101)
            ax.grid(True)
            ax.set_xlabel('DENV percentile')
            if ax == ax_row[0]:
                ax.set_ylabel('ZIKV percentile')

        ax_row[0].set_title('Top VEEV 3/5 ratio genes')
        ax_row[1].set_title('Random genes')
    fig1.tight_layout()
    fig2.tight_layout()

    #fig1.savefig('../figures/comparison_VEEV35ratio_flavis.png')
    #fig1.savefig('../figures/comparison_VEEV35ratio_flavis.pdf')
    #fig1.savefig('../figures/comparison_VEEV35ratio_flavis.svg')
    #fig2.savefig('../figures/comparison_VEEV35ratio_flavis_kde.png')
    #fig2.savefig('../figures/comparison_VEEV35ratio_flavis_kde.pdf')
    #fig2.savefig('../figures/comparison_VEEV35ratio_flavis_kde.svg')

    print('Compare the bottom left corner with the cluster 3 from figure 5')
    genes_lowerleft = df.loc[(df['percentile_VEEV_35'] > 90) & (df['percentile_DENV'] < 10) & (df['percentile_ZIKV'] < 10)].index.values
    genes_cl3 = pd.read_csv('../data/all_viruses/gene_cluster_3.csv', sep='\t', header=None, squeeze=True).values

    plt.ion()
    plt.show()

