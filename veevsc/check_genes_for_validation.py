# vim: fdm=indent
'''
author:     Fabio Zanini
date:       18/06/19
content:    Check gene list from Zhiyuan for validation
'''
import os
import sys
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
import singlet


def load_gene_list():
    df = pd.read_excel('../data/Zhiyuan/Final_list.xlsx', 'Sheet1')
    return df


if __name__ == '__main__':

    print('Load gene list')
    gene_list = list(load_gene_list()['genename'].values)
    # Rename a few
    gene_list[gene_list.index('SURF4')] = ('SURF4_1', 'SURF4')

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

    print('Get average expression')
    ge_mean = ds.counts.get_statistics(metrics=('mean',)).iloc[:, 0]

    print('Exclude controls')
    dsnc = ds.query_samples_by_metadata('(MOI > 0) & (VEEV_counts >= 2)')

    print('Correlate in cells with 2+ vRNA reads')
    corr = dsnc.correlation.correlate_features_phenotypes(['virus_reads_per_million']).fillna(0)

    print('Load 5p and 3p viral coverage')
    cov53 = pd.read_csv('../data/dataset/vRNA_5p_3p.tsv', sep='\t', index_col=0)

    print('Correlate host gene expression with both')
    dsh = ds.query_samples_by_name(cov53.index.values)
    dsh.samplesheet['cov_5p'] = cov53['vRNA_5p']
    dsh.samplesheet['cov_3p'] = cov53['vRNA_3p']
    # NOTE: this is already normalized as number of viral reads in those regions
    # over 1e6 uniquely mapped host reads
    corr53 = dsh.correlation.correlate_features_phenotypes(
            phenotypes=['cov_5p', 'cov_3p']).fillna(0)
    corr['5p'] = corr53['cov_5p']
    corr['3p'] = corr53['cov_3p']

    print('Plot correlations')
    fig, axs = plt.subplots(2, 1, figsize=(10, 4.8), sharey=True)
    tiers = [np.arange(13), np.arange(13, 26)]
    colors = sns.color_palette('colorblind', n_colors=3)
    for ir, (ax, ind) in enumerate(zip(axs, tiers)):
        genes = [gene_list[i] for i in ind]
        titles = []
        for ig, gene in enumerate(genes):
            if isinstance(gene, str):
                gene, title = gene, gene
            else:
                gene, title = gene
            titles.append(title)
            datum = corr.loc[gene]
            for ic, cat in enumerate(['virus_reads_per_million', '5p', '3p']):
                left = ig + 0.25 * ic - 0.25
                width = 0.25
                height = datum[cat]
                if ig == 0:
                    ax.bar(left, height, width, color=colors[ic], label=cat)
                else:
                    ax.bar(left, height, width, color=colors[ic])

        ax.set_xticks(np.arange(13))
        ax.set_xticklabels(titles)
        ax.set_yticks([-0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4])
        ax.set_xlim(-0.5, 12.5)
        ax.grid(True)
        if ir == 0:
            ax.legend(loc='upper right')
    fig.suptitle('Correlations between gene expression and vRNA - VEEV final gene list')
    fig.tight_layout(rect=(0, 0, 1, 0.98))

    plt.ion()
    plt.show()
