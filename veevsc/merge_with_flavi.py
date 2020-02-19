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

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
import singlet


if __name__ == '__main__':

    print('Load biomart gene name conversion')
    fn_biomart = '../data/dengue_Zika/mart_export.tsv'
    bm = pd.read_csv(fn_biomart, sep='\t', index_col=0).to_dict()['Gene name']

    print('Load VEEV')
    fn_veev = '../data/dataset/veev.loom'
    with loompy.connect(fn_veev) as dsl:
        samplenames = dsl.ca['samplename']
        featurenames = dsl.ra['featurename']

        ss = pd.DataFrame([], index=samplenames)
        for col in dsl.ca:
            ss[col] = dsl.ca[col]

        shape = dsl.shape
        mat = np.zeros(shape, np.float32)
        mat[:, :] = dsl.layers[''][:, :]

    ss_veev = singlet.SampleSheet(ss)
    ss_veev['virus'] = 'VEEV'
    ss_veev['virus_reads'] = ss_veev['VEEV_counts']
    ss_veev['virus_reads_per_million'] = 1e6 * ss_veev['virus_reads'] / (ss_veev['virus_reads'] + ss_veev['coverage'])
    counts = singlet.CountsTable(
        pd.DataFrame(
            mat,
            index=featurenames,
            columns=samplenames,
        ))
    counts._spikeins = featurenames[-103:-7]
    counts._otherfeatures = featurenames[-7:]
    counts.name = 'counts'
    counts_veev = counts.normalize()
    del counts

    print('Load DENV')
    counts_denv = singlet.CountsTable(pd.read_csv(
        '../data/dengue_Zika/counts_dengue.tsv',
        sep='\t',
        index_col=0,
        ).astype(np.float32))
    ss_denv = pd.read_csv(
        '../data/dengue_Zika/cell_metadata_dengue.tsv',
        sep='\t',
        index_col='name',
        )
    ss_denv['virus_reads'] = ss_denv['numberDengueReads']
    ss_denv['virus'] = 'DENV'
    ss_denv['coverage'] = counts_denv.values[:-102].sum(axis=0)
    ss_denv['virus_reads_per_million'] = 1e6 * ss_denv['virus_reads'] / (ss_denv['virus_reads'] + ss_denv['coverage'])
    fs_denv = pd.DataFrame([], index=counts_denv.index)
    fs_denv.index.name = 'EnsemblID'
    fs_denv['GeneName'] = [bm.get(x, x) for x in fs_denv.index]
    counts_denv._spikeins = fs_denv.index[-102:-6]
    counts_denv._otherfeatures = fs_denv.index[-6:]
    counts_denv.normalize(inplace=True)
    fs_denv = fs_denv.loc[counts_denv.index]

    print('Load ZIKV')
    counts_zika = singlet.CountsTable(pd.read_csv(
        '../data/dengue_Zika/counts_Zika.tsv',
        sep='\t',
        index_col=0,
        ).astype(np.float32))
    ss_zika = pd.read_csv(
        '../data/dengue_Zika/cell_metadata_Zika.tsv',
        sep='\t',
        index_col='name',
        )
    ss_zika['virus_reads'] = ss_zika['numberZikaReads']
    ss_zika['virus'] = 'ZIKV'
    ss_zika['coverage'] = counts_zika.values[:-102].sum(axis=0)
    ss_zika['virus_reads_per_million'] = 1e6 * ss_zika['virus_reads'] / (ss_zika['virus_reads'] + ss_zika['coverage'])
    fs_zika = pd.DataFrame([], index=counts_zika.index)
    fs_zika.index.name = 'EnsemblID'
    fs_zika['GeneName'] = [bm.get(x, x) for x in fs_zika.index]
    counts_zika._spikeins = fs_zika.index[-102:-6]
    counts_zika._otherfeatures = fs_zika.index[-6:]
    counts_zika.normalize(inplace=True)
    fs_zika = fs_zika.loc[counts_zika.index]

    print('Convert flavi into gene names')
    from collections import Counter
    genescou = Counter(fs_denv['GeneName'].values)
    genesu = np.sort([g for g, c in genescou.items() if c == 1])

    genes_int = np.intersect1d(counts_veev.index, genesu)

    ensids = fs_denv.index[fs_denv['GeneName'].isin(genes_int)]
    mat_denv = counts_denv.loc[ensids]
    mat_denv.index = fs_denv['GeneName'][ensids]
    mat_zika = counts_zika.loc[ensids]
    mat_zika.index = fs_zika['GeneName'][ensids]

    # Align all
    mat_denv = mat_denv.loc[genes_int]
    mat_zika = mat_zika.loc[genes_int]
    mat_veev = counts_veev.loc[genes_int]

    L = len(genes_int)
    N1 = mat_denv.shape[1]
    N2 = mat_zika.shape[1]
    N3 = mat_veev.shape[1]
    N = N1 + N2 + N3
    columns = mat_denv.columns.tolist() + mat_zika.columns.tolist() + mat_veev.columns.tolist()
    mat = np.zeros((L, N), np.float32)
    mat[:, :N1] = mat_denv.values
    mat[:, N1: N1 + N2] = mat_zika.values
    mat[:, N1 + N2:] = mat_veev.values
    mat = pd.DataFrame(mat, index=genes_int, columns=columns).fillna(0).astype(np.float32)

    cols = ['virus', 'virus_reads', 'virus_reads_per_million', 'coverage', 'MOI', 'time [h]']
    ss = pd.concat([ss_denv[cols], ss_zika[cols], ss_veev[cols]])
    ss['CellID'] = ss.index
    cols += ['CellID']

    fn_merge = '../data/dengue_Zika/all3.loom'
    loompy.create(
        fn_merge,
        layers={'': mat.values},
        row_attrs={'GeneName': mat.index.values},
        col_attrs={key: ss[key].values for key in cols},
        )
