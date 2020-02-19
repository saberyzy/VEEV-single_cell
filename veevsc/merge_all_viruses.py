# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/06/19
content:    Load various datasets and align features
'''
import os
import sys
import numpy as np
import scipy as sp
import pandas as pd
import loompy


data_dict = {
        'dengue': {
            'filename': 'dengue/dengue.loom',
            'host': 'human',
            'featurename': 'featurename',
            'vRNA': {
                'table': 'meta',
                'index': 'numberDengueReads',
                },
            },
        'zika': {
            'filename': 'Zika/zika.loom',
            'host': 'human',
            'featurename': 'featurename',
            'vRNA': {
                'table': 'meta',
                'index': 'numberZikaReads',
                },
            },
        'wnv': {
            'filename': 'wnv/wnv.loom',
            'host': 'mouse',
            'featurename': 'featurename',
            'vRNA': {
                'table': 'counts',
                'index': 'WNV',
                },
            },
        'flu': {
            'filename': 'flu/extreme/flu_extreme.loom',
            'host': 'human',
            'featurename': 'gene_long_name',
            'vRNA': {
                'table': 'meta',
                'index': ['flu-wt', 'flu-syn'],
                },
            },
        'mCMV': {
            'filename': 'mCMV/mcmv.loom',
            'host': 'mouse',
            'featurename': 'featurename',
            'vRNA': {
                'table': 'meta',
                'index': 'n_reads_virus',
                },
            },
        'veev': {
            'filename': 'veev/veev.loom',
            'host': 'human',
            'featurename': 'featurename',
            'vRNA': {
                'table': 'meta',
                'index': 'n_viral_reads',
                }
            }
        }


if __name__ == '__main__':

    print('Get genes and cell numbers')
    genes_by_virus = {}
    n_cells_by_virus = {}
    for virus, dd in data_dict.items():
        print(virus)
        fn = '../data/'+data_dict[virus]['filename']
        with loompy.connect(fn) as dsl:
            if dd['host'] == 'mouse':
                col = 'HumanGeneName'
            else:
                col = 'GeneName'
            genes_by_virus[virus] = dsl.ra[col][dsl.ra['has_homologue'].astype(bool)]
            n_cells_by_virus[virus] = dsl.shape[1]

    print('Intersect genes')
    from collections import Counter
    genes_intersect = None
    for genes in genes_by_virus.values():
        # Kill things that are there twice
        gc = Counter(genes)
        genes = [g for g in genes if gc[g] == 1]
        if genes_intersect is None:
            genes_intersect = genes
        else:
            genes_intersect = np.intersect1d(genes_intersect, genes)

    print('Fill output matrices')
    n_genes = len(genes_intersect)
    n_cells = sum(n_cells_by_virus.values())
    counts = np.zeros((n_genes, n_cells), np.float32)
    cells = np.zeros(n_cells, 'U50')
    cell_virus = np.zeros(n_cells, 'U10')
    nc = 0
    for virus, dd in data_dict.items():
        print(virus)
        fn = '../data/'+data_dict[virus]['filename']
        with loompy.connect(fn) as dsl:
            if dd['host'] == 'mouse':
                col = 'HumanGeneName'
            else:
                col = 'GeneName'
            ser = pd.DataFrame(
                    np.vstack([dsl.ra[col], np.arange(dsl.shape[0])]).T,
                    columns=['GeneName', 'index']).set_index('GeneName')['index']
            ind = ser.loc[genes_intersect].values.astype(int)
            tmp = np.asarray(dsl[:, :])
            counts[:, nc: nc+dsl.shape[1]] = tmp[ind]
            cells[nc: nc+dsl.shape[1]] = [virus+'_'+sn for sn in dsl.ca['samplename']]
            cell_virus[nc: nc+dsl.shape[1]] = virus
            nc += dsl.shape[1]

    print('Get metadata on number of viral reads')
    n_viral_reads = np.zeros(n_cells, np.float32)
    nc = 0
    for virus, dd in data_dict.items():
        print(virus)
        fn = '../data/'+data_dict[virus]['filename']
        vrna = data_dict[virus]['vRNA']
        with loompy.connect(fn) as dsl:
            if vrna['table'] == 'meta':
                if isinstance(vrna['index'], str):
                    inds = [vrna['index']]
                else:
                    inds = vrna['index']
                for ind in inds:
                    n_viral_reads[nc: nc+dsl.shape[1]] += dsl.ca[ind]
            else:
                ind = list(dsl.ra['featurename']).index(vrna['index'])
                n_viral_reads[nc: nc+dsl.shape[1]] = dsl[ind, :]
            nc += dsl.shape[1]

    print('Save merged output')
    fn_out = '../data/all/all_viruses.loom'
    meta_genes = {
            'featurename': genes_intersect,
            }
    meta_cells = {
            'samplename': cells,
            'virus': cell_virus,
            'n_viral_reads': n_viral_reads,
            }
    loompy.create(fn_out, counts, meta_genes, meta_cells)
