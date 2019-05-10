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
    plt.subplots_adjust(0.05, 0.08, 0.95, 0.95)

    # Summary
    # position 9226 is polymorphic within our population AND across phylogeny
    # and is around 50 bases after the end of GFP. It is first codon position and makes for a nonsynonymous mutation Ala -> Thr at codon 292 (1-index) - 556 including GFP - of the structural polyprotein which starts at 7561 in this reference. So it is inside E3, a protein part of envelope
    # other polymorphic positions: 400, 1612, 1615, 1618, 2594, 8301, 10551, 11691

    # pos 400 is usually a T, but we have mostly G here and a C in the reference. It is also a synonymous mutation in a Leucin at third position, i.e. all 4 nucleotides are fine.

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



    plt.ion()
    plt.show()


