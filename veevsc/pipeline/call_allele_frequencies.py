# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/05/19
content:    Call allele frequencies for VEEV reads.
'''
import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import xarray as xr
from Bio import SeqIO
import pysam


if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='Call alleles from viral reads')
    pa.add_argument('--maxreads', type=int, default=-1)

    args= pa.parse_args()

    fdn = '/oak/stanford/groups/quake/fzanini/sequencing_data/veev/'

    print('Alphabet')
    alphal = ['A', 'C', 'G', 'T', '-', 'N']
    alpha = np.array(alphal)
    la = len(alpha)

    print('Load reference sequence')
    fn_ref = fdn+'genome/VEEVGFP.fa'
    refseq = SeqIO.read(fn_ref, 'fasta')
    l = len(refseq)

    print('Load samplenames')
    fns = glob.glob(fdn+'bamfiles/*/Aligned.out.bam')
    sns = [x.split('/')[-2] for x in fns]
    n = len(sns)

    print('Prepare output data')
    mat = np.zeros((la, n, l), dtype=np.float32)
    reads_per_cell = pd.Series(np.zeros(n), index=sns)
    old_sn = ''

    print('Call allele frequencies from merged BAM file')
    fn_bam = fdn+'viral_reads/viral_reads.bam'
    with pysam.AlignmentFile(fn_bam, 'rb') as bamfile:
        for ir, read in enumerate(bamfile):
            if ((ir + 1) % 1000) == 0:
                print(ir+1, end='\r')

            # Find sample the read is coming from
            sn = read.qname.split(':')[-1]
            ni = sns.index(sn)

            reads_per_cell[sn] += 1
            if (args.maxreads != -1) and (reads_per_cell[sn] >= args.maxreads):
                continue

            if sn != old_sn:
                print(sn)
                old_sn = sn

            # Find positions and alleles throughout the whole read
            poss = read.get_aligned_pairs()
            qual = read.query_qualities
            seq = read.query_sequence
            for pos_read, pos_ref in poss:
                # Deletion
                if pos_read is None:
                    mat[-1, ni, pos_ref] += 1
                    continue

                # Ignore insertions for now
                if pos_ref is None:
                    continue

                # Ignore low quality
                if qual[pos_read] < 30:
                    continue

                # Match/mismatch
                mat[alphal.index(seq[pos_read]), ni, pos_ref] += 1

    print('Wrap matrix in Dataset')
    acs = xr.DataArray(
        mat,
        coords={'nucleotide': alpha, 'sample': sns, 'position': np.arange(l)},
        dims=['nucleotide', 'sample', 'position'],
        name='allele_counts',
        )
    afs = acs / (acs.sum(axis=0) + 1e-10)
    afs.name = 'allele_frequencies'
    ref = xr.DataArray(
        np.array(refseq),
        coords={'position': np.arange(l)},
        dims=['position'],
        )

    data = xr.Dataset({
        'allele_counts': acs,
        'allele_frequencies': afs,
        'reference': ref,
        })

    print('Save to file')
    if args.maxreads == -1:
        fn = fdn+'viral_reads/allele_frequencies.cdf'
    else:
        fn = fdn+'viral_reads/allele_frequencies_{:}_reads.cdf'.format(args.maxreads)
    data.to_netcdf(fn)
