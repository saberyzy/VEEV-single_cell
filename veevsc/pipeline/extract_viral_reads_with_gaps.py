# vim: fdm=indent
'''
author:     Fabio Zanini
date:       21/05/19
content:    Check gaps within reads
'''
import os
import sys
import pysam


if __name__ == '__main__':

    #n_max = 1000000
    n_max = 'all'

    fdn = '/oak/stanford/groups/quake/fzanini/sequencing_data/veev/viral_reads/'
    fn_bam = fdn+'viral_reads.bam'
    fn_bam_gaps = fdn+'viral_reads_with_gaps_{:}.bam'.format(n_max)
    with pysam.AlignmentFile(fn_bam, 'rb') as f:
        n_reads = sum(1 for read in f)

    with pysam.AlignmentFile(fn_bam, 'rb') as f:
        with pysam.AlignmentFile(fn_bam_gaps, 'wb', template=f) as fout:
            n_with_gaps = 0
            for ir, read in enumerate(f):
                if ((ir+1) % 10000) == 0:
                    print('Read: {:}/{:} read pairs'.format(ir+1, n_reads))

                if 'N' in read.cigarstring:
                    print(read.cigarstring)
                    fout.write(read)
                    n_with_gaps += 1

                    if (n_with_gaps % 100) == 0:
                        print('Reads with gaps: {:}'.format(n_with_gaps))

                    if n_with_gaps == n_max:
                        break
