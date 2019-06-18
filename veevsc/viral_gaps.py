# vim: fdm=indent
'''
author:     Fabio Zanini
date:       21/05/19
content:    Check reads mapped against VEEV with N cigars (intron)
'''
import os
import sys
import numpy as np
import pandas as pd
import pysam

import matplotlib.pyplot as plt


if __name__ == '__main__':

    fdn = '../data/dataset/'
    fn_bam = fdn+'viral_reads_with_gaps_10000.bam'

    # Python convention, 0-index and left excluded
    # the first 4 bases before the gap are often mismapped
    gap = [10762, 10798]

    # That part is very GC-rich, check for RNA structures
    from Bio import SeqIO
    seq = SeqIO.read(fdn+'VEEV_reference.fasta', 'fasta')

    print('GAP reference:')
    print(str(seq.seq[gap[0]: gap[1]]))
    # It's pretty clearly an RNA hairpin


    gap2 = [274, 319]
    print('GAP2 reference:')
    print(str(seq.seq[gap2[0]: gap2[1]]))


    print('Check random stretches')
    import subprocess as sp
    rfdn = '/home/fabio/programs/university/RNAstructure_cli/'
    my_env = os.environ.copy()
    my_env['DATAPATH'] = rfdn+'data_tables/'
    nreps = 1000
    res = []
    for rep in range(nreps):
        print('{:}\r'.format(rep+1))
        gap_ctrl = [np.random.randint(len(seq))]
        gap_ctrl.append(gap_ctrl[0] + gap[1] - gap[0])
        gapseq = str(seq.seq[gap_ctrl[0]: gap_ctrl[1]])

        with open('/tmp/test.fasta', 'wt') as f:
            f.write('>test\n'+gapseq+'\n')

        sp.run(
            rfdn+'exe/Fold /tmp/test.fasta /tmp/test.ct',
            shell=True,
            env=my_env,
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL,
            )
        with open('/tmp/test.ct', 'rt') as f:
            header = f.readline().rstrip('\n')
            if 'ENERGY' in header:
                fe = float(header.split()[3])
            else:
                fe = None

        #sp.run(
        #    rfdn+'exe/Fold /tmp/test.ct /tmp/test.efn2',
        #    shell=True,
        #    env=my_env,
        #    )
        #with open('/tmp/test.efn2', 'rt') as f:
        #    fe = float(f.read().rstrip('\n'))

        res.append({
            'sequence': gapseq,
            'gap': gap_ctrl,
            'free_energy': fe,
            })
    print()
    res = pd.DataFrame(res)

    fig, ax = plt.subplots()
    bins = np.array([-25, -22.5, -20, -17.5, -15, -12.5, -10, -7.5, -5, -2.5, 0.01])
    h = np.histogram(res['free_energy'].fillna(0).values, bins=bins)[0]
    ax.bar(bins[:-1], h, bins[1:]-bins[:-1], align='edge')
    ax.set_xlabel('Free energy [kcal/mol]')
    ax.grid(True)
    ax.set_title('Free energies of random stretches in VEEV genome')

    fig.tight_layout()
    plt.ion()
    plt.show()
