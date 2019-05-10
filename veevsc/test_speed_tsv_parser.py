# vim: fdm=indent
'''
author:     Fabio Zanini
date:       09/05/19
content:    Test the speed of a trivial TSV parser.
'''
import os
import sys
import numpy as np
import pandas as pd
import time


if __name__ == '__main__':


    fn = '../data/dataset/VEEVcounts.tsv'

    # We know the first row is cell names and the first
    # column is gene names, the rest is float32
    t0 = time.time()
    with open(fn, 'rt') as f:
        sns = f.readline().rstrip('\n')
        n = len(sns)
        gns = []
        for line in f:
            gns.append(line[:line.find('\t')])
        l = len(gns)
        f.seek(0)
        mat = np.empty((l, n), dtype=np.float32)
        f.readline()
        for il, line in enumerate(f):
            for ni, field in enumerate(line.split('\t')[1:-1]):
                mat[il, ni] = np.float32(field)
            mat[il, n - 1] = np.float32(field[:-1])

    t1 = time.time()
    print('Manual parser took: {:.1f} seconds'.format(t1 - t0))

    print('Not very fast')

