# vim: fdm=indent
'''
author:     Fabio Zanini
date:       06/05/19
content:    Extract only viral reads from all bamfiles for SNP calling
'''
import os
import sys
import glob
import pysam


if __name__ == '__main__':

    fdn = '/oak/stanford/groups/quake/fzanini/sequencing_data/veev/'
    fns = glob.glob(fdn+'bamfiles/*/Aligned.out.bam')
    fn_header = fdn+'viral_reads/viral_reads.header.bam'
    fn_out = fdn+'viral_reads/viral_reads.bam'

    print('Make header file')
    with pysam.AlignmentFile(fns[0], 'rb') as bamfile0:
        with pysam.AlignmentFile(fn_header, 'wb', template=bamfile0) as bout:
            pass

    print('Merge viral reads')
    with pysam.AlignmentFile(fn_header, 'rb') as bheader:
        with pysam.AlignmentFile(fn_out, 'wb', template=bheader) as bout:
            n = len(fns)
            for i, fn in enumerate(fns):
                sn = fn.split('/')[-2]
                print('Sample {:}/{:}: {:}'.format(i+1, n, sn))
                with pysam.AlignmentFile(fn, 'rb') as bamfile:
                    for read in bamfile:
                        if read.reference_name == 'VEEVGFP':
                            read.qname = read.qname+':'+sn
                            bout.write(read)

    #print('Check output file')
    #with pysam.AlignmentFile(fn_out, 'rb') as bout:
    #    for read in bout:
    #        print(read.qname)
