#!/usr/bin/env python
"""
Convert fasta, bed to GTF format
"""



import os
import sys
import pysam
import pybedtools

# f = '/home/data/genome/dm3/dm3_transposon/dm3.transposon.fa.fai'

def bed2gtf(x):
    """Convert BED to GTF
    source: BED
    feature: exon
    gene_id:
    """
    assert isinstance(x, list)
    des = 'gene_id "%s"; gene_name "%s";' % (x[3], x[3])
    start = str(int(x[1]) + 1)
    tabs = [x[0], 'bed', 'exon', start, x[2], '.', x[5], '.', des]
    return tabs



# with open(f, 'rt') as ff:
#     for line in ff:
#         name, length = line.strip().split('\t')[0:2]
#         fb, label = name.split('_')
#         tabs = [name, '0', length, label, '100', '+']
#         gtf = bed2gtf(tabs)
#         print('\t'.join(gtf))

f = sys.argv[1]

with open(f, 'rt') as ff:
    for line in ff:
        p = line.strip().split('\t')[0:6]
        gtf = bed2gtf(p)
        print('\t'.join(gtf))


