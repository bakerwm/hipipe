#!/usr/bin/env python3

"""Extract N-bases from both ends of SE reads, generate the "Paired-end" reads
Date: 2019-06-08
"""

import os
import sys
import re
import gzip
import subprocess
import argparse
import logging

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout,
    level = logging.DEBUG)
log = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(prog='te_ins.py',
        description='')
    parser.add_argument('-i', type=str, required=True,
        help='FASTQ file')
    parser.add_argument('-o', required=True,
        help='Output directory')
    parser.add_argument('-n', nargs='+', type=int, default=[25],
        help='Length of the PE reads, default: 25')
    parser.add_argument('--prefix', default=None,
        help='The prefix of new PE files, \
        defualt, the same as input fastq file')
    args = parser.parse_args()
    return args


def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def extract_fq(fq_in, fq_out, n=25, tail=False):
    """Extract N-nt from reads"""
    n, slen, qlen = 0, 0, 0
    # open file
    if fq_in.endswith('.gz'):
        f1open = gzip.open
    else:
        f1open = open

    # open file
    if fq_out.endswith('.gz'):
        f2open = gzip.open
    else:
        f2open = open

    # suffix
    suffix = '_1'
    if tail:
        suffix = '_2'
    suffix = ''

    flag = 0
    with f1open(fq_in, 'rt') as fi, f2open(fq_out, 'wt') as fo:
        for name, seq, qual in readfq(fi):
            nlen = len(seq)
            if nlen < 2 * n:
                continue # skip
            if tail:
                s1 = seq[:n]
                q1 = qual[:n]
            else:
                s1 = seq[-n:]
                q1 = qual[-n:]
            flag += 1
            fo.write('\n'.join(['@' + str(flag) + suffix, s1, '+', q1]) + '\n')


def se2pe(fq_in, outdir, n_base=25, gzipped=True, prefix=None):
    """Generate PE reads from SE reads"""
    n, slen, qlen = 0, 0, 0
    # open file
    if fq_in.endswith('.gz'):
        f1open = gzip.open
    else:
        f1open = open

    ## fq anme
    fq_name = re.sub('.gz$', '', os.path.basename(fq_in))
    fq_name = os.path.splitext(fq_name)[0]
    if prefix:
        fq_name = prefix

    ## output
    fq1 = os.path.join(outdir, fq_name + '_1.fq')
    fq2 = os.path.join(outdir, fq_name + '_2.fq')
    f2open = open

    if gzipped:
        fq1 += '.gz'
        fq2 += '.gz'
        f2open = gzip.open

    flag = 0
    with f1open(fq_in, 'rt') as fi, f2open(fq1, 'wt') as f1, f2open(fq2, 'wt') as f2:
        for name, seq, qual in readfq(fi):
            nlen = len(seq)
            if nlen < 2 * n_base:
                continue # skip
            # fq1
            s1 = seq[:n_base]
            q1 = qual[:n_base]
            # fq2
            s2 = seq[-n_base:]
            q2 = qual[-n_base:]
            # print(n)
            if flag % 1000000 == 0:
                log.info('%d reads processed' % flag)
            flag += 1
            f1.write('\n'.join(['@' + str(flag) + ' 1', s1, '+', q1]) + '\n')
            f2.write('\n'.join(['@' + str(flag) + ' 2', s2, '+', q2]) + '\n')


def main():
    args = parse_arguments()

    for n in args.n:
        log.info('cut bases: ' + str(n))
        subout = os.path.join(args.o, 'cut' + str(n))
        if not os.path.exists(subout):
            os.makedirs(subout)
        se2pe(args.i, subout, n_base=n, prefix=args.prefix)


if __name__ == '__main__':
    main()

