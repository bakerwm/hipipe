#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert bam to bigWig 
"""

import os
import pathlib
import argparse
import logging
import shlex
import subprocess

from helper import *


def get_args():
    parser = argparse.ArgumentParser(prog='hipipe-bam2bw', 
                                     description='convert bam to bigWig')
    parser.add_argument('-i', nargs='+', required=True, metavar='INPUT', 
        type=argparse.FileType('r'),
        help='bam files')
    parser.add_argument('-o', default=None, 
        metavar='OUTPUT',  help='The directory to save results, default:\
        current directory.')
    parser.add_argument('-g', required=True, default='hg19', 
        metavar='GENOME', choices=['dm3', 'hg19', 'hg38', 'mm10', 'mm9'],
        help='Reference genome : dm3, hg19, hg39, mm10, default: hg19')
    parser.add_argument('-s', metavar='strandness', type=int, default=1,
        help='Strandness, 0=not, 1=forward, 2=reverse, default: 0')
    parser.add_argument('-b', metavar='binsize', type=int, default=50,
        help='binsize of bigWig, default: 50')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    args = parser.parse_args()
    return args


## utilities
def bam2bigwig(bam, genome, path_out, strandness=0, binsize=1, overwrite=False):
    """Convert bam to bigWig using deeptools
    https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
    history:
    1. Mappable sequence of a genome, see Table 1 in 
       url: https://www.nature.com/articles/nbt.1518.pdf
    2. effective genome size:
        - non-N bases
        - regions (of some size) uniquely mappable
    3. UCSC
    http://genomewiki.ucsc.edu/index.php/Hg19_100way_Genome_size_statistics
    http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics

    !!! strandness
    default: dUTP-based library (read2 is sense strand, read1 is anti-sense strand)
    general RNA library: (NSR, small-RNA-library), read1 is sense, read2 is antisense
    """
    assert os.path.exists(bam)
    assert isinstance(genome, str)
    assert is_path(path_out)
    assert isinstance(strandness, int)
    assert isinstance(binsize, int)
    assert isinstance(overwrite, bool)
    bamcov = which('bamCoverage')
    if bamcov is None:
        raise ValueError('%10s | program not found: bamCoverage' % 'failed')

    effsize = {'dm3': 162367812,
               'dm6': 142573017,
               'mm9': 2620345972,
               'mm10': 2652783500,
               'hg19': 2451960000,
               'hg38': 2913022398,}
    gsize = effsize[genome]

    # prefix = os.path.basename(os.path.splitext(bam)[0])
    prefix = file_prefix(bam)[0]
    bw_log = os.path.join(path_out, prefix + '.deeptools.log')
    logging.info('create bigWig for: %s' % prefix)
    if strandness > 0:
        # strandness
        bw_fwd = os.path.join(path_out, prefix + '.fwd.bigWig')
        bw_rev = os.path.join(path_out, prefix + '.rev.bigWig')
        if strandness == 2:
            bw_fwd, bw_rev = [bw_rev, bw_fwd]
        
        # file existence
        if os.path.exists(bw_fwd) and os.path.exists(bw_rev) and overwrite is False:
            logging.info('bigWig file exists, skipped: %s' % prefix)
        else:
            # attention; bamCoverage using dUTP-based library
            # reverse, forward
            c1 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand forward \
                  --normalizeTo1x {}'.format(bam, bw_fwd, binsize, gsize)
            c2 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand reverse \
                  --normalizeTo1x {}'.format(bam, bw_rev, binsize, gsize)
            if os.path.exists(bw_fwd) and os.path.exists(bw_rev) and not overwrite:
                logging.info('file exists, bigWig skipped ...')
            else:
                with open(bw_log, 'wt') as fo:
                    p1 = subprocess.run(shlex.split(c1), stdout=fo, stderr=fo)
                    p2 = subprocess.run(shlex.split(c2), stdout=fo, stderr=fo)
            if not os.path.exists(bw_fwd) or not os.path.exists(bw_rev):
                raise ValueError('output file is missing, check log file: %s' % bw_log)
    else:
        # strandless
        bw = os.path.join(path_out, prefix + '.bigWig')
        if os.path.exists(bw) and overwrite is False:
            logging.info('bigWig file exists, skipping: %s' % prefix)
        else:
            c3 = 'bamCoverage -b {} -o {} --binSize {} \
                  --normalizeTo1x {}'.format(bam, bw, binsize, gsize)
            if os.path.exists(bw) and not overwrite:
                logging.info('file exists, bigWig skipped ...')
            else:
                with open(bw_log, 'wt') as fo:
                    subprocess.run(shlex.split(c3), stdout=fo, stderr=fo)
            if not os.path.exists(bw):
                raise ValueError('output file is missing, check log file: %s' % bw_log)



def main():
    args = get_args()
    # current dir
    if args.o is None:
        args.o = str(pathlib.Path.cwd())
    bam_files = [f.name for f in args.i]
    for bam_file in bam_files:
        bam2bigwig(bam=bam_file, genome=args.g, path_out=args.o, 
            strandness=args.s, binsize=args.b, overwrite=args.overwrite)

if __name__ == '__main__':
    main()


## EOF