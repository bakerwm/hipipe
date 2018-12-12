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
from helper import BAM, bam2bigwig, bam2bigwig2


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
    parser.add_argument('-c', dest='scale', metavar='scale', type=float, default=1,
        help='Normalization scale, default: 1')
    parser.add_argument('-b', metavar='binsize', type=int, default=50,
        help='binsize of bigWig, default: 50')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    # current dir
    if args.o is None:
        args.o = str(pathlib.Path.cwd())
    bam_files = [f.name for f in args.i]
    for bam_file in bam_files:
        # bam2bigwig(bam=bam_file, genome=args.g, path_out=args.o, 
        #     strandness=args.s, binsize=args.b, overwrite=args.overwrite)
        if args.scale == 0:
            cnt = int(BAM(bam_file).count())
            args.scale = 1e6 / cnt
        # scale
        # bam2bigwig2(bam, path_out, scale=1, binsize=1, overwrite=False)
        bam2bigwig2(bam=bam_file, path_out=args.o, scale=args.scale, binsize=args.b,
            overwrite=args.overwrite)

if __name__ == '__main__':
    main()


## EOF