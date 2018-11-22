#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Mapping reads to reference genome
1. RNAseq (bowtie2, STAR)
2. filt unique reads (--unique-reads)
3. mismatch (0-2)
4. report unmap.fq, map.bam, map.log, map.json
"""


__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"

import os
import sys
import pathlib
import argparse
from alignment import Alignment, Alignment_log, Alignment_stat


def get_args():
    """
    Mapping SE read or one of PE reads to reference genome
    using bowtie, STAR, ... (universal)
    """
    parser = argparse.ArgumentParser(prog='aligner', 
                                     description='mapping reads')
    parser.add_argument('-i', nargs='+', required=True, metavar='INPUT', 
        type=argparse.FileType('r'),
        help='CLIP reads in FASTQ format, (not *.gz), 1-4 files.')
    parser.add_argument('-o', default=None, 
        metavar='OUTPUT',  help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-n', required=True, metavar='NAME',
        help='Name of the experiment')
    parser.add_argument('-g', required=True, default='hg19', 
        metavar='GENOME', choices=['dm3', 'hg19', 'hg38', 'mm10', 'mm9'],
        help='Reference genome : dm3, hg19, hg39, mm10, default: hg19')
    parser.add_argument('-k', default=None, 
        metavar='Spike-in', choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('-x', nargs='+', metavar='align_index',
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')
    parser.add_argument('--threads', default=8, 
        metavar='THREADS', type=int, 
        help='Number of threads to launch, default: 8.')
    parser.add_argument('--unique-only', action='store_true',
        dest='unique_only',
        help='if specified, keep unique mapped reads only')
    parser.add_argument('--n-map', dest='n_map', type=int, default=0,
        help='Report up to N alignments per read. use -k for bowtie and \
        bowtie2 (default 1), --outFilterMultimapNmax for STAR \
        (default 20).')
    parser.add_argument('--aligner', default='bowtie', 
        choices=['bowtie', 'bowtie2', 'STAR'],
        help='Choose which aligner to use. default: bowtie')
    parser.add_argument('--align-to-rRNA', dest='align_to_rRNA',
        action='store_true',
        help='if specified, align to rRNA before genome')
    parser.add_argument('--path_data', 
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    fqs = [f.name for f in args.i]
    if args.o is None:
        args.o = str(pathlib.Path.cwd())
    # tmp = Alignment(fqs, args.o, 
    #                 smp_name=args.n, 
    #                 genome=args.g,
    #                 spikein=args.k, 
    #                 index_ext=args.x,
    #                 multi_cores=args.threads, 
    #                 unique_only=args.unique_only,
    #                 n_map=args.n_map,
    #                 aligner=args.aligner,
    #                 align_to_rRNA=args.align_to_rRNA,
    #                 path_data=args.path_data,
    #                 overwrite=args.overwrite).run()
    tmp = Alignment(
        fqs, args.o, 
        smp_name=args.n,
        genome=args.g,
        spikein=args.k, 
        index_ext=args.x,
        threads=args.threads, 
        unique_only=args.unique_only,
        n_map=args.n_map,
        aligner=args.aligner,
        align_to_rRNA=args.align_to_rRNA,
        path_data=args.path_data,
        overwrite=args.overwrite).run()

if __name__ == '__main__':
    main()


# EOF