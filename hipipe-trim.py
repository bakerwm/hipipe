#!/usr/bin/env python3
#-*- encoding utf-8 -*-
"""
Trimming reads
## Required
1. 3' adapter, (default: TruSeq, optional)
2. low-quality (q=20) at 3' end
3. trim-n
4. remove sequences do not contain adapters (--rm-untrim)
## Optional
1. trim N-bases at either ends
2. trim sliding window of 3'-adapter
3. trim adapter by multiple times (--times=N)
4. remove PCR duplicates (--rm-PCR-dup)

to-do:
1. trim reads from 3', keep maxmium N-nt


## SE
3-adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

## PE
3-adapter1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
3-adapter2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT


"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


# import os
# import sys
# import re
# import shlex
# import subprocess
# import logging
import argparse
from trimmer import Trimmer
# from helper import *



def get_args():
    """
    processing SE read 
    - remove 3' adapter(s) (default: TruSeq RNA-Seq)
    - trim low-quality bases on both 5 and 3 end
    - trim N reads
    - cut N-bases at either end of read
    """
    parser = argparse.ArgumentParser(prog='trimmer', 
                                     description='trimming reads')
    parser.add_argument('-i', nargs='+', required=True, metavar='file', 
        type=argparse.FileType('r'),
        help='reads in FASTQ format, support (*.gz), 1-4 files.')
    parser.add_argument('-a', default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 
        metavar='adapter', type=str,
        help='3-Adapter, default: [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC].')
    parser.add_argument('-o', default=None, metavar='output', 
        help='The directory to save results.')
    parser.add_argument('-g', default='', metavar='adapter-5', type=str,
        help='5-Adapter, default: None')
    parser.add_argument('-m', default=15, metavar='len_min', 
        type=int, help='Minimum length of reads after trimming, defualt [15]')
    parser.add_argument('--read12', type=int, default=1, metavar='read12',
        help='which one of PE reads, 1=read1, 2=read2, default: 1')
    
    ## global arguments    
    parser.add_argument('-q', default=20, metavar='quality', type=int,
        help='The cutoff of base quality, default [20]')    
    parser.add_argument('-e', default=0.1, metavar='err_rate', type=float,
        help='Maximum allowed error rate, default [0.1]')
    parser.add_argument('-O', default=3, metavar='overlap', type=int,
        help='Required N bases overlap between reads and adapter, default [3]')
    parser.add_argument('-p', default=80, metavar='percent', type=int,
        help='minimum percent of bases that must have -q quality, default [80]')
    parser.add_argument('--rm-untrim', action='store_true', dest='rm_untrim',
        help='if specified, discard reads without adapter')
    parser.add_argument('--threads', default=1, metavar='threads', type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('--keep-name', action='store_true', dest='keep_name',
        help='if specified, do not change file names')
    
    ## extra arguments
    parser.add_argument('--adapter-sliding', dest='adapter_sliding', 
        action='store_true',
        help='Trim reads by sliding windows on adapter')
    parser.add_argument('--trim-times', dest='trim_times', type=int,
        default=1, help='Trim adapter from reads by N times, default:1')
    parser.add_argument('--double-trim', action='store_true', 
        dest='double_trim', help='if specified, trim adapters twice')
    parser.add_argument('--rm-dup', action='store_true', dest='rm_dup',
        help='if specified, remove duplicated reads' )
    parser.add_argument('--cut-before-trim', default='0', metavar='cut1', 
        dest='cut_before_trim',
        help='cut bases before trimming adapter, Number of bases to cut \
              from each read, plus on 5-prime end, minus on 3-prime end, \
              could be single, or double numbers, eg: 3 or -4 or 3,-4, \
              default [0]')
    parser.add_argument('--cut-after-trim', default='0', metavar='cut2', 
        dest='cut_after_trim',
        help='cut bases after trimming adapter, Number of bases to cut \
              from each read, plus on 5-prime end, minus on 3-prime end, \
              could be single, or double numbers, eg: 3 or -4 or 3,-4, \
              default [0]')
    parser.add_argument('--trim-to-length', default=0, metavar='max-length',
        dest='trim_to_length', type=int,
        help='trim reads from right, save the specific length of reads. \
              default: [0], 0=the full length')

    ## PE arguments
    parser.add_argument('--fq2', nargs='+', default=None, 
        help='The read2 of pair-end reads')
    parser.add_argument('-A', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
        help='The 3 adapter of read2, default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
    parser.add_argument('-G', default=None,
        help='The 5 adapter of read1, default: None')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    fq_files = [f.name for f in args.i]

    ## SE mode
    if args.fq2 is None: 
        for fq in fq_files:
            tmp = Trimmer(fq, args.a, args.o, args.m,
                adapter5=args.g, read12=args.read12, 
                qual_min=args.q, error_rate=args.e, overlap=args.O,
                rm_untrim=args.rm_untrim, threads=args.threads, 
                overwrite=args.overwrite, keep_name=args.keep_name,
                adapter_sliding=args.adapter_sliding, 
                trim_times=args.trim_times,
                double_trim=args.double_trim,
                rm_dup=args.rm_dup,
                cut_before_trim=args.cut_before_trim,
                cut_after_trim=args.cut_after_trim,
                trim_to_length=args.trim_to_length).run()
    ## PE mode
    else:
        # fq2_files = [f.name for f in args.fq2]
        fq2_files = args.fq2
        for fq1, fq2 in zip(fq_files, fq2_files):
            tmp = Trimmer(fq1, args.a, args.o, args.m,
                adapter5=args.g, read12=args.read12, 
                fq2=fq2, AD3=args.A, AD5=args.G,
                qual_min=args.q, error_rate=args.e, overlap=args.O,
                rm_untrim=args.rm_untrim, threads=args.threads, 
                overwrite=args.overwrite, keep_name=args.keep_name,
                adapter_sliding=args.adapter_sliding, 
                trim_times=args.trim_times,
                double_trim=args.double_trim,
                rm_dup=args.rm_dup,
                cut_before_trim=args.cut_before_trim,
                cut_after_trim=args.cut_after_trim,
                trim_to_length=args.trim_to_length).run()

if __name__ == '__main__':
    main()


## EOF
