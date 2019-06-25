#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check contamination of input fastq files

bacteria classification

"""

import os
import argparse
import pathlib
from alignment_bacteria import Kraken2

def get_args():
    parser = argparse.ArgumentParser(prog='hipipe-bacteria.py', 
                                     description='check bacteria content of input files')
    parser.add_argument('-i', '--fq', nargs='+', required=True, dest='fq',
        type=str, help='fastq files')
    parser.add_argument('-o', '--out', default=None, dest='path_out',
        metavar='OUTPUT',  help='The directory to save results. (default:\
        current directory)')
    parser.add_argument('--db', required=True, help='Path to Kraken2 database')
    parser.add_argument('--top-n', default=10, type=int, dest='topN',
        help='Show topN species. (default: 10)')
    parser.add_argument('--save-out', dest='save_out', action='store_true',
	help='Save alignment output')
    parser.add_argument('--unmap-file', default=None, dest='unmap_file',
        help='Save unclassified reads to file')
    parser.add_argument('--kraken2', default=None, 
        help='Path to the Kraken2 script. default, search in $PATH')
    parser.add_argument('-p', '--threads', default=16, type=int, dest='threads',
        help='Number of processors to use. (default: 16)')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    ## output
    if args.path_out is None:
        args.path_out = str(pathlib.Path.cwd())
    ## run
    for fq in args.fq:
        k = Kraken2(fq, args.path_out, args.db, args.kraken2, args.threads, args.unmap_file, args.overwrite)
        k.report(args.topN)


if __name__ == '__main__':
    main()
