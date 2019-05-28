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


import logging
import argparse
from trimmer import Trimmer
from arguments import args_init, args_default


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
    parser.add_argument('-i', '--fq1', nargs='+', required=True, 
        help='reads in FASTQ files, support (*.gz), 1-4 files.')
    parser.add_argument('-o', '--path_out', default=None, 
        help='The directory to save results.')

    ## use pre-defined arguments
    parser.add_argument('--library-type', default=None, dest='library_type',
        choices=['rnaseq', 'chipseq', 'iclip', 'eclip', 'clipnsr', 'atacseq', 'smrna'],
        help='use pre-defined arguments for specific arguments')

    ## global arguments
    parser.add_argument('-m', '--len_min', default=15, metavar='len_min', 
        type=int, help='Minimum length of reads after trimming. (defualt: 15)')
    parser.add_argument('-a', '--adapter3', nargs='+', metavar='adapter', 
        type=str, default=['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'],         
        help='3-Adapter, \
        TruSeq RNAseq: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        TruSeq smallRNA: TGGAATTCTCGGGTGCCAAGG \
        Nextera: CTGTCTCTTATACACATCT')
    parser.add_argument('--update-name', action='store_false', dest='keep_name',
         help='if specified, add .clean to filename')
    parser.add_argument('--not-gzip', action='store_false', dest='gzip',
        help='if specified, gzip fastq file output')
    parser.add_argument('--threads', default=1, type=int, dest='threads',
        help='Number of threads to launch. (default: 1)')
    parser.add_argument('--overwrite', action='store_true', dest='overwrite',
        help='if spcified, overwrite exists file')

    ## trim options   
    parser.add_argument('-q', '--qual-min', default=20, type=int,
        dest='qual_min',
        help='The cutoff of base quality, default [20]')    
    parser.add_argument('-e', '--error-rate', default=0.1, type=float,
        dest='error_rate',
        help='Maximum allowed error rate, default [0.1]')
    parser.add_argument('-O', '--overlap', default=3, type=int,
        help='Required N bases overlap between reads and adapter. (default: 3)')
    parser.add_argument('-p', '--percent', default=80, type=int,
        help='minimum percent of bases that must have -q quality. (default: 80)')
    parser.add_argument('--trim-times', dest='trim_times', type=int,
        default=1, help='Trim adapter from reads by N times. (default: 1)')

    ## output
    parser.add_argument('--rm-untrim', action='store_true', dest='rm_untrim',
        help='if specified, discard reads without adapter')
    parser.add_argument('--save-untrim', action='store_true', dest='save_untrim',
        help='if specified, save untrimmed reads files to file')
    parser.add_argument('--save-too-short', action='store_true', 
        dest='save_too_short',
        help='if specified, save too-short reads files to file')
    parser.add_argument('--save-too-long', action='store_true', 
        dest='save_too_long',
        help='if specified, save too-long reads files to file')

    ## extra arguments
    parser.add_argument('--adapter-sliding', dest='adapter_sliding', 
        action='store_true', help='Create a list of adatpers using sliding \
        windows on adapter')
    parser.add_argument('--cut-after-trim', default='0', metavar='cut1', 
        dest='cut_after_trim',
        help='cut bases after trimming adapter, Number of bases to cut \
              from each read, plus on 5-prime end, minus on 3-prime end, \
              could be single, or double numbers, eg: [3 | -4 | 3,-4]. \
              (default [0])')
    parser.add_argument('--rmdup', action='store_true', dest='rmdup',
        help='if specified, remove duplicated reads' )
    parser.add_argument('--cut-after-rmdup', default='0', metavar='cut2', 
        dest='cut_after_rmdup',
        help='cut bases after rmdup, Number of bases to cut \
              from each read, plus on 5-prime end, minus on 3-prime end, \
              could be single, or double numbers, eg: [3 | -4 | 3,-4]. \
              (default [0])')
    parser.add_argument('--cut-to-length', default=0, metavar='cut3',
        dest='cut_to_length', type=int,
        help='cut reads from right, save the specific length of reads. \
              0=the full length. (default: [0])')

    ## PE arguments
    parser.add_argument('--fq2', nargs='+', default=None, 
        help='The read2 of pair-end reads')
    parser.add_argument('-A', '--AD3', nargs='+', type=str,
        default=['AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
        help='The 3 adapter of read2, default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
    parser.add_argument('-G', '--AD5', nargs='+',
        default=None,
        help='The 5 adapter of read1, default: None')

    parser.add_argument('-g', '--adapter5', nargs='+', default=None,
        help='5-Adapter, default: None')
    parser.add_argument('--not-trim-adapter', action='store_true', 
        dest='not_trim_adapter', help='Do not trim adapters, just trim ends by \
        specific length.')
    parser.add_argument('--keep-temp-files', action='store_true',
        dest='keep_temp_files', help='Save temporal files.')
    # parser.add_argument('--read12', type=int, default=1,
    #     help='which one of PE reads, 1=read1, 2=read2, default: 1')
    # parser.add_argument('--double-trim', action='store_true', 
    #     dest='double_trim', help='if specified, trim adapters twice')
    args = parser.parse_args()
    return args


def main():
    # args = vars(get_args()) # save as dictionary
    args = args_init(vars(get_args()), trim=True, align=False, call_peak=False) # save as dictionary
    fq1_files = args.pop('fq1', None) # remove  'fq1'

    ## update arguments
    args_lib = args_default(args['library_type'])
    args = {**args, **args_lib}

    logging.info('trimming start')
    ## SE mode
    if args['fq2'] is None:
        for fq1 in fq1_files:
            tmp = Trimmer(fq1=fq1, **args).trimmer()

    ## PE mode
    else:
        fq2_files = args.pop('fq2', None) # remove 'fq2'
        for fq1, fq2 in zip(fq1_files, fq2_files):
            tmp = Trimmer(fq1=fq1, fq2=fq2, **args).trimmer()

    logging.info('trimming finish')


if __name__ == '__main__':
    main()


## EOF
