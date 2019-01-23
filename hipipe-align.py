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
from arguments import args_init
from alignment import Alignment, Alignment_log, Alignment_stat
from hipipe_reporter import Alignment_reporter


def get_args():
    """
    Mapping SE read or one of PE reads to reference genome
    using bowtie, STAR, ... (universal)
    """
    parser = argparse.ArgumentParser(prog='aligner', 
                                     description='mapping reads')
    parser.add_argument('-i', '--fq1', nargs='+', required=True,
        help='CLIP reads in FASTQ format, (not *.gz), 1-4 files.')
    parser.add_argument('-o', '--path_out', default=None, 
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-n', '--smp_name', required=True,
        help='Name of the experiment')
    parser.add_argument('-g', '--genome', required=True, default='hg19', 
        choices=['dm3', 'hg19', 'hg38', 'mm10', 'mm9'],
        help='Reference genome : dm3, hg19, hg39, mm10, default: hg19')
    parser.add_argument('-k', '--spikein', default=None, 
        choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('-x', '--index_ext', nargs='+',
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')
    parser.add_argument('--threads', default=8, type=int, 
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
    parser.add_argument('--not-merge-replicate', action='store_false',
        dest='not_merge_rep',
        help='if specified, do not merge replicates')
    parser.add_argument('--repeat-masked-genome', dest='repeat_masked_genome',
        action='store_true',
        help='map to repeat masked reference genome, data from EnsEMBL')
    parser.add_argument('--genome-path', dest='genome_path',
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    args = parser.parse_args()
    return args


def main():
    args = args_init(vars(get_args()), trim=False, align=True, call_peak=False) # save as dictionary

    # specific arguments
    args['align_to_rRNA'] = True # force mapping to rRNA
    if args['not_merge_rep']:
        args['merge_rep'] = False # for replicates

    tmp1 = Alignment(**args).run()
    tmp2 = Alignment_reporter(input=args['path_out'], output=args['path_out']).run()


if __name__ == '__main__':
    main()


# EOF