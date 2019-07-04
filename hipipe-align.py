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
from helper import *

def get_args():
    """
    Mapping SE read or one of PE reads to reference genome
    using bowtie, STAR, ... (universal)
    """
    parser = argparse.ArgumentParser(prog='aligner', 
                                     description='mapping reads')
    parser.add_argument('-i', '--fq1', nargs='+', required=True,
        help='Reads in FASTQ format, 1-4 files.')
    parser.add_argument('-I', '--fq2', nargs='+', required=False,
        default=None,
        help='The second mate of paired end reads in fastq format.')
    parser.add_argument('-o', '--path-out', default=None, dest='path_out',
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-n', '--smp-name', dest='smp_name',
        help='Name of the experiment')
    parser.add_argument('-g', '--genome', required=True, default='hg19', 
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10', 'GRCh38'],
        help='Reference genome : dm3, hg19, hg39, mm10, default: hg19')
    parser.add_argument('-k', '--spikein', default=None, 
        choices=[None, 'dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10', 'GRCh38'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('-x', '--extra-index', nargs='+', dest='extra_index',
        help='Provide extra alignment index(es) for alignment, support multiple\
        indexes. eg. Transposon, tRNA, rRNA and so on. if specified, ignore -g, -k')
    parser.add_argument('--align-to-te', dest='align_to_te', action='store_true',
        help='Align to TE consensus sequences (only for dm3), ignore -g, -k, -x')
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
    parser.add_argument('--small-genome', dest='small_genome',
        action='store_true',
        help='if align to small genome, over 90 percent of the reads are not aligned,\
        should use this option for STAR alignment, to reduce the alignment time.')
    parser.add_argument('--simple-name', dest='simple_name', 
        action='store_true',
        help='use simple name for fastq prefix, remove .not_..., .map_to..., ')
    parser.add_argument('--log-level', default='INFO', dest='log_level',
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)

    return args


def main():
    args = args_init(vars(get_args()), align=True) # save as dictionary

    log.info('aaaaa')

    # args['align_to_te'] = True

    ## run alignment
    map_bam_list = Alignment(**args).run()

    # logging.info('aaaaaa')

    # # specific arguments
    # args['align_to_rRNA'] = True # force mapping to rRNA
    # if args['not_merge_rep']:
    #     args['merge_rep'] = False # for replicates

    # tmp1 = Alignment(**args).run()
    # tmp2 = Alignment_reporter(input=args['path_out'], output=args['path_out']).run()


if __name__ == '__main__':
    main()


# EOF
