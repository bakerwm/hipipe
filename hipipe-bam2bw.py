#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert bam to bigWig 
"""

import sys
import os
import pathlib
import argparse
import logging
import shlex
import subprocess
from arguments import args_init
from helper import BAM, bam2bigwig, bam2bigwig2, bam2bw


def get_args():
    parser = argparse.ArgumentParser(prog='hipipe-bam2bw', 
                                     description='convert bam to bigWig')
    parser.add_argument('-b', '--bam', nargs='+', required=True, dest='bam',
        type=str, help='bam files')
    parser.add_argument('-o', '--out', default=None, dest='path_out',
        metavar='OUTPUT',  help='The directory to save results. (default:\
        current directory)')
    parser.add_argument('-g', '--genome', required=False, default=None, 
        dest='genome', choices=[None, 'dm3', 'dm6', 'hg19', 'hg38', 'mm10', 'mm9'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm10. (default: hg19)')
    parser.add_argument('-s', '--scale', nargs='+', required=False, \
        default=['1.0'], type=float, dest='scale', 
        help='The computed scaling factor will be multiplied by this. \
        (default: 1.0)')
    parser.add_argument('--library-type', dest='library_type',
        choices=[1, 2], type=int, default=1, help='The library type, \
        1: minus strand reads map on gene forward strand; dUTP, NSR, ... \
        2: plus strand reads map on gene forwrad strand; smRNA, CLIP, ... \
        (default: 1)')
    parser.add_argument('--normalizeUsing', default=None,
        choices=[None, 'RPKM', 'CPM', 'BPM', 'RPGC'],
        help='Use one of the entered methods to normalize the number of \
        reads per bin. By default, no normalization is performed. RPKM = Reads \
        Per Kilobase per Million mapped reads; CPM = Counts Per Million mapped \
        reads; BPM = Bins Per Million mapped reads; RPGC = reads per genomic \
        content (1x normalization); None = the default and equivalent to not \
        setting this option at all. default: None')
    parser.add_argument('--filterRNAstrand', required=False, 
        dest='filterRNAstrand', choices=[None, 'forward', 'reverse', 'both'],
        help='This is the option bamCoverage command from deeptools: \
        Selects RNA-seq reads originating from genes on the given \
        strand. This option assumes a strandard dUTP-based library preparation \
        (that is, --filterRNAstrand=forward keeps minus-strand reads, which \
        originally came from genes on the forward strand using a dUTP-based \
        methods). Consider using --samExcludeFlag instead for filtering by \
        strand in other contexts. (default: None)')
    parser.add_argument('--samFlagExclude', required=False, 
        dest='samFlagExclude', type=int, default=None,
        help='Exlude reads based on the SAM flag. For example, to get only\
        reads that map to the forward strand, use --samFlagExclude 16, \
        which is equal to samtools view -F, where 16 is the SAM flag for \
        reads that map to the reverse strand. (default: None)')
    parser.add_argument('--samFlagInclude', required=False,
        dest='samFlagInclude', type=int, default=None,
        help='Include reads based on the SAM falg. For example, to get only\
        reads that are the first mate, using a flag of 64. This is useful to \
        count properly paired reads only once, as otherwise the second mate \
        will be also considered for the coverage. (default: None)')
    parser.add_argument('--binsize', metavar='binsize', type=int, 
        default=10, dest='binsize', help='binsize of bigWig. (default: 10)')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('-p', type=int, default=1,
        help='Number of processors to use. (default: 1)')
    args = parser.parse_args()
    return args


def main():
    ## prepare arguments
    args = args_init(vars(get_args()), bam2bw=True)

    # print(args)
    if isinstance(args['bam'], str):
        if isinstance(args['scale'], float):
            pass
        else:
            raise Exception('--bam, --scale, not in the same length')
    elif isinstance(args['bam'], list):
        if isinstance(args['scale'], list) and len(args['scale']) == 1:
            pass
        elif isinstance(args['scale'], list) and len(args['bam']) == len(args['scale']):
            pass
        else:
            raise Exception('--bam, --scale, not in the same length')
    else:
        pass

    ## functions
    bam_files = args.pop('bam', None)
    path_out = args.pop('path_out', None)
    scale_list = args.pop('scale', None)    
    if len(scale_list) == 1 and len(bam_files) > 1:
        scale_list = scale_list * len(bam_files)
    elif len(scale_list) == len(bam_files):
        pass
    else:
        sys.exit('lenght of --scale and -b are not equal')
 #   scale_list = args.pop('scale', None)

    ## program
    for bam_file in bam_files:
        i = bam_files.index(bam_file)
        scale = scale_list[i]
        ## strand
        if args['filterRNAstrand'] == 'both':
            ## forward strand
            args['filterRNAstrand'] = 'forward'
            bam2bw(bam_file, path_out, scale, **args)
            ## reverse strand
            args['filterRNAstrand'] = 'reverse'
            bam2bw(bam_file, path_out, scale, **args)
            ## reset
            args['filterRNAstrand'] = 'both'
        else:
            bam2bw(bam_file, path_out, scale, **args)



if __name__ == '__main__':
    main()


## EOF
