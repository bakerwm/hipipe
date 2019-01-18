#!/usr/bin/env python3
"""
Create fastqc report for multiple fastq files
"""

__author__ = 'Ming Wang'
__email__ = 'wangm08@hotmail.com'
__date__ = '2018-12-25'
__version__ = '0.2'


import argparse
from hipipe_reporter import QC_reporter


def get_args():
    """
    Mapping SE read or one of PE reads to reference genome
    using bowtie, STAR, ... (universal)
    """
    parser = argparse.ArgumentParser(prog='aligner', 
                                     description='mapping reads')
    parser.add_argument('-i', '--input', nargs='+', required=True, 
        help='fastq files or directories contain fastq files')
    parser.add_argument('-o', '--output', default=None, 
        metavar='OUTPUT',  help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-t', '--template', default=None,
        help='the template Rmarkdown for qc report')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    QC_reporter(args.input, args.output, args.template).run()


if __name__ == '__main__':
    main()


## EOF