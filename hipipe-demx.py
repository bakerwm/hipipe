#!/usr/bin/env python3

__author__ = 'Ming Wang'
__email__ = 'wangm08@hotmail.com'
__date__ = '2018-12-25'
__version__ = '0.3'

import argparse
import logging
from arguments import args_init
from demx import Demx

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

def get_args():
    """get args"""
    parser = argparse.ArgumentParser(
        prog = 'demultiplexing', description='demultiplexing Illumina reads',
        epilog='Example: ')
    parser.add_argument('-r1', '--fq1', required=True, dest='fq1', 
        help='SE read or read1 of PE reads')
    parser.add_argument('-r2', '--fq2', required=False, dest='fq2',
        default=None, 
        help='read2 of PE reads')
    parser.add_argument('--sample-info', required=True, 
        dest='sample_info',
        help='a list of sample names and index, barcode sequences\
        <p7> <barcode> <sample>')
    parser.add_argument('--demx-type', choices=['p7', 'barcode', 'both'],
        default='p7', dest='demx_type',
        help='Demultiplexing type, P7-index or inline-barcode, default: [p7]')
    parser.add_argument('-o', '--path-out', required=True, dest='path_out',
        help='The directory to save results.')
    parser.add_argument('--n-mismatch', dest='n_mismatch', default=0, 
        type=int,
        help='Number of mismatches allowed in barcode/index, default: 0')
    parser.add_argument('--bc-n-left', dest='bc_n_left', default=3, 
        metavar='N-LEFT', type=int,
        help='N-bases at the left of barcode. default: 3')
    parser.add_argument('--bc-n-right', dest='bc_n_right', default=2, 
        metavar='N-RIGHT', type=int,
        help='N-bases at the right of barcode. default: 2')
    parser.add_argument('--bc-in-read', dest='bc_in_read', choices=[1, 2],
        default=1, type=int, 
        help='barcode in read1/2 of PE reads, default: 1')
    parser.add_argument('--cut', action='store_true',
        help='cut the barcode from read if specified')
    args = parser.parse_args()
    return args


def main():
    logging.info('Demultiplexing: start')
    args = args_init(vars(get_args()), demx=True, trim=False, align=False, call_peak=False) # save as dictionary
    Demx(**args).run()
    logging.info('Delumtiplexing: finish')


if __name__ == '__main__':
    main()


## EOF
