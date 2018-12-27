#!/usr/bin/env python3

__author__ = 'Ming Wang'
__email__ = 'wangm08@hotmail.com'
__date__ = '2018-12-25'
__version__ = '0.3'

import os, sys
import argparse
import logging
from demx import Demx

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

def get_args():
    """get args"""
    parser = argparse.ArgumentParser(
        prog = 'demultiplexing', description='demultiplexing Illumina reads',
        epilog='Example: ')
    parser.add_argument('-r1', required=True, metavar='read1', 
        type=argparse.FileType('r'),
        help='SE read or read1 of PE reads')
    parser.add_argument('-r2', required=False, metavar='read2', 
        type=argparse.FileType('r'),  default=None, 
        help='read2 of PE reads')
    parser.add_argument('--sample-info', required=True, 
        dest='sample_info',
        help='a list of sample names and index, barcode sequences\
        <p7> <barcode> <sample>')
    parser.add_argument('--demx-type', choices=['p7', 'barcode', 'both'],
        default='p7', dest='demx_type',
        help='Demultiplexing type, P7-index or inline-barcode, default: [p7]')
    parser.add_argument('-o', '--out', required=True, metavar='OUTDIR', 
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
    args = get_args()
    read1 = args.r1.name
    read2 = args.r2.name if args.r2 else None
    sample_info = args.sample_info
    path_out = args.out
    demx_type = args.demx_type
    n_mismatch = args.n_mismatch
    bc_in_read = args.bc_in_read
    n_left = args.bc_n_left
    n_right = args.bc_n_right
    cut = args.cut # whether cut the barcode

    # Demx class
    demx = Demx(fq1=read1, 
        fq2=read2, 
        sample_info=sample_info, 
        path_out=path_out, 
        n_mismatch=n_mismatch, 
        bc_in_read = bc_in_read, 
        n_left = n_left, 
        n_right = n_right)

    #####################
    ## both p7 and bc  ##
    #####################
    if demx_type == 'both':
        # split p7 and barcode
        index_list = demx.sampleinfo_spliter() # [p7_list, bc_list1, bc_list2, ...]
        p7_list = index_list[0] # p7_list
        index_list.remove(p7_list)
        bc_list = index_list[0] # bc_list

        #####################
        ## run p7 1st      ##
        #####################
        path_p7 = os.path.join(path_out, 'p7_index')
        if not os.path.exists(path_p7):
            os.makedirs(path_p7)
        p7_info_file = os.path.join(path_p7, 'p7_sampleinfo.txt')
        with open(p7_info_file, 'wt') as fo:
            fo.write('\n'.join(p7_list) + '\n')

        demx_p7 = Demx(fq1=read1, 
            fq2=read2, 
            sample_info=p7_info_file, 
            path_out=path_p7, 
            n_mismatch=n_mismatch)

        logging.info('Demultiplexing: P7-index')
        if os.path.exists(read2):
            logging.info('PE mode')
            demx_p7.p7_demx_pe()
        else:
            logging.info('SE mode')
            demx_p7.p7_demx_se()

        #####################
        ## run barcode 2nd ##
        #####################
        for bc_sublist in bc_list:
            if len(bc_sublist) == 0:
                continue
            p7, bc, name = bc_sublist[0].split('\t')
            path_bc = os.path.join(path_out, 'barcode_%s' % p7)
            if not os.path.exists(path_bc):
                os.makedirs(path_bc)
            bc_info_file = os.path.join(path_bc, 'barcode_sampleinfo.txt')
            with open(bc_info_file, 'wt') as fo:
                fo.write('\n'.join(bc_sublist) + '\n')

            # output of P7-index
            bc_fq1 = os.path.join(path_p7, p7 + '_1.fq.gz')
            bc_fq2 = os.path.join(path_p7, p7 + '_1.fq.gz')

            demx_bc = Demx(fq1=bc_fq1, 
                fq2=bc_fq2, 
                sample_info=bc_info_file, 
                path_out=path_bc, 
                n_mismatch=n_mismatch, 
                bc_in_read = bc_in_read, 
                n_left = n_left, 
                n_right = n_right)

            logging.info('Demultiplexing: inline barcode, %s' % p7)
            if os.path.exists(bc_fq2):
                logging.info('PE mode')
                demx_bc.bc_demx_pe()
            else:
                logging.info('SE mode')
                demx_bc.bc_demx_se()

    ###################
    ## both p7 only  ##
    ###################
    elif demx_type == 'p7':
        logging.info('Demultiplexing: P7-index')
        if os.path.exists(read2):
            logging.info('PE mode')
            demx.p7_demx_pe()
        else:
            logging.info('SE mode')
            demx.p7_demx_se()

    ########################
    ## both barcode only  ##
    ########################
    elif demx_type == 'barcode':
        logging.info('Demultiplexing: inline barcode')
        if os.path.exists(read2):
            logging.info('PE mode')
            demx.bc_demx_pe()
        else:
            logging.info('SE mode')
            demx.bc_demx_se()
    else:
        raise Exception('either --demx-p7 or --demx-barcode is required')
    logging.info('Delumtiplexing: finish')


if __name__ == '__main__':
    main()


## EOF
