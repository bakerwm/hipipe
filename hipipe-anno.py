#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
get annotation of bed

# homer annotation
$ hipipe-anno.py -i demo.bed -g dm3

# basic annotation
$ hipipe-anno.py -i demo.bed -g dm3 -t basic

# multiple BAM/BED files, save stat and pdf
$ hipipe-anno.py -i demo-*.bed -g dm3 -o demo.anno -f demo.pdf

"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


import os
import argparse
import shlex
import subprocess
import pandas as pd
from bed_annotation import Bed_anno



def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'bed_annotation',
        description = 'annotation for bed file',
        epilog = 'Example: hipipe-anno.py -i demo.bed -g dm3 -o demo.anno -f demo.pdf ')
    parser.add_argument('-i', nargs='+', required = True, metavar='BED', 
        type=argparse.FileType('r'),
        help = 'BED or BAM file(s)')
    parser.add_argument('-g', required = True, default = 'dm3', 
        metavar = 'GENOME', choices = ['dm3', 'hg19', 'hg38', 'mm10'],
        help = 'Reference genome : dm3, hg19, hg38, mm10, default: dm3')
    parser.add_argument('-o', metavar='output', default=None,
        help='Save results to file')
    parser.add_argument('-t', default = 'homer', 
        choices = ['homer', 'basic'], metavar = 'Type', 
        help = 'Type of the annotation, [basic|homer], default: homer')
    parser.add_argument('-f', dest='f', metavar='pdf', 
        default=None, help='Save plot to pdf flie')
    parser.add_argument('--path_data', default=None,
        help='The directory of genome files, default: $HOME/data/genome')
    args = parser.parse_args()
    return args



def ann2pdf(fn_anno, fn_pdf):
    """Create PDF for given annotation"""
    dir_pwd = os.path.dirname(os.path.realpath(__file__))
    a2p = os.path.join(dir_pwd, 'bed_annotation_plot.R')
    assert os.path.exists(a2p) #
    assert os.path.exists(fn_anno)
    assert isinstance(fn_pdf, str)

    c1 = 'Rscript %s %s %s' % (a2p, fn_anno, fn_pdf)
    p1 = subprocess.run(shlex.split(c1), stderr=subprocess.PIPE)
    if p1.stderr.decode() is '':
        logging.info('save plot to file: %s' % fn_pdf)
    return fn_pdf



def main():
    args = get_args()
    frames = []
    for fn in [f.name for f in args.i]:
        dx = Bed_anno(fn, genome=args.g, group=args.t, path_data=args.path_data).anno
        frames.append(dx)
    df = pd.concat(frames) # DataFrame
    if isinstance(args.o, str):
        fn_dir = os.path.abspath(os.path.dirname(args.o))
        if not os.path.exists(fn_dir):
            os.makedirs(fn_dir)
        Bed_anno(df, args.g).saveas(args.o)
        if isinstance(args.f, str):
            ann2pdf(args.o, args.f)
    else:
        print(df)



if __name__ == '__main__':
    main()



## EOF
