#!/usr/bin/env python3
"""
Create alignment statistics report
"""

__author__ = 'Ming Wang'
__email__ = 'wangm08@hotmail.com'
__date__ = '2018-12-25'
__version__ = '0.2'

import logging
import argparse
from hipipe_reporter import Alignment_reporter

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

def get_args():
    """Alignment statistics
    """
    parser = argparse.ArgumentParser(prog='hipipe-align-stat', 
                                     description='mapping statistics')
    parser.add_argument('-i', '--input', nargs='+', required=True, 
        help='stat.csv file or directories contain stat.csv file')
    parser.add_argument('-o', '--output', default=None, 
        metavar='OUTPUT',  help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-t', '--template', default=None,
        help='the template Rmarkdown for align report')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    Alignment_reporter(args.input, args.output, args.template).run()


if __name__ == '__main__':
    main()


## EOF