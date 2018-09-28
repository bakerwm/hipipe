#!/usr/bin/env python
"""
filter bad records in BED file

# bug
print dDataFrame to to screen, formats

"""

import os
import sys
import io
import argparse
import logging
import pathlib
import tempfile
import binascii
import gzip
import numpy as np
import pandas as pd

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)


def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'bed_fixer',
        description = 'filter illegal bed records',
        epilog = 'Example: bed_fixer.py -o out.bed input.bed')
    parser.add_argument('bed', nargs='?', type=argparse.FileType('r'),
        default=sys.stdin, metavar='bed',
        help='BED file or or STDIN')
    parser.add_argument('-o', 
        default=sys.stdout, metavar='output',
        help='save output BED records, default: STDOUT')
    args = parser.parse_args()
    return args


def is_gz(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def file_row_counter(fn):
    """
    count the file rows
    count '\n' 
    from @glglgl on stackoverflow, modified
    https://stackoverflow.com/a/9631635/2530783
    """
    def blocks(files, size = 1024 * 1024):
        while True:
            b = files.read(size)
            if not b: break
            yield b
    freader = gzip.open if is_gz(fn) else open
    with freader(fn, 'rt', encoding="utf-8", errors='ignore') as fi:
        return sum(bl.count('\n') for bl in blocks(fi))


class Bed_parser(object):


    def __init__(self, bed):
        """
        parsing BED records from file
        """
        self.bed = bed
        if isinstance(bed, Bed_parser):
            self.bed = bed.bed
        elif isinstance(bed, pd.DataFrame):
            self.bed = bed
        elif isinstance(bed, io.TextIOWrapper):
            self.bed = self._bed_parser()
        elif os.path.exists(bed):
            self.bed = self._bed_parser()
        else:
            raise ValueError('not supported file')


    def _tmp(self):
        """
        Create a temp file
        """
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.tmp',
                                            delete=False)
        tmpfn = tmpfn.name
        return tmpfn


    def _collapse(self, df=None, fn=None):
        """
        Collapses an DataFrame into file fn
        Returns the newly created filename.
        """
        if fn is None:
            fn = self._tmp()

        if df is None:
            df = self.bed

        default_kwargs = dict(sep='\t', header=False, index=False)
        if isinstance(fn, io.TextIOWrapper):
            print(df.to_string(index=False, header=False, justify='left'))
        else:
            df.to_csv(fn, **default_kwargs)
        return fn


    def saveas(self, _out=None):
        """
        Make a copy of the BED records
        """
        if _out is None:
            _out = self._tmp()

        _out = self._collapse(fn=_out)
        return _out



    def count(self):
        """
        count number of records
        """
        if isinstance(self.bed, Bed_parser):
            df = self.bed.bed
        elif isinstance(self.bed, pd.DataFrame):
            df = self.bed
        else:
            raise ValueError('unknown type of values')
        n = len(df.index)
        return n


    
    def _bed_parser(self, usecols=None):
        """
        read BED file as pandas DataFrame
        select specific columns, default all, (None)
        require at least 6 columns
        """
        bed = self.bed
        df = pd.DataFrame(columns = ['chr', 'start', 'end', 'name', 'score',
                                     'strand'])
        if isinstance(bed, Bed_parser):
            return bed
        elif isinstance(bed, pd.DataFrame):
            return Bed_parser(bed)
        elif isinstance(bed, io.TextIOWrapper):
            df = pd.read_table(bed, '\t', usecols=usecols, header=None,
                    dtype = {'0': np.str, '1': np.int64, '2': np.int64, \
                             '3': np.str, '4': np.int64, '5': np.str})
            df = df.rename(index = str, columns = {0: 'chr', 1: 'start', \
                           2: 'end', 3: 'name', 4: 'score', 5: 'strand'})
            return df
        elif os.path.exists(bed):
            if file_row_counter(bed) == 0:
                logging.error('empty bed file: %s' % bed)
            else:
                df = pd.read_table(bed, '\t', usecols=usecols, header=None) #,
                        # dtype = {'0': np.str, '1': np.int64, '2': np.int64, '3': np.str, '4': np.int64, '5': np.str})
                df = df.rename(index = str, columns = {0: 'chr', 1: 'start', 2: 'end', 3: 'name', 4: 'score', 5: 'strand'})
            return df
        else:
            raise ValueError('unknown values')



    def bed_fixer(self, usecols=None):
        """
        filt BED records
        1. start, end both are int
        2. start < end
        """
        bed = self.bed
        if isinstance(bed, Bed_parser):
            df = bed.bed
        elif isinstance(bed, pd.DataFrame):
            df = bed
        elif bed is None:
            bed = self.bed
            if isinstance(bed, pd.DataFrame):
                df = bed
            else:
                df = self._bed_parser().bed
        else:
            raise ValueError('unknown values')

        dx = df[['start', 'end']].apply(pd.to_numeric)
        a = dx['start'] >= 0
        b = dx['end'] > 0
        c = dx['start'] < dx['end']
        v = a.multiply(b).multiply(c)
        df_filted = df.loc[v, :]
        return Bed_parser(df_filted)
        


def main():
    args = get_args()
    # b = Bed_parser(args.bed).bed_fixer().saveas(args.o)
    a = Bed_parser(args.bed)
    b = a.bed_fixer()
    b.saveas(args.o)
    logging.error('input %s, report: %s (%.2f%%)' % (a.count(), b.count(), b.count()/a.count()*100))

if __name__ == '__main__':
    main()

## EOF
