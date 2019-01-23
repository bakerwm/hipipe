#!/usr/bin/env python3

"""
Processing unitas output

1. summary statistics
2. miRNA expression
3. tRF expression

"""

import os
import re
from functools import reduce
import argparse
import pandas as pd
from helper import is_path


def get_args():
    parser = argparse.ArgumentParser(prog='unitas_run', 
                                     description='working with unitas')
    parser.add_argument('-i', nargs='+', required=True, metavar='unitas_path',
        help='path of unitas output')
    parser.add_argument('-o', required=True, metavar='output',
        help='directory to save results')
    args = parser.parse_args()
    return args


class Unitas(object):
    """Parsing the output of unitas
    small RNA annotation, quantification
    """

    def __init__(self, path_out):
        """the output directory"""
        self.path_out = path_out


    def check_unitas(self):
        """Require 4 files"""



    def leading_spaces(self, x):
        """Count the leading spaces"""
        return len(x) - len(x.lstrip())


    def summary(self, depth=1):
        """annotation summary
        level, the group levels
        only save depth=1 counts
        """
        f = os.path.join(self.path_out, 'unitas.annotation_summary.txt')
        if os.path.exists(f):
            dd = {}
            with open(f, 'rt') as fi:
                for line in fi:
                    level = int(self.leading_spaces(line) / 3)
                    if level > 0:
                        continue
                    group, count = line.strip().split('\t')
                    group = re.sub(' ', '_', group)
                    dd[group] = dd.get(group, int(float(count)))
            ##
            df = pd.DataFrame.from_dict(dd, orient='index', columns=['count'])
            df = df.transpose()
            return df
        else:
            return None


    def full_annotation_matrix(self):
        """Read full annotation matrix as pandas"""
        f = os.path.join(self.path_out, 'unitas.full_annotation_matrix.txt')
        if os.path.exists(f):
            pass
        else:
            return None


    def hits_per_target(self):
        """Read matrix for hits per targets"""
        f = os.path.join(self.path_out, 'unitas.hits_per_target.txt')
        header = ['class', 'name', 'count']
        dtype = {'class': 'str', 'name': 'str', 'count': 'float'}
        if os.path.exists(f):
            df = pd.read_csv(f, sep='\t', header=None, skiprows=1, 
                names=header, dtype=dtype)
            return df
        else:
            return None


    def miR_mod(self):
        """unitas.miR-modifications_Human.txt"""
        f = os.path.join(self.path_out, 'unitas.miR-modifications_Human.txt')

        with open(f, 'rt') as fi:
            for line in fi:
                if ' ' in line:
                    # head line
                    group = line.strip().split(' ')[0]
                    group = re.sub('[^a-zA-Z0-9]', '_', group)
                    pass


    def miR_exp(self):
        """miRNA expression
        unitas.miR-table_Human.txt
        """
        f = os.path.join(self.path_out, 'unitas.miR-table_Human.txt')

        if os.path.exists(f):
            df = pd.read_csv(f, '\t', header=0, usecols=['miR-name', 'total_reads'],
                index_col=False)
            df['miR-name'] = df['miR-name'].apply(lambda x: x.split('(')[0])
            df2 = df.groupby('miR-name').sum()
            return df2
        else:
            return None


    def tRF_exp(self):
        """tRF expression
        unitas.tRF-table.simplified.txt
        """
        f = os.path.join(self.path_out, 'unitas.tRF-table.simplified.txt')

        if os.path.exists(f):
            header = ['name', 'count', 'count_pre']
            dtype = {'name': 'str', 'count': 'float', 'count_pre': 'float'}
            df = pd.read_csv(f, sep='\t', header=None, skiprows=1,
                names=header, dtype=dtype, index_col=[0])
            return df
        else:
            return None


def main():
    args = get_args()
    assert is_path(args.o)

    ## 1. summary 
    s = []
    for d in args.i:
        # filename #
        n = d.split('_')[2]
        n = os.path.splitext(n)[0]
        a = Unitas(d).summary()
        a.index.names = ['sample'] # rename index
        a.index = [n] # rename index value
        s.append(a)
    df = pd.concat(s, axis=0, sort=True)
    # convert NaN to 0
    df = df.fillna(0)
    df2 = df.T

    summary = os.path.join(args.o, 'summary.txt')
    df2.to_csv(summary, '\t', index=True)
    # print(df2)
    

    ## 2. miRNA exp
    s = []
    for d in args.i:
        # filename #
        n = d.split('_')[2]
        n = os.path.splitext(n)[0]
        a = Unitas(d).miR_exp()
        a.columns = [n]
        s.append(a)
    df = pd.concat(s, axis=1, sort=False)
    df = df.fillna(0)

    ## save
    mir = os.path.join(args.o, 'miRNA_expression.txt')
    df.to_csv(mir, '\t', index=True)
    # print(df.head())


    ## 3. tRF exp
    s = []
    for d in args.i:
        # filename #
        n = d.split('_')[2]
        n = os.path.splitext(n)[0]
        a = Unitas(d).tRF_exp()
        a.columns = [n, n + '_pre']
        s.append(a)

    df = pd.concat(s, axis=1, sort=True)
    df = df.fillna(0)

    ## save
    trf = os.path.join(args.o, 'tRF_expression.txt')
    df.to_csv(trf, '\t', index=True)
    # print(df.head())

    pass



if __name__ == '__main__':
    main()





# p = '/home/wangming/work/others/smallRNA/results/04.annotation/unitas/UNITAS_20-01-2019_a1.fq_#1'
# a = Unitas(p).summary()
# a = Unitas(p).hits_per_target()
# a = Unitas(p).tRF_exp()
# print(a)