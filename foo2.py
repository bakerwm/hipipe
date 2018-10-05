#!/usr/bin/env python

import os, sys, re
import pandas as pd
from alignment import Alignment_log, Alignment_stat
from helper import str_common


## parsing mapping stat
def align_log_reader(fn):
    """Extract file, group, aligner"""
    assert isinstance(fn, str)
    fn_name = os.path.splitext(os.path.basename(fn))[0]
    group, aligner = fn.split('.')[-3:-1]
    group = re.sub('^map_', '', group)
    # stat
    d = Alignment_log(fn).stat
    df = pd.DataFrame({'total' : d['total'],
                       'unique' : d['unique'], 
                       'multiple' : d['multiple'], 
                       'unmap' : d['unmap']}, 
                       index=[fn_name])
    return df


def align_log_dirs(path):
    """Extract alignment log from one directory"""
    assert os.path.isdir(path)
    log_files = [os.path.join(path, f) for f in os.listdir(path) if f.endswith('.log')]

    if len(log_files) == 0:
        return None
    elif len(log_files) == 1:
        # map_genome
        df = align_log_reader(log_files[0])
        return df
    elif len(log_files) > 1:
        log_files = sorted(log_files, key=len)
        # map to genome
        dfg = align_log_reader(log_files[-1])
        dfg = dfg.loc[:, ['unique', 'multiple', 'unmap']]

        # map to others
        frames = []
        for log_file in log_files[:-1]:
            # exccept map_genome
            df = align_log_reader(log_file)
            log_name = df.index.tolist()[0] # first one
            group = log_name.split('.')[-2] # map_group
            group = re.sub('map_', '', group)
            df = df.assign(mapped = df.unique + df.multiple)
            df.columns.values[-1] = group # rename column
            frames.append(df.loc[:, group])
        # all rows
        dfx = pd.concat(frames, axis=1)
        dfx.index = dfg.index.values
       
        # all stat
        df = pd.concat([dfx, dfg], axis=1)
        df = df.assign(total=pd.Series(df.sum(axis=1)).values)
        return df




def main():
    p = '/home/wangming/work/yu_2018/projects/20180720_RNAseq_nxf2/demo/results/genome_mapping/treatment/'
    df = map_stat(p)
    print(df)


if __name__ == '__main__':
    main()


# EOF