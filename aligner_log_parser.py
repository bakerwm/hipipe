#!/usr/bin/env python
"""
Processing STAR log file
saving stat
"""


import os
import re
import sys
import glob
import logging
import argparse
import pandas as pd
from alignment import Alignment_log

def extract_group(fn):
    """Extract file, group, aligner"""
    assert isinstance(fn, str)
    fn_name = os.path.basename(fn)
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


def map_stat(path):
    """Stat mapping reads of each file"""
    log_files = glob.glob(os.path.join(path, '*.log'))
    if len(log_files) == 0:
        logging.error('log files not found: %s' % path)
    else:
        log_files = sorted(log_files, key=len) # by length
        frames = [extract_group(f) for f in log_files]
        df = pd.concat(frames)
        return(df)


def main():
    if len(sys.argv) < 2:
        print('Usage: aligner_log_parser.py <path>')
        sys.exit()
    p = sys.argv[1]
    df = map_stat(p)
    print(df)



if __name__ == '__main__':
    main()


## EOF
