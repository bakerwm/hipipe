#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
annotate bed file by annotatePeaks.pl from HOMER

$ hipipe-anno2.py -i in.bed -g mm9 -o out.anno.bed --is-clipper-peak

"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


import os
import re
import pathlib
import logging
import argparse
import shlex
import subprocess
import pandas as pd
from helper import which
from bed_fixer import Bed_parser


def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog='bed_annotation',
        description='annotation for bed file',
        epilog='Example: hipipe-anno2.py -i demo.bed -g dm3 -o demo.anno')
    parser.add_argument('-i', required=True, metavar='BED', 
        help='BED or BAM file(s)')
    parser.add_argument('-g', required=True, default='dm3', 
        metavar='GENOME', choices=['dm3', 'hg19', 'hg38', 'mm10', 'mm9'],
        help='Reference genome : [dm3, hg19, hg38, mm10, mm9], default: dm3')
    parser.add_argument('-o', metavar='output', default=None,
        help='Save results to file, default, <in.bed>.anno')
    parser.add_argument('--is-clipper-peak', dest='k', action='store_true',
        help='the input is CLIPper output, to extract count in name field')
    parser.add_argument('--hub-txt', dest='u', default=None,
        help='the url of hub.txt file')
    args = parser.parse_args()
    return args



def hub_to_url(hub_txt, genome, position=None):
    """
    create url for hub.txt
    position: chr1:1-1000 convert to url encode
    'chr1%3A1-1000'
    """
    base_url = 'http://genome.ucsc.edu/cgi-bin/hgTracks'
    db = '?db=' + genome
    hubUrl = '&hubUrl=' + hub_txt
    if position is None:
        return(base_url + db + hubUrl)
    else:
        pos = '&position=' + re.sub(r':', r'%3A', position)
        return(base_url + db + hubUrl + pos)



def clipper_bed_fixer(fn):
    """Fix CLIPper output bed file
    1. extract read count from name
    """
    df = Bed_parser(fn).bed

    # extract count
    def id2count(s):
        return s.split('_')[-1]

    cnt = df['name'].apply(id2count)
    peak_len = df.apply(lambda row: row['end'] - row['start'], axis=1)
    df = df.assign(count=pd.Series(cnt).values,
        length=pd.Series(peak_len).values)
    return df



def anno_run(fn, fout, flog, genome):
    """Run annotatePeaks.pl from HOMER"""
    assert os.path.exists(fn)
    assert isinstance(fout, str)
    assert isinstance(flog, str)
    assert isinstance(genome, str)

    annotatePeaks = which('annotatePeaks.pl')
    if annotatePeaks is None:
        raise ValueError('%10s | command not found: annotatePeaks.pl' % 'failed')

    logging.info('running annotatePeaks.pl for %s' % os.path.basename(fn))
    c1 = '%s %s %s' % (annotatePeaks, fn, genome)
    with open(fout, 'wt') as ffout, open(flog, 'wt') as fflog:
        p1 = subprocess.run(shlex.split(c1), stdout=ffout, stderr=fflog)
    
    if os.path.exists(fout):
        return fout
    else:
        return None



def anno_append(fn, fanno, fout, clipper_bed=False):
    """Append annotation to bed file"""
    assert os.path.exists(fn)
    assert os.path.exists(fanno)
    assert isinstance(fout, str)

    logging.info('append annotation to bed file')

    if clipper_bed is True:
        df1 = clipper_bed_fixer(fn)
    else:
        df1 = Bed_parser(fn).bed
    df2 = pd.read_csv(fanno, "\t")

    # rename header
    #
    # 1  'PeakID (cmd=annotatePeaks.pl MIWI_piRNA.clipper.count.bed mm9)', 
    # 2  'Chr',
    # 3  'Start', 
    # 4  'End', 
    # 5  'Strand', 
    # 6  'Peak Score', 
    # 7  'Focus Ratio/Region Size',
    # 8  'Annotation', 
    # 9  'Detailed Annotation', 
    # 10 'Distance to TSS',
    # 11 'Nearest PromoterID', 
    # 12 'Entrez ID', 
    # 13 'Nearest Unigene', 
    # 14 'Nearest Refseq',
    # 15 'Nearest Ensembl', 
    # 16 'Gene Name', 
    # 17 'Gene Alias', 
    # 18 'Gene Description',
    # 19 'Gene Type'
    
    # to
    # h0 ... h18

    df2.columns = ['h' + str(i) for i in range(len(df2.columns))]

    # 1  'PeakID'
    # 16 'Gene Name'
    # 19 'Gene Type'

    df2x = df2.iloc[:, [0, 15, 18]]
    df2x.columns = ['name', 'gene', 'type']

    # merge two data.frame
    # default_kwargs = dict(sep='\t', header=False, index=False)
    # df3 = pd.merge(df1, df2x, how='left', on='name').to_csv(fout, **default_kwargs)
    df3 = pd.merge(df1, df2x, how='left', on='name')
    return df3



def hub_append(df, hub_txt, genome):
    """Add UCSC track url for each bed record
    position: chr13:1-1000
    """
    assert isinstance(df, pd.DataFrame)

    # add position column
    def get_pos(chr, start, end, expand=200):
        """return chr:start-end"""
        if abs(expand) > 0:
            start = int(start) - expand
            end = int(end) + expand
        if start < 0:
            start = 1
        return '%s:%s-%s' % (str(chr), str(start), str(end))

    pos = df.apply(lambda row: get_pos(row['chr'], row['start'], row['end']), axis=1)
    df1 = df.assign(position=pd.Series(pos).values)
    url = df1.apply(lambda row: hub_to_url(hub_txt, genome, 
        row['position']), axis=1)
    df2 = df1.assign(pos_url=pd.Series(url).values)
    return df2




def main():
    args = get_args()
    f_in = args.i
    f_out = args.o
    genome = args.g
    clipper_bed = args.k
    hub_txt = args.u
    assert os.path.exists(f_in)
    if f_out is None:
        f_out = f_in + '.annotate.bed'

    logging.info('start annotation')
    # temp file
    # f_tmp = tempfile.NamedTemporaryFile(prefix='tmp',
    #                                     suffix='.bed',
    #                                     delete=False)
    # f_tmp = f_tmp.name
    f_tmp = f_out + '.annotatePeaks.out'
    f_log = os.path.splitext(f_out)[0] + '.log'
    f_anno = anno_run(f_in, f_tmp, f_log, genome)
    if f_anno is None:
        raise ValueError('%10s | annotatePeaks.pl goes wrong' % 'failed')
    df_anno = anno_append(f_in, f_anno, f_out, clipper_bed)

    # add pos_url
    if isinstance(hub_txt, str):
        df_anno = hub_append(df_anno, hub_txt, genome)
    
    # save to file
    default_kwargs = dict(sep='\t', header=False, index=False)
    df_anno.to_csv(f_out, **default_kwargs)

    # remove temp file
    os.remove(f_tmp)

    logging.info('finish')



if __name__ == '__main__':
    main()



## EOF
