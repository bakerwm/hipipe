#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
support functions for goldclip
"""

import os
import sys
import re
import datetime
import json
import glob
import argparse
import shlex
import subprocess
import pathlib
import warnings
import logging
import numpy as np
import pandas as pd
import pysam
import pybedtools
import binascii
from goldclip.bin.bed_fixer import *
from goldclip.configure import goldclip_home
from goldclip.helper import *

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)



##-------------------------------------------##
## formatter
def nested_dict_values(d):
    """
    get all values from nested dict
    """
    for v in d.values():
        if isinstance(v, dict):
            yield from nested_dict_values(v)
        else:
            yield v
            

def is_gz(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def xopen(fn, mode='r', bgzip=False):
    """
    Read / Write regular and gzip file, also support stdin
    """
    assert isinstance(fn, str)
    if fn == '-':
        return sys.stdin if 'r' in mode else sys.stdout
    if fn.endswith('.gz') and mode.startswith('w') or is_gz(fn):
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)


def get_time():
    """
    get current time in this format:
    2006-01-02 14:23:35
    """
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def is_path(path, create = True):
    """
    Check path, whether a directory or not
    if not, create it
    """
    assert isinstance(path, str)
    if os.path.exists(path):
        return True
    else:
        if create:
            try:
                os.makedirs(path)
                return True
            except IOError:
                logging.error('failed to create directories: %s' % path)
        else:
            return False


def seq_type(fn, top_n = 1000):
    """
    Check the top 1000 rows of fn
    identify @ for fastq, > for fasta, * unknown
    """
    assert isinstance(fn, str)
    tag = set()
    with xopen(fn, 'rt') as fi:
        for i, line in enumerate(fi):
            if i > top_n:
                break
            elif i % 4 == 0:
                b = line[0] # the first base
                if b.lower() in 'acgtn':
                    continue
                else:
                    tag.add(line[0])
            else:
                continue
    if tag ==  {'@'}:
        return 'fastq'
    elif tag ==  {'>'}:
        return 'fasta'
    else:
        return None


def is_fastq(fn):
    if seq_type(fn) == 'fastq':
        return True
    else:
        return False


def is_fasta(fn):
    if seq_type(fn) == 'fasta':
        return True
    else:
        return False


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


def str_common(strList, suffix = False):
    # extract longest prefix/suffix from list of strings
    # default: prefix
    # sort strings by len
    def iterStop(exp):
        if exp is False:
            raise StopIteration
        else:
            return True    

    def commonPrefix(s1, s2):
        # prefix
        return ''.join(list(val for i, val in enumerate(s1) 
                       if iterStop(s2[i] is val)))

    def fact(l):
        if len(l) ==  1:
            return l[0]
        else:
            la = l[0:2]
            lb = l[2:]
            s = commonPrefix(la[0], la[1])
            lb.insert(0, s)
            return fact(lb)

    ## empty or single item 
    if len(strList) ==  0:
        return ''
    elif len(strList) ==  1:
        return strList[0]
    else:
        ## save a copy of list
        L2 = sorted(strList, key = len)
        c = fact(L2)
    
    ## suffix, reverse strings
    if suffix is True:
        L2 = [i[::-1] for i in L2]
        c = fact(L2)
        c = c[::-1]

    return c # string 0-index


def file_prefix(fn, with_path = False):
    """
    extract the prefix of a file
    remove extensions
    .gz, .fq.gz
    """
    assert isinstance(fn, str)
    p1 = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    if px.endswith('gz') or px.endswith('.bz'):
        px = os.path.splitext(p1)[1] + px
        p1 = os.path.splitext(p1)[0]
    if not with_path:
        p1 = os.path.basename(p1)
    return [p1, px]


def rm_suffix1(fn):
    """
    simplify the name of bam files
    from: {name}.not_{}.not_{}.....map_{}
    to: {name}
    """
    if '.' in fn:
        p = os.path.splitext(fn)[0]
        px = os.path.splitext(fn)[1]
        if px.startswith('.not_') or px.startswith('.map_'):
            return rm_suffix1(p)
        else:
            return fn
    else:
        return fn


def filename_shorter(fn, with_path=False):
    """
    input: name1.not_spikein.not_mtrRNA.map_genome.bam 
           name2.not_spikein.not_mtrRNA.map_genome.bam
    output: name1.bam
            name2.bam
    """
    p1 = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    p2 = rm_suffix1(p1)
    if not with_path:
        p2 = os.path.basename(p2)
    return p2 + px

##--------------------------------------------##
## virtualenv 
def is_venv():
    """
    determine if inside virtualenv
    """
    return (hasattr(sys, 'real_prefix') or 
        (hasattr(sys, 'base_prefix') and sys.base_prefix !=  sys.prefix))

# current, run
def _base_prefix():
    if hasattr(sys, 'real_prefix'):
        return sys.real_prefix
    elif hasattr(sys, 'base_prefix'):
        return sys.base_prefix
    else:
        return None


# enter env
def _venv_into(venv, out = False):
    venv_bin = os.path.join(venv, 'bin', 'activate_this.py')
    if not out:
        if sys.version_info[0:2] ==  (2, 7):
            execfile(venv_bin, dict(__file__ = venv_bin)) # python2        
        elif sys.version_info[0:1] >=  (3, ):
            exec(open(venv_bin).read(), {}, dict(__file__ = venv_bin)) #python3
        else:
            logging.error('unknown version of python: ' + sys.version)        
    else:
        pass
        #subprocess.run(['deactivate'])


def venv_checker(venv = '~/envs/py27', into_venv = True):
    """
    check virtualenv 
    if not, go into
    into_venv, into env or out env
    """
    venv_in = os.path.expanduser(venv)
    
    if into_venv: # go into env
        if is_venv() and venv_in ==  _base_prefix():
            return ('already in venv: ' + venv_in)
        elif os.path.exists(venv_in):
            _venv_into(venv_in)
        else:
            logging.error('virtualenv not exists - ' + venv)
    else: # exit env
        if is_venv() and venv_in ==  sys.prefix:
            _venv_into(venv_in, out = True)
        else:
            return ('not in venv: ' + venv_in)



##--------------------------------------------##
## config
def bam2bw(bam, genome, path_out, strandness=True, binsize=1, overwrite=False):
    """
    Convert bam to bigWig using deeptools
    https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
    history:
    1. Mappable sequence of a genome, see Table 1 in 
       url: https://www.nature.com/articles/nbt.1518.pdf
    2. effective genome size:
        - non-N bases
        - regions (of some size) uniquely mappable
    3. UCSC
    http://genomewiki.ucsc.edu/index.php/Hg19_100way_Genome_size_statistics
    http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics
    """
    assert is_path(path_out)
    effsize = {'dm3': 162367812,
               'dm6': 142573017,
               'mm9': 2620345972,
               'mm10': 2652783500,
               'hg19': 2451960000,
               'hg38': 2913022398,}
    gsize = effsize[genome]
    # prefix = os.path.basename(os.path.splitext(bam)[0])
    prefix = file_prefix(bam)[0]
    bw_log = os.path.join(path_out, prefix + '.deeptools.log')
    if strandness:
        bw_fwd = os.path.join(path_out, prefix + '.fwd.bigWig')
        bw_rev = os.path.join(path_out, prefix + '.rev.bigWig')
        c1 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand forward \
              --normalizeTo1x {}'.format(bam, bw_fwd, binsize, gsize)
        c2 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand reverse \
              --normalizeTo1x {}'.format(bam, bw_rev, binsize, gsize)
        if os.path.exists(bw_fwd) and os.path.exists(bw_rev) and not overwrite:
            logging.info('file exists, bigWig skipped ...')
        else:
            with open(bw_log, 'wt') as fo:
                subprocess.run(shlex.split(c1), stdout=fo, stderr=fo)
            with open(bw_log, 'wa') as fo:
                subprocess.run(shlex.split(c2), stdout=fo, stderr=fo)
    else:
        bw = os.path.join(path_out, prefix + '.bigWig')
        c3 = 'bamCoverage -b {} -o {} --binSize {} \
              --normalizeTo1x {}'.format(bam, bw, binsize, gsize)
        if os.path.exists(bw) and not overwrite:
            logging.info('file exists, bigWig skipped ...')
        else:
            with open(bw_log, 'wt') as fo:
                subprocess.run(shlex.split(c3), stdout=fo, stderr=fo)



def bam_merge(bam_ins, bam_out):
    """
    merge multiple bam files
    input: list of bam files
    input: out.bam
    """
    # check input files
    bam_flag = []
    for b in bam_ins:
        if not os.path.exists(b) is True:
            bam_flag.append(b)
    if len(bam_flag) > 0:
        sys.exit('BAM files not exists:' + '\n'.join(bam_flag))
    # check output file
    if os.path.exists(bam_out) is True:
        pass
        # sys.exit('BAM exists:' + bam_out)
    else:
        # merge
        pysam.merge('-f', bam_out + '.unsorted.bam', *bam_ins) # overwrite output BAM
        pysam.sort('-o', bam_out, bam_out + '.unsorted.bam')
        pysam.index(bam_out)
        os.remove(bam_out + '.unsorted.bam')



## virtual env
class Genome_info():
    """
    including the information of genome
    index, annotation, ...
    """

    def __init__(self, genome, **kwargs):
        assert isinstance(genome, str)
        self.kwargs = kwargs
        self.kwargs['genome'] = genome
        if not 'path_data' in kwargs:
            self.kwargs['path_data'] = os.path.join(pathlib.Path.home(), 
                                                    'data', 'genome')
        

    def get_fa(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        gfa = os.path.join(path_data, genome, 'bigZips', genome + '.fa')
        assert os.path.exists(gfa)
        return gfa


    def get_fasize(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        gsize = os.path.join(path_data, genome, 'bigZips', genome + '.chrom.sizes')
        assert os.path.exists(gsize)
        return gsize


    def bowtie_index(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        return idx_picker(genome, path_data=path_data, aligner='bowtie')


    def bowtie2_index(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        return idx_picker(genome, path_data=path_data, aligner='bowtie2')


    def hisat2_index(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        return idx_picker(genome, path_data=path_data, aligner='hisat2')


    def star_index(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        return idx_picker(genome, path_data=path_data, aligner='star')


    def phylop100(self):
        """
        only support hg19
        """
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        phylop100 = os.path.join(self.kwargs['path_data'],
                            genome, 'phyloP100way', 
                            genome + '.100way.phyloP100way.bw')
        if not os.path.exists(phylop100):
            phylop100 = None
        return phylop100

        
    def gene_bed(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        gbed = os.path.join(path_data, 
                            genome,
                            'annotation_and_repeats', 
                            genome + '.refseq.bed')
        if not os.path.exists(gbed):
            gbed = None
        return gbed


    def gene_rmsk(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        grmsk= os.path.join(path_data, 
                            genome,
                            'annotation_and_repeats', 
                            genome + '.rmsk.bed')
        if not os.path.exists(grmsk):
            grmsk = None
        return grmsk



def is_idx(path, aligner='bowtie'):
    """
    check aligner index, bowtie, bowtie2, STAR
    """
    # bowtie/bowtie2
    c = [aligner + '-inspect', '-s', path]
    if aligner.lower() == 'star':
        pg = os.path.join(path, 'Genome')
        flag = True if os.path.exists(pg) else False
    else:
        p = subprocess.run(c, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
        flag = True if len(p) > 0 else False
    return flag



def idx_picker(genome, group='genome', path_data=None, aligner='bowtie'):
    """
    return the path of index
    group: genome, rRNA, tRNA, ...
    aligner: bowtie, bowie2
    #
    default: path ~/data/genome/
    """
    assert isinstance(group, str)
    # assert isinstance(genome, str)
    if genome is None:
        return None    
    if path_data is None:
        path_data = os.path.join(pathlib.Path.home(), 'data', 'genome')
    idx = os.path.join(path_data, genome, aligner + '_index', group)
    if aligner.lower() == 'star':
        idx = os.path.join(path_data, genome, 'STAR_index')
    if is_idx(idx, aligner):
        return idx
    else:
        return None
    


def idx_grouper(genome, path_data=None, aligner='bowtie'):
    """
    return a group of indexes for genome mapping
    eg: spikein, MT_trRNA, genome
    """
    group1 = ['viral', 'repeatRNA', 'retroviral', 'MT_trRNA', 'genome']
    idxes = [idx_picker(genome, g, path_data=path_data, aligner=aligner) for g in group1]
    idxes = list(filter(None.__ne__, idxes))
    return idxes




def bed_parser(fn, usecols = None):
    """
    read BED file as pandas DataFrame
    select specific columns, default all, (None)
    require at least 6 columns
    """
    if not pathlib.Path(fn).is_file() or os.path.getsize(fn) ==  0:
        df = pd.DataFrame(columns = ['chr', 'start', 'end', 'name', 'score', 
                                     'strand'])
        logging.warning('empty bed file: %s' % fn)
        return df
    else:
        df = pd.read_table(fn, '\t', usecols = usecols, header = None,
            dtype = {'0': np.str, '1': np.int64, '2': np.int64, '3': np.str, \
                '4': np.int64, '5': np.str})
        df = df.rename(index = str, columns = {0: 'chr', 1: 'start', 2: 'end', \
                3: 'name', 4: 'score', 5: 'strand'})
        return bed_fixer(df)



# def bed_filter(fn, bed_exclude, bed_out, overlap = True, save = True):
#     """
#     remove records from fn that have overlap with bed_exclude, and 
#     save to fn using pybedtools
#     overlap, True: intersect, False: not intersect
#     """
#     assert pathlib.Path(fn).is_file()
#     bed_out_path = os.path.dirname(bed_out)
#     assert is_path(bed_out_path)
#     a = pybedtools.BedTool(fn)
#     b = pybedtools.BedTool(bed_exclude)
#     if overlap is True:
#         a_and_b = a.intersect(b, wa = True, u = True) # intersect with b
#     elif overlap is False:
#         a_and_b = a.intersect(b, wa = True, v = True) # exclude b
#     else:
#         logging.error('unknown overlap: %s' % overlap)
#     if save is True:
#         a_and_b.moveto(bed_out)
#     else:
#         return a_and_b # BedTool object


# def bed_fixer(df):
#     """
#     filt BED records 
#     1. start, end both are int
#     2. start < end
#     """
#     dx = df[['start', 'end']].apply(pd.to_numeric)
#     c = ((dx['start'] >=  0) & dx['end'] >=  0) & (dx['start'] < dx['end'])
#     return df.loc[c, :]

