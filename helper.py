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
import pickle
import warnings
import logging
import numpy as np
import pandas as pd
import pysam
import pybedtools
import binascii
from requests import get  # to make GET request
from goldclip.bin.bed_fixer import *
from goldclip.configure import goldclip_home
from goldclip.helper import *

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)



def findfiles(which, where='.'):
    """Returns list of filenames from `where` path matched by 'which'
    shell pattern. Matching is case-insensitive.
    # findfiles('*.ogg')
    """    
    # TODO: recursive param with walk() filtering
    rule = re.compile(which.translate(which), re.IGNORECASE)
    return [name for name in os.listdir(where) if rule.match(name)]



def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def args_checker(d, x, update=False):
    """Check if dict and x are consitent"""
    assert isinstance(d, dict)
    flag = None
    if os.path.exists(x):
        # read file to dict
        with open(x, 'rb') as fi:
            d_checker = pickle.load(fi)
        if d == d_checker:
            flag = True
        else:
            if update:
                with open(x, 'wb') as fo:
                    pickle.dump(d, fo, protocol=pickle.HIGHEST_PROTOCOL)
    elif isinstance(x, str):
        # save dict to new file
        with open(x, 'wb') as fo:
                    pickle.dump(d, fo, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        logging.error('illegal x= argument: %s' % x)
    return flag


def args_logger(d, x, overwrite=False):
    """Format dict, save to file
        key: value
    """
    assert isinstance(d, dict)
    n = ['%20s : %-40s' % (k, d[k]) for k in list(d.keys())]
    if os.path.exists(x) and overwrite is False:
        return True
    else:
        with open(x, 'wt') as fo:
            fo.write('\n'.join(n) + '\n')
        return '\n'.join(n)


##-------------------------------------------##
## formatter
# def nested_dict_values(d):
#     """
#     get all values from nested dict
#     """
#     for v in d.values():
#         if isinstance(v, dict):
#             yield from nested_dict_values(v)
#         else:
#             yield v
            

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




class Genome_info(object):
    """List related information of specific genome
    1. get_fa(), genome fasta
    2. get_fasize(), genome fasta size
    3. bowtie_index(), bowtie index, optional, rRNA=True
    4. bowtie2_index(), bowtie2 index, optional, rRNA=True
    5. star_index(), STAR index, optional, rRNA=True
    6. gene_bed(), 
    7. gene_rmsk(), 
    8. gene_gtf(), optional, version='ucsc|ensembl|ncbi'
    9. te_gtf(), optional, version='ucsc'
    10. te_consensus(), optional, fruitfly()
    ...

    index, annotation, ...
    """

    def __init__(self, genome, **kwargs):
        assert isinstance(genome, str)
        self.kwargs = kwargs
        self.kwargs['genome'] = genome
        if not 'genome_path' in kwargs:
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
        idx = os.path.join(path_data, genome, 'STAR_index', group)
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

def bam2bigwig2(bam, path_out, scale=1, binsize=1, overwrite=False):
    """Convert BAM to bigWig using deeptools"""
    bamcoverage_exe = which('bamCoverage')
    BAM(bam).is_indexed() # check *.bai file
    is_path(path_out) # create path
    prefix = os.path.basename(os.path.splitext(bam)[0])
    bw_out = os.path.join(path_out, prefix + '.bigWig')
    bw_log = os.path.join(path_out, prefix + '.log')
    c = '%s --bam %s -o %s --scaleFactor %s --binSize %s' % (bamcoverage_exe, 
        bam, bw_out, scale, binsize)
    if os.path.exists(bw_out) and overwrite is False:
       logging.info('bigWig file exists, %s' % bw_out) 
    else:
        with open(bw_log, 'wt') as fo:
            subprocess.run(shlex.split(c), stdout=fo, stderr=fo)
        if not os.path.exists(bw_out):
            raise ValueError('failed to create bigWig file, %s' % bw_out)



## utilities
def bam2bigwig(bam, genome, path_out, strandness=0, binsize=1, overwrite=False):
    """Convert bam to bigWig using deeptools
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

    !!! strandness
    default: dUTP-based library (read2 is sense strand, read1 is anti-sense strand)
    general RNA library: (NSR, small-RNA-library), read1 is sense, read2 is antisense
    """
    assert os.path.exists(bam)
    assert isinstance(genome, str)
    assert is_path(path_out)
    assert isinstance(strandness, int)
    assert isinstance(binsize, int)
    assert isinstance(overwrite, bool)
    bamcov = which('bamCoverage')
    if bamcov is None:
        raise ValueError('%10s | program not found: bamCoverage' % 'failed')

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
    logging.info('create bigWig for: %s' % prefix)
    if strandness > 0:
        # strandness
        bw_fwd = os.path.join(path_out, prefix + '.fwd.bigWig')
        bw_rev = os.path.join(path_out, prefix + '.rev.bigWig')
        if strandness == 2:
            bw_fwd, bw_rev = [bw_rev, bw_fwd]
        
        # file existence
        if os.path.exists(bw_fwd) and os.path.exists(bw_rev) and overwrite is False:
            logging.info('bigWig file exists, skipped: %s' % prefix)
        else:
            # attention; bamCoverage using dUTP-based library
            # reverse, forward
            c1 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand forward \
                  --normalizeTo1x {}'.format(bam, bw_fwd, binsize, gsize)
            c2 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand reverse \
                  --normalizeTo1x {}'.format(bam, bw_rev, binsize, gsize)
            if os.path.exists(bw_fwd) and os.path.exists(bw_rev) and not overwrite:
                logging.info('file exists, bigWig skipped ...')
            else:
                with open(bw_log, 'wt') as fo:
                    p1 = subprocess.run(shlex.split(c1), stdout=fo, stderr=fo)
                    p2 = subprocess.run(shlex.split(c2), stdout=fo, stderr=fo)
            if not os.path.exists(bw_fwd) or not os.path.exists(bw_rev):
                raise ValueError('output file is missing, check log file: %s' % bw_log)
    else:
        # strandless
        bw = os.path.join(path_out, prefix + '.bigWig')
        if os.path.exists(bw) and overwrite is False:
            logging.info('bigWig file exists, skipping: %s' % prefix)
        else:
            c3 = 'bamCoverage -b {} -o {} --binSize {} \
                  --normalizeTo1x {}'.format(bam, bw, binsize, gsize)
            if os.path.exists(bw) and not overwrite:
                logging.info('file exists, bigWig skipped ...')
            else:
                with open(bw_log, 'wt') as fo:
                    subprocess.run(shlex.split(c3), stdout=fo, stderr=fo)
            if not os.path.exists(bw):
                raise ValueError('output file is missing, check log file: %s' % bw_log)



################################################################################
## split line
################################################################################

def download(url, file_name):
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = get(url)
        # write to file
        file.write(response.content)


def index_validator(index, aligner='bowtie'):
    """Check the index
    search the aligner in $PATH
    1. bowtie: bowtie-inspect -s <index>
    2. bowtie2: bowtie2-inspect -s <index>
    3. STAR: <index>/Genome, file exists
    4. ...
    """
    # check command
    aligner_exe = which(aligner)
    if not aligner_exe:
        raise ValueError('aligner not detected in $PATH: %s' % aligner)
    # else:
    #     aligner = aligner_exe

    flag = False
    if aligner.lower().startswith('bowtie'):
        # bowtie, bowtie2
        c = [aligner_exe + '-inspect', '-s', index]
        p = subprocess.run(c, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if len(p.stdout) > 0:
            flag = True
    elif aligner.lower() == 'star':
        p = os.path.join(index, 'Genome')
        if os.path.exists(p):
            flag = True
    else:
        raise ValueError('unknown aligner: %s' % aligner)

    return flag


def index_finder(genome, aligner='bowtie', rRNA=False, genome_path=None,
    repeat_masked_genome=False):
    """Find aligner index
    if rRNA is True, return the rRNA index only
    if return None, the index not found, or file not exists
    check if repeat masked genome required (for non-TE mapping)
    
    structure of genome_path:
    default: {HOME}/data/genome/{genome_version}/{aligner}/

    /genome_path/
        |- genome
        |- rRNA
        |- MT_trRNA
        |- 

    """
    # determine the path of genome data
    if genome_path is None:
        genome_path = os.path.join(pathlib.Path.home(), 'data', 'genome')

    # rRNA
    if rRNA:
        # choose the first one
        rRNA_list = ['MT_trRNA', 'rRNA']
        idx_rRNA = [os.path.join(genome_path, genome, aligner + '_index', i) for i in rRNA_list]
        idx_rRNA = [i for i in idx_rRNA if index_validator(i, aligner)]
        if len(idx_rRNA) >= 1:
            idx = idx_rRNA[0]
        else:
            idx = None
    # genome
    else:
        if repeat_masked_genome:
            tag = 'genome_rm'
        else:
            tag = 'genome'
        idx = os.path.join(genome_path, genome, aligner + '_index', tag)

    # validate
    if not index_validator(idx, aligner):
        idx = None
    return idx


class Genome(object):
    """List related information of specific genome
    1. get_fa(), genome fasta
    2. get_fasize(), genome fasta size
    3. bowtie_index(), bowtie index, optional, rRNA=True
    4. bowtie2_index(), bowtie2 index, optional, rRNA=True
    5. star_index(), STAR index, optional, rRNA=True
    6. gene_bed(), 
    7. gene_rmsk(), 
    8. gene_gtf(), optional, version='ucsc|ensembl|ncbi'
    9. te_gtf(), optional, version='ucsc'
    10. te_consensus(), optional, fruitfly()
    ...

    directory structure of genome should be like this:
    /path-to-data/{genome}/
        |- bigZips  # genome fasta, fasize, chromosome 
        |- annotation_and_repeats  # gtf, bed, rRNA, tRNA, annotation
        |- bowtie_index
        |- bowtie2_index
        |- STAR_index
        |- hisat2_index
        |- phylop100
        |- ...

    default: $HOME/data/genome/{genome}

    """

    def __init__(self, genome, genome_path=None, repeat_masked_genome=False, **kwargs):
        assert isinstance(genome, str)
        self.genome = genome
        if not genome_path:
            genome_path = os.path.join(pathlib.Path.home(), 'data', 'genome')
        self.genome_path = genome_path
        self.repeat_masked_genome = repeat_masked_genome
        self.kwargs = kwargs


    def get_fa(self):
        """Get the fasta file of specific genome
        {genome}/bigZips/{genome}.fa
        also check ".gz" file
        """
        fa = os.path.join(self.genome_path, self.genome, 'bigZips', self.genome + '.fa')
        fa_gz = fa + '.gz'
        if not os.path.exists(fa):
            if os.path.exists(fa_gz):
                logging.error('require to unzip the fasta file: %s' % fa_gz)
            else:
                logging.error('fasta file not detected: %s' % fa)
            return None
        else:
            return fa


    def get_fasize(self):
        """Get the fasta size file, chromosome size
        optional, fetch chrom size from ucsc
        http://hgdownload.cse.ucsc.edu/goldenPath/<db>/bigZips/<db>.chrom.sizes

        or using UCSC tool: fetchChromSizes
        fetchChromSizes hg39 > hg38.chrom.sizes
        """
        fa_size = os.path.join(self.genome_path, self.genome, 'bigZips', 
            self.genome + '.chrom.sizes')
        if not os.path.exists(fa_size):
            logging.info('Downloading chrom.sizes from UCSC')
            url = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes' % self.genome
            download(url, fa_size)
        return fa_size


    def bowtie_index(self, rRNA=False):
        """Return the bowtie index for the genome
        optional, return rRNA index
        """
        index = index_finder(self.genome, aligner='bowtie', rRNA=rRNA,
            genome_path=self.genome_path, 
            repeat_masked_genome=self.repeat_masked_genome)
        return index


    def bowtie2_index(self, rRNA=False):
        """Return the bowtie2 index for the genome
        optional, return rRNA index
        """
        index = index_finder(self.genome, aligner='bowtie2', rRNA=rRNA,
            genome_path=self.genome_path, 
            repeat_masked_genome=self.repeat_masked_genome)
        return index


    def star_index(self, rRNA=False):
        """Return the STAR index for the genome
        optional, return rRNA index
        """
        index = index_finder(self.genome, aligner='STAR', rRNA=rRNA,
            genome_path=self.genome_path, 
            repeat_masked_genome=self.repeat_masked_genome)
        return index


    def hisat2_index(self, rRNA=False):
        """Return the Hisat2 index for the genome
        optional, return rRNA index
        """
        index = index_finder(self.genome, aligner='hisat2', rRNA=rRNA,
            genome_path=self.genome_path, 
            repeat_masked_genome=self.repeat_masked_genome)
        return index


    def phylop100(self):
        """Return the phylop100 bigWig file of hg19, only
        for conservation analysis
        """
        p = os.path.join(self.genome_path, self.genome, 'phyloP100way',
            self.genome + '.100way.phyloP100way.bw')
        if not os.path.exists(p):
            p = None
        return p


    def gene_bed(self, version='ucsc', rmsk=False):
        """Return the gene annotation in BED format
        support UCSC, ensembl version
        """
        if rmsk:
            suffix = '.rmsk.bed'
        else:
            suffix = '.refseq.bed'
        g = os.path.join(self.genome_path, self.genome, 'annotation_and_repeats',
            self.genome + suffix)
        if not os.path.exists(g):
            g = None
        return g


    def gene_gtf(self, version='ucsc'):
        """Return the gene annotation in GTF format
        support UCSC, ensembl version
        """
        if version.lower() == 'ucsc':
            suffix = '.refseq.gtf'
        elif version.lower() == 'ensembl':
            suffix = '.ensembl.gtf'
        else:
            suffix = '.gtf'
        g = os.path.join(self.genome_path, self.genome, 'annotation_and_repeats',
            self.genome + suffix)
        if not os.path.exists(g):
                g = None
        return g


    def te(self):
        """Return TE annotation of the genome
        or return TE consensus sequence for the genome (dm3)
        """
        pass
        # additional Class


class BAM(object):
    """Operation for BAM files
    sort, index, ...
    """

    def __init__(self, fn):
        self.fn = fn


    def index(self):
        """Create index for BAM file"""
        bam = self.fn
        bai = self.fn + '.bai'
        if not os.path.exists(bai):
            pysam.index(bam)


    def sort(self):
        """Sort bam fle"""
        pass


    def count(self):
        """Count reads in BAM file"""
        return pysam.view('-c', self.fn)


    def to_bed(self, bed=None):
        """Convert BAM to bed using bedtools"""
        bam = self.fn
        if bed is None:
            bed = os.path.splitext(bam)[0] + '.bed'
        if not os.path.exists(bed):
            pybedtools.BedTool(bam).bam_to_bed().saveas(bed)
        return bed
        

    def is_indexed(self, overwrite=False):
        """Check if *.bai file exists"""
        bai = self.fn + '.bai'
        if os.path.exists(bai) and overwrite is False:
            return True
        else:
            pysam.index(self.fn)
            if os.path.exists(bai):
                return True
            else:
                return False