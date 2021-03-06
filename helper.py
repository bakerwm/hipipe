#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
support functions for hipipe
"""

import os
import sys
import re
import gzip
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


# from goldclip.bin.bed_fixer import *
# from goldclip.configure import goldclip_home
# from goldclip.helper import *


# logging.basicConfig(format = '[%(asctime)s] %(message)s', 
#                     datefmt = '%Y-%m-%d %H:%M:%S', 
#                     level = logging.DEBUG)


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)


def supportedGenome(x):
    """Check if input genome is supported by current script
    saved in file: supported_genome.txt
    """
    scriptpath = os.path.realpath(__file__)
    s = os.path.join(os.path.dirname(scriptpath), 'supported_genomes.txt')
    hit = None
    with open(s, 'rt') as fi:
        for line in fi:
            if line.strip().startswith('#') and not line.strip():
                continue
            g = line.strip()
            if x == g or x.lower() == g.lower():
                hit = g

    if not hit:
        log.warning('genomes not supported: {}'.format(x))

    return hit


def eprint(*args, **kwargs):
    """Print to stderr"""
    print(*args, file=sys.stderr, **kwargs)


def gzip_file(x, decompress=False, rm=True):
    """Compress file using gzip module
    remove original file or not
    """
    assert isinstance(x, str)

    tag = False
    if decompress:
        ## decompress
        x_ungz = os.path.splitext(x)[0] # remove extensiion
        if is_gz(x):
            if os.path.exists(x_ungz):
                log.info('target file exists, skip ungzip - %s' % x)
            else:
                with gzip.open(x, 'rb') as fi, open(x_ungz, 'wb') as fo:
                    fo.writelines(fi)

                if rm:
                    os.remove(x)
            tag = x_ungz
        else:
            log.info('expect a gzip file, skip unzip- %s' % x)
    else:
        ## compress
        xgz = x + '.gz'
        if os.path.exists(xgz):
            log.info('file exists , skip gzip - %s' % xgz)
            tag = xgz
        elif is_gz(x):
            log.info('file is gzipped, skip gzip - %s' % x)
            tag = x
        else:
            with open(x, 'rb') as fi, gzip.open(xgz, 'wb') as fo:
                fo.writelines(fi)

            if rm:
                os.remove(x)
            tag = xgz
    return tag


def gzip_file2(x, decompress=False, rm=True):
    """Compress file using gzip command
    """
    assert isinstance(x, str)

    tag = False
    if decompress:
        ## decompress
        x_ungz = os.path.splitext(x)[0] # remove extensiion
        if is_gz(x):
            if os.path.exists(x_ungz):
                log.info('target file exists, skip ungzip - %s' % x)
            else:
                cmd = 'gunzip {}'.format(x)
                run_shell_cmd(cmd)
                # with gzip.open(x, 'rb') as fi, open(x_ungz, 'wb') as fo:
                #     fo.writelines(fi)

                # if rm:
                #     os.remove(x)
            tag = x_ungz
        else:
            log.info('expect a gzip file, skip unzip- %s' % x)
    else:
        ## compress
        xgz = x + '.gz'
        if os.path.exists(xgz):
            log.info('file exists , skip gzip - %s' % xgz)
            tag = xgz
        elif is_gz(x):
            log.info('file is gzipped, skip gzip - %s' % x)
            tag = x
        else:
            cmd = 'gzip {}'.format(x)
            run_shell_cmd(cmd)
            # with open(x, 'rb') as fi, gzip.open(xgz, 'wb') as fo:
            #     fo.writelines(fi)

            # if rm:
            #     os.remove(x)
            tag = xgz
    return tag


def run_shell_cmd(cmd):
    """This command is from 'ENCODE-DCC/atac-seq-pipeline'
    https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_common.py
    """
    p = subprocess.Popen(['/bin/bash','-o','pipefail'], # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid) # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(pid, pgid, rc,
        stderr.strip(), stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')


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
    """Check if dict and x are consitent
    d is dict
    x is pickle file
    """
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
        log.error('illegal x= argument: %s' % x)

    return flag


def args_logger(d, x, overwrite=False):
    """Format dict, save to file
        key: value
    """
    assert isinstance(d, dict)
    n = ['%20s : %-40s' % (k, d[k]) for k in sorted(list(d.keys()))]
    if os.path.exists(x) and overwrite is False:
        return True
    else:
        with open(x, 'wt') as fo:
            fo.write('\n'.join(n) + '\n')
        return '\n'.join(n)


def is_gz(filepath):
    if os.path.exists(filepath):
        with open(filepath, 'rb') as test_f:
            return binascii.hexlify(test_f.read(2)) == b'1f8b'
    else:
        if filepath.endswith('.gz'):
            return True
        else:
            return False


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
                log.error('failed to create directories: %s' % path)
        else:
            return False


def seq_type(fn, top_n = 1000):
    """
    Check the top 1000 rows of fn
    identify @ for fastq, > for fasta, * unknown
    """
    assert isinstance(fn, str)
    if os.path.exists(fn):
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
            fx_type = 'fastq'
        elif tag ==  {'>'}:
            fx_type = 'fasta'
        else:
            fx_type = None
    else:
        fn_name = os.path.basename(fn).lower()
        if '.fa' in fn_name or '.fasta' in fn_name:
            fx_type = 'fasta'
        elif '.fq' in fn_name or '.fastq' in fn_nmae:
            fx_type = 'fastq'
        else:
            fx_type = None

    return fx_type


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


def fx_counter(fn):
    """Count fastq and fasta file
    only support:
    sequence in one line:
    fastq: 1 record = 4-line
    fasta: 1 record = 2-line
    """
    fx_type = seq_type(fn)
    fx_lines = file_row_counter(fn)
    if is_empty_file(fn):
        return 0
    elif fx_type == 'fasta':
        return int(fx_lines / 2)
    elif fx_type == 'fastq':
        return int(fx_lines / 4)
    else:
        return None


def str_common(x, suffix=False, longest=False):

    from difflib import SequenceMatcher

    def common_start(sa, sb):
        """ returns the longest common substring from the beginning of sa and sb """
        def _iter():
            for a, b in zip(sa, sb):
                if a == b:
                    yield a
                else:
                    return

        return ''.join(_iter())

    # see: https://stackoverflow.com/a/39404777/2530783
    def common_longest(s1, s2):
        match = SequenceMatcher(None, s1, s2).find_longest_match(0, len(s1), 0, len(s2))
        return s1[match.a: match.a + match.size]


    def common_multi(s, longest=False):
        if isinstance(s, str):
            return s
        elif len(s) == 2:
            if longest:
                return common_longest(s[0], s[1])
            else:
                return common_start(s[0], s[1])
        elif len(s) > 2:
            s1 = s.pop()
            s2 = s.pop()
            if longest:
                sx = common_longest(s1, s2)
            else:
                sx = common_start(s1, s2)
            s.append(sx)
            # return
            return common_multi(s)


    if suffix:
        # print('suffix')
        x = [f[::-1] for f in x]
        x_common = common_multi(x, longest)
        x_common = x_common[::-1]
    else:
        # print('prefix')
        x_common = common_multi(x, longest)

    return(x_common)


# def str_common(strList, suffix = False):
#     # extract longest prefix/suffix from list of strings
#     # default: prefix
#     # sort strings by len
#     def iterStop(exp):
#         if exp is False:
#             raise StopIteration
#         else:
#             return True    

#     def commonPrefix(s1, s2):
#         sys.exit('{} : {}'.format(s1, s2))
#         # prefix
#         return ''.join(list(val for i, val in enumerate(s1) 
#                        if iterStop(s2[i] is val)))

#     def fact(l):
#         if len(l) ==  1:
#             return l[0]
#         else:
#             la = l[0:2]
#             lb = l[2:]
#             s = commonPrefix(la[0], la[1])
#             lb.insert(0, s)
#             return fact(lb)

#     ## empty or single item 
#     if len(strList) ==  0:
#         return ''
#     elif len(strList) ==  1:
#         return strList[0]
#     else:
#         ## save a copy of list
#         L2 = sorted(strList, key = len)
#         c = fact(L2)
    
#     ## suffix, reverse strings
#     if suffix is True:
#         L2 = [i[::-1] for i in L2]
#         c = fact(L2)
#         c = c[::-1]

#     return c # string 0-index


def file_prefix(fn, with_path = False):
    """
    extract the prefix of a file
    remove extensions
    .gz, .fq.gz
    """
    assert isinstance(fn, str)
    p1 = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    if px.endswith('gz') or px.endswith('.bz2'):
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


def rm_file(x):
    """remove files"""
    if isinstance(x, str):
        if os.path.exists(x):
            log.info('removing file: %s' % x)
            os.remove(x)
        else:
            log.info('file does not exists: %s' % x)
    elif isinstance(x, list):
        tmp = [rm_file(i) for i in x]


def is_empty_file(x):
    """Test if file is empty
    Return True fi the uncompressed data in x have zero length
    or x itself has zero length
    """
    with xopen(x, 'rb') as fi:
        data = fi.read(1)
    return len(data) == 0


################################################################################
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
            log.error('unknown version of python: ' + sys.version)        
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
            log.error('virtualenv not exists - ' + venv)
    else: # exit env
        if is_venv() and venv_in ==  sys.prefix:
            _venv_into(venv_in, out = True)
        else:
            return ('not in venv: ' + venv_in)


################################################################################
## config
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


def bed_parser(fn, usecols = None):
    """
    read BED file as pandas DataFrame
    select specific columns, default all, (None)
    require at least 6 columns
    """
    if not pathlib.Path(fn).is_file() or os.path.getsize(fn) ==  0:
        df = pd.DataFrame(columns = ['chr', 'start', 'end', 'name', 'score', 
                                     'strand'])
        log.warning('empty bed file: %s' % fn)
        return df
    else:
        df = pd.read_table(fn, '\t', usecols = usecols, header = None,
            dtype = {'0': np.str, '1': np.int64, '2': np.int64, '3': np.str, \
                '4': np.int64, '5': np.str})
        df = df.rename(index = str, columns = {0: 'chr', 1: 'start', 2: 'end', \
                3: 'name', 4: 'score', 5: 'strand'})
        return bed_fixer(df)


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
       log.info('bigWig file exists, %s' % bw_out) 
    else:
        with open(bw_log, 'wt') as fo:
            subprocess.run(shlex.split(c), stdout=fo, stderr=fo)
        if not os.path.exists(bw_out):
            raise ValueError('failed to create bigWig file, %s' % bw_out)


def bam2bigwig(bam, genome, path_out, strandness=0, binsize=1, overwrite=False, **kwargs):
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
               'hg38': 2913022398,
               'GRCh38': 2913022398}
    gsize = effsize[genome]

    # create bam index
    BAM(bam).index()

    # prefix = os.path.basename(os.path.splitext(bam)[0])
    prefix = file_prefix(bam)[0]
    bw_log = os.path.join(path_out, prefix + '.deeptools.log')
    log.info('create bigWig for: %s' % prefix)
    if strandness > 0:
        # strandness
        bw_fwd = os.path.join(path_out, prefix + '.fwd.bigWig')
        bw_rev = os.path.join(path_out, prefix + '.rev.bigWig')
        # print(bw_fwd)
        # if strandness == 2:
        #     bw_fwd, bw_rev = [bw_rev, bw_fwd]
        
        # file existence
        if os.path.exists(bw_fwd) and os.path.exists(bw_rev) and overwrite is False:
            log.info('file exists : %s' % prefix)
        else:
            # attention; bamCoverage using dUTP-based library
            # reverse, forward
            c1 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand forward \
                  --effectiveGenomeSize {}'.format(bam, bw_fwd, binsize, gsize)
            c2 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand reverse \
                  --effectiveGenomeSize {}'.format(bam, bw_rev, binsize, gsize)
            if os.path.exists(bw_fwd) and os.path.exists(bw_rev) and not overwrite:
                log.info('file exists, bigWig skipped ...')
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
            log.info('bigWig file exists, skipping: %s' % prefix)
        else:
            c3 = 'bamCoverage -b {} -o {} --binSize {} \
                  --effectiveGenomeSize {}'.format(bam, bw, binsize, gsize)
            if os.path.exists(bw) and not overwrite:
                log.info('file exists, bigWig skipped ...')
            else:
                with open(bw_log, 'wt') as fo:
                    subprocess.run(shlex.split(c3), stdout=fo, stderr=fo)
            if not os.path.exists(bw):
                raise ValueError('output file is missing, check log file: %s' % bw_log)


# def bam2bw(bam, genome, path_out, strandness=True, binsize=1, overwrite=False):
#     """
#     Convert bam to bigWig using deeptools
#     https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
#     history:
#     1. Mappable sequence of a genome, see Table 1 in 
#        url: https://www.nature.com/articles/nbt.1518.pdf
#     2. effective genome size:
#         - non-N bases
#         - regions (of some size) uniquely mappable
#     3. UCSC
#     http://genomewiki.ucsc.edu/index.php/Hg19_100way_Genome_size_statistics
#     http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics
#     """
#     assert is_path(path_out)
#     effsize = {'dm3': 162367812,
#                'dm6': 142573017,
#                'mm9': 2620345972,
#                'mm10': 2652783500,
#                'hg19': 2451960000,
#                'hg38': 2913022398,}
#     gsize = effsize[genome]
#     # prefix = os.path.basename(os.path.splitext(bam)[0])
#     prefix = file_prefix(bam)[0]
#     bw_log = os.path.join(path_out, prefix + '.deeptools.log')
#     if strandness:
#         bw_fwd = os.path.join(path_out, prefix + '.fwd.bigWig')
#         bw_rev = os.path.join(path_out, prefix + '.rev.bigWig')
#         c1 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand forward \
#               --normalizeTo1x {}'.format(bam, bw_fwd, binsize, gsize)
#         c2 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand reverse \
#               --normalizeTo1x {}'.format(bam, bw_rev, binsize, gsize)
#         if os.path.exists(bw_fwd) and os.path.exists(bw_rev) and not overwrite:
#             log.info('file exists, bigWig skipped ...')
#         else:
#             with open(bw_log, 'wt') as fo:
#                 subprocess.run(shlex.split(c1), stdout=fo, stderr=fo)
#             with open(bw_log, 'wa') as fo:
#                 subprocess.run(shlex.split(c2), stdout=fo, stderr=fo)
#     else:
#         bw = os.path.join(path_out, prefix + '.bigWig')
#         c3 = 'bamCoverage -b {} -o {} --binSize {} \
#               --normalizeTo1x {}'.format(bam, bw, binsize, gsize)
#         if os.path.exists(bw) and not overwrite:
#             log.info('file exists, bigWig skipped ...')
#         else:
#             with open(bw_log, 'wt') as fo:
#                 subprocess.run(shlex.split(c3), stdout=fo, stderr=fo)


def bam2bw(bam, path_out, scale=1, **args):
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

    --scaleFactor
    --filterRNAstrand
    --samFlagExclude
    --samFlagInclude
    --genome


    """
    assert os.path.exists(bam)
    assert is_path(path_out)
    assert isinstance(scale, float)

    effsize = {'dm3': 162367812,
               'dm6': 142573017,
               'mm9': 2620345972,
               'mm10': 2652783500,
               'hg19': 2451960000,
               'hg38': 2913022398,
               'GRCh38': 2913022398}
    gsize = effsize.get(args['genome'], None)

    ## command
    bamcoverage_exe = which('bamCoverage')
    if bamcoverage_exe is None:
        raise ValueError('%10s | program not found: bamCoverage' % 'failed')

    # create bam index
    BAM(bam).index()

    prefix = file_prefix(bam)[0]
    bw_log = os.path.join(path_out, prefix + '.deeptools.log')
    log.info('create bigWig for: %s' % prefix)
    
    ## genome
    if gsize is None:
        ## no model organism
        opt1 = ''
    else:
        opt1 = '--effectiveGenomeSize ' + str(gsize)

    ## filt RNA strand
    if args['filterRNAstrand'] is None:
        bw_file = os.path.join(path_out, prefix + '.bigWig')
        opt2 = ''
    else:
        if args['filterRNAstrand'] == 'forward':
            if args['library_type'] == 2:
                bw_file = os.path.join(path_out, prefix + '.rev.bigWig')
            else:
                bw_file = os.path.join(path_out, prefix + '.fwd.bigWig')
        elif args['filterRNAstrand'] == 'reverse':
            if args['library_type'] == 2:
                bw_file = os.path.join(path_out, prefix + '.fwd.bigWig')
            else:
                bw_file = os.path.join(path_out, prefix + '.rev.bigWig')
        else:
            log.error('argument --filterRNAstrand: invalid choice: (\'forward\', \'reverse\')')
        opt2 = '--filterRNAstrand ' + str(args['filterRNAstrand'])

    ## opt
    opt3 = '--scaleFactor %s -b %s -o %s --binSize %s -p %s' % (scale, bam, bw_file, args['binsize'], args['p'])
  
    ## opt
    if args['normalizeUsing']:
        opt3 += ' --normalizeUsing %s' % (args['normalizeUsing'])

    ## cmd
    cmd = ' '.join([bamcoverage_exe, opt1, opt2, opt3])
    
    ## run
    if os.path.exists(bw_file) and not args['overwrite']:
        log.info('file exists, skipped: %s' % bam)
    else:
        bw_log = os.path.splitext(bw_file)[0] + '.log'
        with open(bw_log, 'wt') as fo:
            p = subprocess.run(shlex.split(cmd), stdout=fo, stderr=fo)

    ## check rresults
    if not os.path.exists(bw_file):
        raise Exception('bamCoverage failed, output not found: %s' % bw_file)

    ## return
    return bw_file


################################################################################
def download(url, file_name):
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = get(url)
        # write to file
        file.write(response.content)


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
        self.repeat_masked_genome = repeat_masked_genome
        self.kwargs = kwargs

        if genome_path is None:
            genome_path = os.path.join(str(pathlib.Path.home()), 'data', 'genome')
        self.genome_path = genome_path

        if not supportedGenome(genome):
            log.error('genome not supported: {}'.foramt(genome))


    def get_fa(self):
        """Get the fasta file of specific genome
        {genome}/bigZips/{genome}.fa
        also check ".gz" file
        """
        fa = os.path.join(self.genome_path, self.genome, 'bigZips', self.genome + '.fa')
        if not os.path.exists(fa):
            # gencode version
            fa = os.path.join(self.genome_path, self.genome, 'fasta', self.genome + '.fa')

        fa_gz = fa + '.gz'
        if not os.path.exists(fa):
            if os.path.exists(fa_gz):
                log.error('require to unzip the fasta file: %s' % fa_gz)
            else:
                log.error('fasta file not detected: %s' % fa)
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
        fa = self.get_fa()
        fa_size = fa + '.chrom.sizes'

        if not os.path.exists(fa_size):
            # log.info('Downloading chrom.sizes from UCSC')
            # url = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes' % self.genome
            # download(url, fa_size)
            log.warning('file not exists, run samtools faidx to generate it')
            pysam.faidx(fa) # create *.fa.fai
            os.rename(fa + '.fai', fa_size)

        return fa_size


    def index_validator(self, index, aligner):
        """Validate the index for aligner
        search the aligner in $PATH
        1. bowtie: bowtie-inspect -s <index>
        2. bowtie2: bowtie2-inspect -s <index>
        3. STAR: <index>/Genome, file exists
        4. ...
        """
        # check command
        aligner_exe = which(aligner)

        if not aligner_exe:
            raise Exception('aligner not detected in $PATH: %s' % aligner)
        # else:
        #     aligner = aligner_exe

        # if index is None:
        #     return None

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


    def index_finder(self, aligner='bowtie', rRNA=False,
        repeat_masked_genome=False):
        """Find the index for aligner: STAR, bowtie, bowtie2
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
        # rRNA
        if rRNA:
            # choose the first one
            rRNA_list = ['MT_trRNA', 'rRNA']
            index_rRNA = [os.path.join(self.genome_path, self.genome, aligner + '_index', i) for i in rRNA_list]
            index_rRNA = [i for i in index_rRNA if self.index_validator(i, aligner)]
            if len(index_rRNA) >= 1:
                index = index_rRNA[0]
            else:
                index = None
        # genome
        else:
            if repeat_masked_genome:
                tag = 'genome_rm'
            else:
                tag = 'genome'
            index = os.path.join(self.genome_path, self.genome, aligner + '_index', tag)

        # validate
        if not self.index_validator(index, aligner):
            index = None
        return index


    def bowtie_index(self, rRNA=False):
        """Return the bowtie index for the genome
        optional, return rRNA index
        """
        index = self.index_finder(aligner='bowtie', rRNA=rRNA,
            repeat_masked_genome=self.repeat_masked_genome)
        return index


    def bowtie2_index(self, rRNA=False):
        """Return the bowtie2 index for the genome
        optional, return rRNA index
        """
        index = self.index_finder(aligner='bowtie2', rRNA=rRNA,
            repeat_masked_genome=self.repeat_masked_genome)
        return index


    def star_index(self, rRNA=False):
        """Return the STAR index for the genome
        optional, return rRNA index
        """
        index = self.index_finder(aligner='STAR', rRNA=rRNA,
            repeat_masked_genome=self.repeat_masked_genome)
        return index


    def hisat2_index(self, rRNA=False):
        """Return the Hisat2 index for the genome
        optional, return rRNA index
        """
        index = self.index_finder(aligner='hisat2', rRNA=rRNA,
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


    def gene_bed(self, version='refseq', rmsk=False):
        """Return the gene annotation in BED format
        support UCSC, ensembl, gencode
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


    def gene_gtf(self, version='refseq'):
        """Return the gene annotation in GTF format
        support refseq, ensembl, gencode
        """
        version = version.lower() #


        gtf = os.path.join(
            self.genome_path, 
            self.genome, 
            'annotation_and_repeats',
            self.genome + '.' + version + '.gtf')
        # print('AAAA1')
        # print(gtf)

        if not os.path.exists(gtf):
            gtf = os.path.join(
            self.genome_path, 
            self.genome, 
            'gtf',
            self.genome + '.' + version + '.gtf')
        # print('AAAA2')
        # print(gtf)

        if not os.path.exists(gtf):
            gtf = None
        # print('AAAA3')
        # print(gtf)

        # if version.lower() == 'ucsc':
        #     suffix = '.refseq.gtf'
        # elif version.lower() == 'ensembl':
        #     suffix = '.ensembl.gtf'
        # else:
        #     suffix = '.gtf'
        # g = os.path.join(self.genome_path, self.genome, 'annotation_and_repeats',
        #     self.genome + suffix)
        # if not os.path.exists(g):
        #         g = None
        return gtf


    def te_gtf(self, format='gtf'):
        """Return TE annotation of the genome
        or return TE consensus sequence for the genome (dm3)
        """
        # only dm3 supported
        te_gtf = os.path.join(self.genome_path, self.genome, 
            self.genome + '_transposon', 
            self.genome + '_transposon.gtf')
        if not os.path.exists(te_gtf):
            te_gtf = None

        return te_gtf


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

        if os.path.exists(bai):
            return True
        else:
            return False


    def sort(self):
        """Sort bam fle"""
        pass


    def merge(self):
        """Merge multiple BAM files"""
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
                log.error('empty bed file: %s' % bed)
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
        

