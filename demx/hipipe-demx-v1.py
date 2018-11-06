#!/usr/bin/env python
"""
demultiplex P7-index
illumina i7-index: The first 6-nt index 

# Name
@ST-E00310:586:HJH7JCCXY:7:1101:27275:1555 1:N:0:NCAACAAT
@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>

url1: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
url2: https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html
"""

__author__ = 'Ming Wang'
__email__ = 'wangm08@hotmail.com'
__date__ = '2018-05-21'
__version__ = '0.1'


import os
import sys
import argparse
import glob
import datetime
import gzip
import json
import binascii
# import shlex
# import subprocess
# import shutil
import logging
# import pandas as pd


def get_args():
    """get args"""
    parser = argparse.ArgumentParser(
        prog = 'demultiplexing', description='demultiplexing reads',
        epilog='Example: ')
    parser.add_argument('--r1', required=True, metavar='READ1', 
        type=argparse.FileType('r'),
        help='SE read or read1 of PE reads')
    parser.add_argument('--r2', required=False, default=None, 
        metavar='READ2', type=argparse.FileType('r'),
        help='read2 of PE reads')
    parser.add_argument('--sample-list', required=True, 
        help='a list of sample names and index, barcode sequences\
        <p7> <barcode> <sample>')
    parser.add_argument('--p7-index', action='store_true',
        help='demultiplexing P7-index if specified')
    parser.add_argument('--barcode', action='store_true',
        help='demultiplexing barcode if specified')
    parser.add_argument('-o', '--out', required=True, metavar='OUTDIR', 
        help='The directory to save results.')
    parser.add_argument('--bc-n-left', dest='bc_n_left', default=3, 
        metavar='N-LEFT', type=int,
        help='N nt at the left of barcode. default: 3')
    parser.add_argument('--bc-n-right', dest='bc_n_right', default=2, 
        metavar='N-RIGHT', type=int,
        help='N nt at the right of barcode. default: 2')
    parser.add_argument('--n-mismatch', dest='n_mismatch', default=0, 
        metavar='Mismatches', type=int,
        help='Number of mismatches allowed in barcode/index, default: 0')
    parser.add_argument('--cut', action='store_true',
        help='cut the barcode from read if specified')
    parser.add_argument('--bc-in-read', dest='bc_in_read', choices=[1, 2],
        default=1, type=int, 
        help='barcode in read1/2 of PE reads, default: 1')
    parser.add_argument('--bioawk', required=False, action='store_true',
        help='use bioawk to demulplex only read1.')
    args = parser.parse_args()
    return args


##----------------------------------------------------------------------------##
def fileRowCounter(fn):
    """Calculate the number of rows for a file
    # cite: @Michael Bacon on stackoverflow
    # url: https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    """
    def _make_gen(reader):
        b = reader(1024 * 1024)
        while b:
            yield b
            b = reader(1024 * 1024)

    def rawgencount(fn):
        try:
            with open(fn, 'rb') as f:
                f_gen = _make_gen(f.raw.read)
                return sum(buf.count(b'\n') for buf in f_gen)
        except IOError:
            print('failed to process file: ' + fn)

    return rawgencount(fn)


def is_gz(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def str_mismatch(a, b):
    """Calculate the mismatches between two strings
    query: a
    subject: b
    if b == null, return 0 mismatch
    """
    assert isinstance(a, str)
    assert isinstance(b, str)
    m = sum(list(map(lambda x, y : 0 if x == y else 1, a, b)))
    if b.lower() == 'null':
        m = 0
    return m


class Json_file(object):
    """Parsing Json and dict file"""
    def __init__(self, fn):
        self.fn = fn
        if isinstance(fn, Json_file):
            self.stat = fn.stat
        elif isinstance(fn, dict):
            self.stat = fn
        elif os.path.exists(fn):
            self.stat = self.json_reader()
        else:
            raise ValueError('unknown file format:')


    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.json',
                                            delete=False)
        return tmpfn.name


    def json_reader(self):
        """Load json file as dict"""
        fn = self.fn
        if os.path.isfile(fn) and os.path.getsize(fn) > 0:
            with open(fn, 'rt') as ff:
                return json.load(ff)


    def json_writer(self, to=None):
        """Write dict to file in json format"""
        fn = self.fn

        if to is None:
            to = self._tmp()

        if isinstance(fn, Json_file):
            fn = fn.fn
        elif isinstance(fn, dict):
            fn = fn
        with open(to, 'wt') as ff:
            json.dump(fn, ff, indent=4, sort_keys=True)


def bc_parser(fn, p7=True, barcode=True):
    """Parse the barcode and sample name from sample list
    return dict"""
    bc_dict = {}
    bc_dup = 0
    d1 = d2 = {}
    with open(fn, 'rt') as fi:
        for line in fi:
            p = line.strip().split('\t') # <p7> <barcode> <name>
            if not len(p) == 3:
                sys.exit('require <p7> <barcode> <name> in sample list')
            p[0] = p[0].upper()
            p[1] = p[1].upper()
            # p7, bc length
            if not p[0] == 'NULL':
                d1[p[0]] = 1
            if not p[1] == 'NULL':
                d2[p[1]] = 1
            # id to name
            if p7 is True and barcode is True:
                tag = p[0] + ',' + p[1]
            elif p7 is True:
                tag = p[0]
            elif barcode is True:
                tag = p[1]
            else:
                sys.exit('either or both p7 and barcode are required')
            name = p[2] # sample name
            if tag in bc_dict:
                bc_dup += 1
            bc_dict[tag] = name
    if bc_dup > 0:
        sys.exit('[barcode] duplicated, exiting...')
    # barcode in same length
    len1 = set(list(map(len, list(d1.keys()))))
    len2 = set(list(map(len, list(d2.keys()))))
    if len(len1) > 1 or len(len2) > 1:
        raise ValueError('p7, barcode not in same length')
    # sample name
    _name = sorted(list(bc_dict.values()))
    _name_uniq = sorted(list(set(_name)))
    if not _name == _name_uniq:
        sys.exit('[sample name] duplicated, exiting...')

    return bc_dict


# checkout barcode exists in list
def seq_validator(s, bc_dict, mm=0):
    """Check the barcode|index for each read
    criteria: mismatches <= mm
    return sample name
    """
    n = [x for x in list(bc_dict.keys()) if 
          str_mismatch(s, x) <= mm]
    if len(n) == 1:
        return n[0] # bc (from list)
    elif len(n) > 1:
        return 'multi'
    else:
        return None # undemx


def p7_validator(seq_unit, p7_dict, p7_len, mm=0):
    """Extract P7 index from comment field of fastq, eg: 1:N:0:NCAACAAT
    seq_unit : [name, seq, +, qual]
    return sample name
    """
    _name, _seq, _flag, _qual = seq_unit
    s = _name.split(':')[-1]
    p7_query = s[:p7_len]
    m = seq_validator(p7_query, p7_dict, mm=mm)
    return m


def bc_validator(seq_unit, bc_dict, bc_length, n_left=3, n_right=2, cut=True, mm=1):
    """
    extract barcode and randomer from read
    seq_unit : [name, seq, +, qual]
    append the randomer to the name
    """
    _name, _seq, _flag, _qual = seq_unit
    _x = (n_left + bc_length)
    _y = (n_left + bc_length + n_right)
    bc_query = _seq[n_left:_x]
    bc_random = _seq[0:n_left] + _seq[_x:_y]
    # check bc_dict
    m = seq_validator(bc_query, bc_dict, mm=mm) # barcode in list
    _name_list = _name.split(' ')
    _name_list.insert(1, bc_random)
    _name = ' '.join(_name_list)
    if cut and m:
        _seq = _seq[_y:]
        _qual = _qual[_y:]
    return [m, [_name, _seq, _flag, _qual]]


# def p7_bc_checker(p7, bc, p7_bc_dict, mm=0):
#     """Check whether p7 exists, barcode exists
#     allow no more than {mm} mismatche(s)
#     return the name
#     """
#     p7_list = []
#     bc_list = []
#     for k in list(p7_bc_dict.keys()):
#         if ',' in k:
#             p7, bc = k.upper().split(',')[0:2]
#             if not p7.lower() == 'null':
#                 p7_list.append(p7)
#             if not bc.lower() == 'null':
#                 bc_list.append(bc)
#             p7_list.append(p7)
#             bc_list.append(bc)
#         else:
#             continue
#     # unique records
#     p7_list = list(set(p7_list))
#     bc_list = list(set(bc_list))
#     # check mismatch
#     a = [x for x in p7_list if str_mismatch(p7, x) <= mm]
#     b1 = [x for x in bc_list if str_mismatch(bc, x) <= mm]
#     b2 = [x for x in bc_list if str_mismatch(bc, x) > 1000]
#     # check if barcode is null, return 1000
#     b = b2 if len(b2) == 1 else b1
#     # tag
#     tag = a[0] + ',' + b[0]
#     if len(a) == 1 and len(b) == 1:
#         return tag # index sequence
#     elif len(a) > 1 or len(b) > 1:
#         return 'multi'
#     else:
#         return None # undemx


def seq_in_array(s, array, mm=0):
    """Determine, whether the sequence s present in array
    compare two sequences, allow {mm} mismatches
    consider 2 situations:
    [1] null in array, check non-null values first, then null
    [2] null not in array, check non-null values
    return the hits in array
    return None, if 0 hits found
    """
    assert isinstance(s, str)
    assert isinstance(array, list)

    # hits
    hit_null = []
    hit_non_null = []
    for x in array:
        n = str_mismatch(s, x)
        if n <= mm:
            if x.lower() == 'null':
                hit_null.append(x.lower())
            else:
                hit_non_null.append(x)
    # unique hits
    hit_null = list(set(hit_null))
    hit_non_null = list(set(hit_non_null))
    # return number of hits
    if len(hit_non_null) > 0:
        return hit_non_null
    else:
        return hit_null


def p7_bc_validator(seq_unit, p7_bc_dict, n_left=3, n_right=2, cut=False, mm=0):
    """Demultiplexing both p7 and inline barcode at the same time
    allow no more than {mm} mismatche(s) for each
    return sample name
    """
    p7_list = []
    bc_list = []
    for k in list(p7_bc_dict.keys()):
        if ',' in k:
            p7, bc = k.upper().split(',')[0:2]
            if not p7.lower() == 'null':
                p7_list.append(p7)
            if not bc.lower() == 'null':
                bc_list.append(bc)
        else:
            continue
    # unique records
    p7_list = list(set(p7_list))
    bc_list = list(set(bc_list))
    # lenght
    p7_len = int(sum(list(map(len, p7_list))) / len(p7_list))
    bc_len = int(sum(list(map(len, bc_list))) / len(bc_list))
    # fastq sequence
    _name, _seq, _flag, _qual = seq_unit
    # p7
    s = _name.split(':')[-1]
    p7_query = s[:p7_len]
    # bc
    _x = (n_left + bc_len)
    _y = (n_left + bc_len + n_right)
    bc_query = _seq[n_left:_x]
    bc_random = _seq[0:n_left] + _seq[_x:_y]
    # check mismatch
    # !!! bug: both p7 and bc are 'null', 
    a = seq_in_array(p7_query, p7_list, mm)
    b = seq_in_array(bc_query, bc_list, mm)
    if len(a) == 1 and len(b) == 1:
        tag = a[0] + ',' + b[0]
        if tag in p7_bc_dict:
            m = tag
        else:
            m = None
    elif len(a) > 1 or len(b) > 1:
        m = 'multi'
    else:
        m = None # undemx
    ## update name
    _name_list = _name.split(" ")
    _name_list.insert(1, bc_random)
    _name = " ".join(_name_list)
    ## update sequence
    if cut is True and not m is None:
        _seq = _seq[_y:]
        _qual = _qual[_y:]
    return [m, [_name, _seq, _flag, _qual]]


##----------------------------------------------------------------------------##
## demultiplex illumina HiSeq TruSeq P7-index, 6-nt
## located in comment field of fastq file
def p7_demx_se(fn1, smp_list, path_out, n_mismatch=0, bc_in_read=1, n_left=3,
    n_right=2, cut=False, fn2=None):
    """Demultiplex reads by P7 index
    """
    p7_dict = bc_parser(smp_list, p7=True, barcode=False)
    p7_len = int(sum(len(x) for x in list(p7_dict.keys())) / len(p7_dict))
    p7_dict['undemx'] = 'undemx' # undemx file
    p7_dict['multi'] = 'multi' # barcode match multi hits
    
    ####################
    ## prepare writer ##
    ####################
    p7_count = {}
    p7_writer = {}
    for idx in p7_dict:
        p7_count[idx] = {}
        p7_count[idx]['count'] = 0
        p7_count[idx]['name'] = p7_dict[idx]
        p7_writer[idx] = gzip.open(os.path.join(path_out, 
                                    p7_dict[idx] + '.fq.gz'), 'wt')

    ####################
    ## prepare reader ##
    ####################
    fq_reader = gzip.open if is_gz(fn1) else open
    with fq_reader(fn1, 'rt') as fi:
        while True:
            try:
                seq_unit = [next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(),]
                tag = p7_validator(seq_unit, p7_dict, p7_len, n_mismatch)
                if tag is None:
                    tag = 'undemx'
                p7_writer[tag].write('\n'.join(seq_unit) + '\n')
                p7_count[tag]['count'] += 1
            except StopIteration:
                break
    # close file
    for idx in p7_writer: 
        p7_writer[idx].close() # close writers
    # bug: multi could not close
    p7_writer['multi'].close()
    # save report
    report_file = os.path.join(path_out, 'report_demx.json')
    Json_file(p7_count).json_writer(report_file)


def p7_demx_pe(fn1, fn2, smp_list, path_out, n_mismatch=0, bc_in_read=1, 
    n_left=3, n_right=2, cut=False):
    """Demultiplex PE reads, index in NAME
    """
    p7_dict = bc_parser(smp_list, p7=True, barcode=False)
    p7_len = int(sum(len(x) for x in list(p7_dict.keys())) / len(p7_dict))
    p7_dict['undemx'] = 'undemx' # undemx file
    p7_dict['multi'] = 'multi' # barcode match multi hits

    # # p7-index not correlate to read1 or read2
    # #####################
    # ## prepare rread12 ##
    # #####################
    r1_suffix = '_1.fq.gz'
    r2_suffix = '_2.fq.gz'
    # if bc_in_read == 2:
    #     (fn1, fn2) = (fn2, fn1) # switch read 1/2
    #     (r1_suffix, r2_suffix) = (r2_suffix, r1_suffix) # switch read 1/2


    ####################
    ## prepare writer ##
    ####################
    p7_count = {}
    p7_writer = {}
    for idx in p7_dict:
        p7_count[idx] = {}
        p7_count[idx]['count'] = 0
        p7_count[idx]['name'] = p7_dict[idx]
        p7_writer[idx] = [gzip.open(os.path.join(path_out, 
                              p7_dict[idx] + r1_suffix), 'wt'),
                          gzip.open(os.path.join(path_out, 
                              p7_dict[idx] + r2_suffix), 'wt'),]

    ####################
    ## prepare reader ##
    ####################
    fq_reader1 = gzip.open if is_gz(fn1) else open
    fq_reader2 = gzip.open if is_gz(fn2) else open
    with fq_reader1(fn1, 'rt') as f1, fq_reader2(fn2, 'rt') as f2:
        while True:
            try:
                seq_unit1 = [next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(),]
                seq_unit2 = [next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(),]
                # check read names from read1
                if not seq_unit1[0].split(' ')[0] == seq_unit2[0].split(' ')[0]:
                    raise NameError(seq_unit1[0] + '\n' + seq_unit2[0])
                tag = p7_validator(seq_unit1, p7_dict, p7_len, n_mismatch)
                if tag is None:
                    tag = 'undemx'
                p7_writer[tag][0].write('\n'.join(seq_unit1) + '\n')
                p7_writer[tag][1].write('\n'.join(seq_unit2) + '\n')
                p7_count[tag]['count'] += 1
            except StopIteration:
                break
    for idx in p7_writer: 
        p7_writer[idx][0].close() # close writers
        p7_writer[idx][1].close() # close writers
    # bug: multi could not close
    p7_writer['multi'][0].close()
    p7_writer['multi'][1].close()
    # save report
    report_file = os.path.join(path_out, 'report_demx.json')
    Json_file(p7_count).json_writer(report_file)


##----------------------------------------------------------------------------##
## in-line barcode dmultiplexing
## nnn-bc-nn
def bc_demx_se(fn1, smp_list, path_out, n_mismatch=0, bc_in_read=1, 
    n_left=3, n_right=2, cut=False, fn2=None):
    """Demultiplex SE reads
    inline-barcode at the beginning of read
    """
    bc_dict = bc_parser(smp_list, p7=False, barcode=True)
    bc_len = int(sum(len(x) for x in list(bc_dict.keys())) / len(bc_dict))
    bc_dict['undemx'] = 'undemx' # add undemx
    bc_dict['multi'] = 'multi' # barcode match multi hits

    ####################
    ## prepare writer ##
    ####################
    bc_count = {}
    bc_writer = {}
    for bc in bc_dict:
        bc_count[bc] = {}
        bc_count[bc]['count'] = 0
        bc_count[bc]['name'] = bc_dict[bc]
        bc_writer[bc] = gzip.open(os.path.join(path_out, 
                                               bc_dict[bc] + '.fq.gz'), 'wt')

    ####################
    ## prepare reader ##
    ####################
    fq_reader = gzip.open if is_gz(fn1) else open
    with fq_reader(fn1, 'rt') as fi:
        while True:
            try:
                seq_unit = [next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(),]
                tag, seq_unit_new = bc_validator(seq_unit, 
                                                 bc_dict, 
                                                 bc_len, 
                                                 n_left=n_left, 
                                                 n_right=n_right, 
                                                 cut=cut, 
                                                 mm=n_mismatch)
                if tag is None:
                    tag = 'undemx'
                bc_writer[tag].write('\n'.join(seq_unit_new) + '\n')
                bc_count[tag]['count'] += 1
            except StopIteration:
                break
    for bc in bc_writer: 
        bc_writer[bc].close() # close writers
    # bug: multi could not close
    bc_writer['multi'].close()
    # save report
    report_file = os.path.join(path_out, 'report_demx.json')
    Json_file(bc_count).json_writer(report_file)


def bc_demx_pe(fn1, fn2, smp_list, path_out, n_mismatch=0, bc_in_read=1, 
    n_left=3, n_right=2, cut=False):
    """Demultiplex PE reads
    barcode located in 5-prime end of read1/2
    with_bc: read1 | read2
    """
    bc_dict = bc_parser(smp_list, p7=False, barcode=True)
    bc_len = int(sum(len(x) for x in list(bc_dict.keys())) / len(bc_dict))
    bc_dict['undemx'] = 'undemx' # add undemx
    bc_dict['multi'] = 'multi' # barcode match multi hits

    ####################
    ## prepare read12 ##
    ####################
    r1_suffix = '_1.fq.gz'
    r2_suffix = '_2.fq.gz'
    if bc_in_read == 2:
        (fn1, fn2) = (fn2, fn1) # switch read 1/2
        (r1_suffix, r2_suffix) = (r2_suffix, r1_suffix) # switch read 1/2

    ####################
    ## prepare writer ##
    ####################
    bc_count = {}
    bc_writer = {}
    for bc in bc_dict:
        bc_count[bc] = {}
        bc_count[bc]['count'] = 0
        bc_count[bc]['name'] = bc_dict[bc]
        bc_writer[bc] = [gzip.open(os.path.join(path_out, 
                                   bc_dict[bc] + r1_suffix), 'wt'),
                         gzip.open(os.path.join(path_out, 
                                   bc_dict[bc] + r2_suffix), 'wt'),]

    ####################
    ## prepare reader ##
    ####################
    fq_reader1 = gzip.open if is_gz(fn1) else open
    fq_reader2 = gzip.open if is_gz(fn2) else open
    with fq_reader1(fn1, 'rt') as f1, fq_reader2(fn2, 'rt') as f2:
        while True:
            try:
                seq_unit1 = [next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(),]
                seq_unit2 = [next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(),]
                if not seq_unit1[0].split(' ')[0] == seq_unit2[0].split(' ')[0]:
                    raise NameError(seq_unit1[0] + '\n' + seq_unit2[0])
                tag, seq_unit1 = bc_validator(seq_unit1, 
                                              bc_dict, 
                                              bc_len, 
                                              n_left=n_left, 
                                              n_right=n_right, 
                                              cut=cut, 
                                              mm=n_mismatch)
                if tag is None:
                    tag = 'undemx'
                bc_writer[tag][0].write('\n'.join(seq_unit1) + '\n')
                bc_writer[tag][1].write('\n'.join(seq_unit2) + '\n')
                bc_count[tag]['count'] += 1
            except StopIteration:
                break
    for bc in bc_writer: 
        bc_writer[bc][0].close() # close writers
        bc_writer[bc][1].close() # close writers
    # bug: multi could not close
    bc_writer['multi'][0].close()
    bc_writer['multi'][1].close()
    # save report
    report_file = os.path.join(path_out, 'report_demx.json')
    Json_file(bc_count).json_writer(report_file)


##----------------------------------------------------------------------------##
## both p7-index and in-line barcode dmultiplexing
## nnn-bc-nn, p7-index
def p7_bc_demx_se(fn1, smp_list, path_out, n_mismatch=0, bc_in_read=1,
    n_left=3, n_right=2, cut=False, fn2=None):
    """Demultiplex reads by P7 index and barcode
    """
    p7_bc_dict = bc_parser(smp_list, p7=True, barcode=True)
    p7_bc_dict['undemx'] = 'undemx' # undemx file
    p7_bc_dict['multi'] = 'multi' # barcode match multi hits

    ##################
    # prepare writer #
    ##################
    p7_bc_count = {}
    p7_bc_writer = {}
    for p7_bc in p7_bc_dict:
        p7_bc_count[p7_bc] = {}
        p7_bc_count[p7_bc]['count'] = 0
        p7_bc_count[p7_bc]['name'] = p7_bc_dict[p7_bc]
        p7_bc_writer[p7_bc] = gzip.open(os.path.join(path_out,
                                        p7_bc_dict[p7_bc] + '.fq.gz'), 'wt')

    ##################
    # prepare reader #
    ##################
    fq_reader = gzip.open if is_gz(fn1) else open
    with fq_reader(fn1, 'rt') as fi:
        while True:
            try:
                seq_unit = [next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(),]
                tag, seq_unit_new = p7_bc_validator(seq_unit, 
                                                    p7_bc_dict, 
                                                    n_left=n_left,
                                                    n_right=n_right,
                                                    cut=cut,
                                                    mm=n_mismatch) 
                if tag is None:
                    tag = 'undemx'
                p7_bc_writer[tag].write('\n'.join(seq_unit_new) + '\n')
                p7_bc_count[tag]['count'] += 1
            except StopIteration:
                break
    for p7_bc in p7_bc_writer: 
        p7_bc_writer[p7_bc].close() # close writers
    # bug: multi could not close
    p7_bc_writer['multi'].close()
    # save report
    report_file = os.path.join(path_out, "report_demx.json")
    Json_file(p7_bc_count).json_writer(report_file)


def p7_bc_demx_pe(fn1, fn2, smp_list, path_out, n_mismatch=0, bc_in_read=1,
    n_left=3, n_right=2, cut=False):
    """Demultiplex PE reads, by p7 index and barcode
    p7 in comment of fastq
    barcode in fn1
    """
    p7_bc_dict = bc_parser(smp_list, p7=True, barcode=True)
    p7_bc_dict['undemx'] = 'undemx' # undemx file
    p7_bc_dict['multi'] = 'multi' # barcode match multi hits

    ####################
    ## prepare read12 ##
    ####################
    r1_suffix = '_1.fq.gz'
    r2_suffix = '_2.fq.gz'
    if bc_in_read == 2:
        (fn1, fn2) = (fn2, fn1) # switch read 1/2
        (r1_suffix, r2_suffix) = (r2_suffix, r1_suffix) # switch read 1/2

    ####################
    ## prepare writer ##
    ####################
    p7_bc_count = {}
    p7_bc_writer = {}
    for p7_bc in p7_bc_dict:
        p7_bc_count[p7_bc] = {}
        p7_bc_count[p7_bc]['count'] = 0
        p7_bc_count[p7_bc]['name'] = p7_bc_dict[p7_bc]
        p7_bc_writer[p7_bc] = [gzip.open(os.path.join(path_out,
            p7_bc_dict[p7_bc] + r1_suffix), 'wt'),
        gzip.open(os.path.join(path_out,
            p7_bc_dict[p7_bc] + r2_suffix), 'wt'),]
        
    ####################
    ## prepare reader ##
    ####################
    fq_reader1 = gzip.open if is_gz(fn1) else open # contain barcode in first file
    fq_reader2 = gzip.open if is_gz(fn2) else open
    with fq_reader1(fn1, 'rt') as f1, fq_reader2(fn2, 'rt') as f2:
        while True:
            try:
                seq_unit1 = [next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(),]
                seq_unit2 = [next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(),]
                if not seq_unit1[0].split(" ")[0] == seq_unit2[0].split(" ")[0]:
                    raise NameError(seq_unit1[0] + "\n" + seq_unit2[0])
                tag, seq_unit1_new = p7_bc_validator(seq_unit1, 
                                                     p7_bc_dict,
                                                     n_left=n_left,
                                                     n_right=n_right,
                                                     cut=cut,
                                                     mm=n_mismatch) 
                if tag is None:
                    tag = 'undemx'
                p7_bc_writer[tag][0].write('\n'.join(seq_unit1_new) + '\n')
                p7_bc_writer[tag][1].write('\n'.join(seq_unit2) + '\n')
                p7_bc_count[tag]['count'] += 1
            except StopIteration:
                break
    for idx in p7_bc_writer: 
        p7_bc_writer[idx][0].close() # close writers
        p7_bc_writer[idx][1].close() # close writers
    # bug: multi could not close
    p7_bc_writer['multi'][0].close()
    p7_bc_writer['multi'][1].close()
    # save report
    report_file = os.path.join(path_out, "report_demx.json")
    Json_file(p7_bc_count).json_writer(report_file)


##----------------------------------------------------------------------------##
def main():
    """Run demx"""
    args = get_args()
    read1 = args.r1.name
    read2 = args.r2.name if args.r2 else None
    smp_list = args.sample_list
    path_out = args.out
    demx_p7 = args.p7_index
    demx_barcode = args.barcode
    n_mismatch = args.n_mismatch
    bc_in_read = args.bc_in_read
    n_left = args.bc_n_left
    n_right = args.bc_n_right
    cut = args.cut # whether cut the barcode
   
    #######################
    ## prepare work dirs ##
    #######################
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    ## determine demx programe
    if demx_p7 is True and demx_barcode is True:
        if read2 is None:
            demx_tool = p7_bc_demx_se
        else:
            demx_tool = p7_bc_demx_pe
    elif demx_p7 is True:
        if read2 is None:
            demx_tool = p7_demx_se
        else:
            demx_tool = p7_demx_pe
    elif demx_barcode is True:
        if read2 is None:
            demx_tool = bc_demx_se
        else:
            demx_tool = bc_demx_pe
    else:
        sys.exit('--p7-index and --barcode, either or both required')

    ## time
    d = datetime.datetime.today()
    t1 = d.strftime('%Y-%m-%d %H:%M:%S')

    ## run demx
    tmp = demx_tool(fn1=read1, 
                    fn2=read2,
                    smp_list=smp_list,
                    path_out=path_out,
                    n_mismatch=n_mismatch,
                    bc_in_read=bc_in_read,
                    n_left=n_left,
                    n_right=n_right,
                    cut=cut)
    ## report
    logging.info('demx finish!')

    ## time
    d = datetime.datetime.today()
    t2 = d.strftime('%Y-%m-%d %H:%M:%S')

    ## save parameters
    run_lib = os.path.join(path_out, "run_demx.lib")
    run_lib_list = [
        '%25s: %s' % ('Programe', 'hipipe-demx'),
        '%25s: %s' % ('Start', t1),
        '%25s: %s' % ('Finish', t2),
        '%25s: %s' % ('Read1', read1),
        '%25s: %s' % ('Read2', read2),
        '%25s: %s' % ('Barcode_list', smp_list),
        '%25s: %s' % ('Output_dir', path_out),
        '%25s: %s' % ('Demx_P7', demx_p7),
        '%25s: %s' % ('Demx_barcode', demx_barcode),
        '%25s: %d' % ('Number_of_mismatch', n_mismatch),
        '%25s: %s' % ('Barcode_in_read', bc_in_read),
        '%25s: %s' % ('Randomer_in_barcode_left', n_left),
        '%25s: %s' % ('Randomer_in_barcode_right', n_right),
        '%25s: %s' % ('Cut_barcode', cut)]
    with open(run_lib, 'wt') as fo:
        fo.write('\n'.join(run_lib_list) + '\n')

if __name__ == '__main__':
    main()


## EOF
