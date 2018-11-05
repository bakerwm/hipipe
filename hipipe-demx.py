#!/usr/bin/env python3




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
import Bio
from Bio import SeqIO

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

def is_gz(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


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


def str_mismatch(q, s):
    """Return the number of mismatches between two strings
    if subject is 'null', return 0
    """
    assert isinstance(q, str)
    assert isinstance(s, str)
    m = sum(c1!=c2 for c1,c2 in zip(q, s))
    if s.lower() == 'null':
        m = 0
    return m


def bc_parser(fn, p7=True, barcode=True):
    """Parse the p7-index, barcode and sample name from sample list
    return a dict
    require: <p7-index> <barcode> <sample_name>
    """
    assert os.path.isfile(fn)
    assert isinstance(p7, bool)
    assert isinstance(barcode, bool)
    p7_dict = {}
    bc_dict = {}
    p7_bc_dict = {}
    n_rec = 0
    # parse file
    with open(fn, 'rt') as fi:
        for line in fi:
            p = line.strip().split('\t') # <p7> <barcode> <name>
            if not len(p) == 3:
                sys.exit('require <p7> <barcode> <name> in sample list')
            p[0] = p[0].upper()
            p[1] = p[1].upper()
            tag = p[0] + ',' + p[1]
            if not p[0] == 'NULL':
                p7_dict[p[0]] = p[2]
            if not p[1] == 'NULL':
                bc_dict[p[1]] = p[2]
            p7_bc_dict[tag] = p[2]
            n_rec += 1 # number of records
    # check duplicates
    if p7 is True and barcode is True:
        if len(p7_bc_dict) < n_rec:
            raise ValueError('p7 + barcode duplicates detected')
    elif p7 is True:
        if len(p7_dict) < n_rec or 'NULL' in p7_dict:
            raise ValueError('p7 duplicates or NULL detected')
    elif barcode is True:
        if len(bc_dict) < n_rec or 'NULL' in bc_dict:
            raise ValueError('barcode duplicates or NULL detected')
    else:
        sys.exit('either or both p7 and barcode are required')
    # check barcode length
    n1 = set(list(map(len, list(p7_dict.keys()))))
    n2 = set(list(map(len, list(bc_dict.keys()))))
    if len(n1) > 1 or len(n2) > 1:
        raise ValueError('p7, barcode not in the same length')
    # check sample name
    _name = sorted(list(p7_bc_dict.values()))
    _name_uniq = sorted(list(set(_name)))
    if not _name == _name_uniq:
        raise ValueError('Sample_name duplicates detected')
    # return values
    return [p7_dict, bc_dict, p7_bc_dict]



def seq_index_finder(s, cut=False, bc_len=0, bc_n_left=0, bc_n_right=0, mode=0):
    """Extract the barcode, p7-index from sequence
    Processing: s = Bio.SeqRecord.SeqRecord
    p7 in the tail of name: ST-E00310:586:HJH7JCCXY:7:1101:27275:1555 1:N:0:NCAACAAT
    barcode in the left end of seq:
    cut: whether cut the barcode
    # add
    mode: 0=p7, 1=barcode, 2=both
    """
    # assert isinstance(s, Bio.SeqRecord.SeqRecord)
    assert isinstance(s, list)
    assert isinstance(cut, bool)
    assert isinstance(bc_len, int) and bc_len >= 0
    assert isinstance(bc_n_left, int) and bc_n_left >= 0
    assert isinstance(bc_n_right, int) and bc_n_right >= 0
    s_description, s_seq, _, s_qual = s
    # extract index
    p7_query = s_description.split(':')[-1]
    bc_query = None
    # extract barcode
    if bc_len > 0:
        _x = (bc_n_left + bc_len)
        _y = (bc_n_left + bc_len + bc_n_right)
        bc_query = str(s_seq[bc_n_left:_x])
        bc_random = str(s_seq[:bc_n_left]) + str(s_seq[_x:_y])
        if cut is True:
            _name = s_description.split(' ')
            _name.insert(1, bc_random)
            s_description = ' '.join(_name)
            s_seq = s_seq[_y:]
            s_qual = s_qual[_y:]

    return [[s_description, s_seq, '+', s_qual], p7_query, bc_query]

# def seq_index_finder(s, cut=False, bc_len=0, bc_n_left=0, bc_n_right=0):
#     """Extract the barcode, p7-index from sequence
#     Processing: s = Bio.SeqRecord.SeqRecord
#     p7 in the tail of name: ST-E00310:586:HJH7JCCXY:7:1101:27275:1555 1:N:0:NCAACAAT
#     barcode in the left end of seq:
#     cut: whether cut the barcode
#     """
#     assert isinstance(s, Bio.SeqRecord.SeqRecord)
#     assert isinstance(cut, bool)
#     assert isinstance(bc_len, int) and bc_len >= 0
#     assert isinstance(bc_n_left, int) and bc_n_left >= 0
#     assert isinstance(bc_n_right, int) and bc_n_right >= 0
#     # extract index
#     p7_query = s.description.split(':')[-1]
#     bc_query = None
#     # extract barcode
#     if bc_len > 0:
#         _x = (bc_n_left + bc_len)
#         _y = (bc_n_left + bc_len + bc_n_right)
#         bc_query = str(s.seq[bc_n_left:_x])
#         bc_random = str(s.seq[:bc_n_left]) + str(s.seq[_x:_y])
#         if cut is True:
#             _name = s.description.split(' ')
#             _name.insert(1, bc_random)
#             s.description = ' '.join(_name)
#             s = s[_y:]

#     return [s, p7_query, bc_query]


def seq_validator(q, s_dict, mm=0):
    """Validate the barcode or P7-index in dictionary
    check number of mismatch
    if null in dict, return 0
    """
    assert isinstance(q, str)
    assert isinstance(s_dict, dict)
    assert isinstance(mm, int)
    # match in dict
    n = [x for x in list(s_dict.keys()) if
        str_mismatch(q, x) <= mm]
    if len(n) == 1:
        return n[0] # the barcode in list
    elif len(n) > 1:
        return 'multi'
    else:
        return None # undemx


def p7_bc_validator(p7, bc, p7_dict, bc_dict, p7_bc_dict, mm=0):
    """Validate the p7-index AND barcode at the same time
    allow mm mismatch
    return: None, multi, index
    """
    # p7_dict, bc_dict
    tag1 = seq_validator(p7, p7_dict, mm=mm)
    tag2 = seq_validator(bc, bc_dict, mm=mm)
    # check tags
    if tag1 == 'multi' or tag2 == 'multi':
        tag = 'multi'
    else:
        if tag1 is None:
            if not tag2 is None:
                if 'NULL,' + tag2 in p7_bc_dict:
                    # print('AAA')
                    tag = 'NULL,' + tag2
                else:
                    # print('BBB')
                    tag = 'undemx'
            else:
                # print('CCC')
                tag = 'undemx'
        else:
            if not tag2 is None:
                if tag1 + ',' + tag2 in p7_bc_dict:
                    # print('DDD')
                    tag = tag1 + ',' + tag2
                else:
                    # print('EEE')
                    tag = 'undemx'
            else:
                if tag1 + ',NULL' in p7_bc_dict:
                    # print('FFF')
                    tag = tag1 + ',NULL'
                else:
                    # print('GGG')
                    tag = 'undemx'
    return tag


def p7_demx_se(fn1, smp_list, path_out, n_mismatch=0, bc_in_read=1, n_left=3,
    n_right=2, cut=False, fn2=None):
    p7_dict = bc_parser(smp_list, p7=True, barcode=False)[0]
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
        p7_writer[idx] = open(os.path.join(path_out, 
            p7_dict[idx] + '.fq'), 'wt')
    ####################
    ## prepare reader ##
    ####################
    fq_reader = gzip.open if is_gz(fn1) else open
    with fq_reader(fn1, 'rt') as f1:
        while True:
            try:
                s1 = [next(f1), next(f1), next(f1), next(f1)]
                s1 = [i.strip() for i in s1]
                s1_new, p7_query, bc_query = seq_index_finder(s1, cut=cut, 
                    bc_len=0, bc_n_left=0, bc_n_right=0)
                tag = seq_validator(p7_query, p7_dict, mm=n_mismatch)
                if tag is None: #undemx
                    tag = 'undemx'
                p7_writer[tag].write('\n'.join(s1_new) + '\n')
                p7_count[tag]['count'] += 1
            except StopIteration:
                break
    # close writers
    for idx in p7_writer: 
        p7_writer[idx].close() # close writers

    # save report
    report_file = os.path.join(path_out, 'report_demx.json')
    Json_file(p7_count).json_writer(report_file)



def p7_demx_pe(fn1, fn2, smp_list, path_out, n_mismatch=0, bc_in_read=1, 
    n_left=3, n_right=2, cut=False):
    """Demultiplex PE reads, index in NAME
    """
    p7_dict = bc_parser(smp_list, p7=True, barcode=False)[0]
    p7_len = int(sum(len(x) for x in list(p7_dict.keys())) / len(p7_dict))
    p7_dict['undemx'] = 'undemx' # undemx file
    p7_dict['multi'] = 'multi' # barcode match multi hits

    # # p7-index not correlate to read1 or read2
    # #####################
    # ## prepare rread12 ##
    # #####################
    r1_suffix = '_1.fq'
    r2_suffix = '_2.fq'

    ####################
    ## prepare writer ##
    ####################
    p7_count = {}
    p7_writer = {}
    for idx in p7_dict:
        p7_count[idx] = {}
        p7_count[idx]['count'] = 0
        p7_count[idx]['name'] = p7_dict[idx]
        p7_writer[idx] = [open(os.path.join(path_out, 
                              p7_dict[idx] + r1_suffix), 'wt'),
                          open(os.path.join(path_out, 
                              p7_dict[idx] + r2_suffix), 'wt'),]

    ####################
    ## prepare reader ##
    ####################
    fq_reader1 = gzip.open if is_gz(fn1) else open
    fq_reader2 = gzip.open if is_gz(fn2) else open
    with fq_reader1(fn1, 'rt') as f1, fq_reader2(fn2, 'rt') as f2:
        while True:
            try:
                s1 = [next(f1), next(f1), next(f1), next(f1)]
                s2 = [next(f2), next(f2), next(f2), next(f2)]
                s1 = [i.strip() for i in s1]
                s2 = [i.strip() for i in s2]
                # check read names from read1
                if not s1[0].split(' ')[0] == s2[0].split(' ')[0]:
                    raise ValueError(s1[0] + '\n' + s2[0])
                # default bc in s1
                s1_new, p7_query, bc_query = seq_index_finder(s1, cut=cut,
                    bc_len=0, bc_n_left=0, bc_n_right=0)
                tag = seq_validator(p7_query, p7_dict, mm=n_mismatch)
                if tag is None: #undemx
                    tag = 'undemx'
                # write fastq
                p7_writer[tag][0].write('\n'.join(s1_new) + '\n')
                p7_writer[tag][1].write('\n'.join(s2) + '\n')
                p7_count[tag]['count'] += 1
            except StopIteration:
                break
    # close writers
    for idx in p7_writer: 
        p7_writer[idx][0].close() # close writers
        p7_writer[idx][1].close() # close writers
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
    bc_dict = bc_parser(smp_list, p7=False, barcode=True)[1]
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
        bc_writer[bc] = open(os.path.join(path_out, 
            bc_dict[bc] + '.fq'), 'wt')

    ####################
    ## prepare reader ##
    ####################
    fq_reader = gzip.open if is_gz(fn1) else open
    with fq_reader(fn1, 'rt') as f1:
        while True:
            try:
                s1 = [next(f1), next(f1), next(f1), next(f1)]
                s1 = [i.strip() for i in s1]
                s1_new, p7_query, bc_query = seq_index_finder(s1, cut=cut, 
                    bc_len=bc_len, bc_n_left=n_left, bc_n_right=n_right)
                tag = seq_validator(bc_query, bc_dict, mm=n_mismatch)
                if tag is None: #undemx
                    tag = 'undemx'
                bc_writer[tag].write('\n'.join(s1_new) + '\n')
                bc_count[tag]['count'] += 1
            except StopIteration:
                break
    # close writers
    for idx in bc_writer: 
        bc_writer[idx].close() # close writers

    # save report
    report_file = os.path.join(path_out, 'report_demx.json')
    Json_file(bc_count).json_writer(report_file)


def bc_demx_pe(fn1, fn2, smp_list, path_out, n_mismatch=0, bc_in_read=1, 
    n_left=3, n_right=2, cut=False):
    """Demultiplex PE reads
    barcode located in 5-prime end of read1/2
    with_bc: read1 | read2
    """
    bc_dict = bc_parser(smp_list, p7=False, barcode=True)[1]
    bc_len = int(sum(len(x) for x in list(bc_dict.keys())) / len(bc_dict))
    bc_dict['undemx'] = 'undemx' # add undemx
    bc_dict['multi'] = 'multi' # barcode match multi hits

    ####################
    ## prepare read12 ##
    ####################
    r1_suffix = '_1.fq'
    r2_suffix = '_2.fq'
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
                s1 = [next(f1), next(f1), next(f1), next(f1)]
                s2 = [next(f2), next(f2), next(f2), next(f2)]
                s1 = [i.strip() for i in s1]
                s2 = [i.strip() for i in s2]
                # check read names from read1
                if not s1[0].split(' ')[0] == s2[0].split(' ')[0]:
                    raise ValueError(s1[0] + '\n' + s2[0])
                # default bc in s1
                s1_new, p7_query, bc_query = seq_index_finder(s1, cut=cut,
                    bc_len=bc_len, bc_n_left=n_left, bc_n_right=n_right)
                tag = seq_validator(bc_query, bc_dict, mm=n_mismatch)
                if tag is None: #undemx
                    tag = 'undemx'
                # write fastq
                bc_writer[tag][0].write('\n'.join(s1_new) + '\n')
                bc_writer[tag][1].write('\n'.join(s2) + '\n')
                bc_count[tag]['count'] += 1
            except StopIteration:
                break
    # close writers
    for bc in bc_writer: 
        bc_writer[bc][0].close() # close writers
        bc_writer[bc][1].close() # close writers
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
    n_dict = bc_parser(smp_list, p7=True, barcode=True)
    p7_dict = n_dict[0]
    bc_dict = n_dict[1]
    bc_len = int(sum(len(x) for x in list(bc_dict.keys())) / len(bc_dict))
    p7_bc_dict = n_dict[2]
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
        p7_bc_writer[p7_bc] = open(os.path.join(path_out,
                                        p7_bc_dict[p7_bc] + '.fq'), 'wt')

    ##################
    # prepare reader #
    ##################
    fq_reader = gzip.open if is_gz(fn1) else open
    with fq_reader(fn1, 'rt') as f1:
        while True:
            try:
                s1 = [next(f1), next(f1), next(f1), next(f1)]
                s1 = [i.strip() for i in s1]
                # check index
                s1_new, p7_query, bc_query = seq_index_finder(s1, cut=cut,
                    bc_len=bc_len, bc_n_left=n_left, bc_n_right=n_right)
                tag = p7_bc_validator(p7_query, bc_query, p7_dict, bc_dict, 
                    p7_bc_dict, mm=n_mismatch)
                # write fastq
                p7_bc_writer[tag].write('\n'.join(s1_new) + '\n')
                p7_bc_count[tag]['count'] += 1
            except StopIteration:
                break
    # close writers
    for p7_bc in p7_bc_writer: 
        p7_bc_writer[p7_bc].close() # close writers
    # save report
    report_file = os.path.join(path_out, "report_demx.json")
    Json_file(p7_bc_count).json_writer(report_file)


def p7_bc_demx_pe(fn1, fn2, smp_list, path_out, n_mismatch=0, bc_in_read=1,
    n_left=3, n_right=2, cut=False):
    """Demultiplex PE reads, by p7 index and barcode
    p7 in comment of fastq
    barcode in fn1
    """
    n_dict = bc_parser(smp_list, p7=True, barcode=True)
    p7_dict = n_dict[0]
    bc_dict = n_dict[1]
    bc_len = int(sum(len(x) for x in list(bc_dict.keys())) / len(bc_dict))
    p7_bc_dict = n_dict[2]
    p7_bc_dict['undemx'] = 'undemx' # undemx file
    p7_bc_dict['multi'] = 'multi' # barcode match multi hits

    ####################
    ## prepare read12 ##
    ####################
    r1_suffix = '_1.fq'
    r2_suffix = '_2.fq'
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
        p7_bc_writer[p7_bc] = [open(os.path.join(path_out,
            p7_bc_dict[p7_bc] + r1_suffix), 'wt'),
                               open(os.path.join(path_out,
            p7_bc_dict[p7_bc] + r2_suffix), 'wt'),]
        
    ####################
    ## prepare reader ##
    ####################
    fq_reader1 = gzip.open if is_gz(fn1) else open 
    fq_reader2 = gzip.open if is_gz(fn2) else open
    with fq_reader1(fn1, 'rt') as f1, fq_reader2(fn2, 'rt') as f2:
        while True:
            try:
                s1 = [next(f1), next(f1), next(f1), next(f1)]
                s2 = [next(f2), next(f2), next(f2), next(f2)]
                s1 = [i.strip() for i in s1]
                s2 = [i.strip() for i in s2]
                # check read names from read1
                if not s1[0].split(' ')[0] == s2[0].split(' ')[0]:
                    raise ValueError(s1[0] + '\n' + s2[0])
                # check name
                s1_new, p7_query, bc_query = seq_index_finder(s1, cut=cut,
                    bc_len=bc_len, bc_n_left=n_left, bc_n_right=n_right)
                tag = p7_bc_validator(p7_query, bc_query, p7_dict, bc_dict, p7_bc_dict, mm=n_mismatch)
                # write fastq
                p7_bc_writer[tag][0].write('\n'.join(s1_new) + '\n')
                p7_bc_writer[tag][1].write('\n'.join(s2) + '\n')
                p7_bc_count[tag]['count'] += 1
            except StopIteration:
                break
    # close writers
    for idx in p7_bc_writer: 
        p7_bc_writer[idx][0].close() # close writers
        p7_bc_writer[idx][1].close() # close writers
    # save report
    report_file = os.path.join(path_out, "report_demx.json")
    Json_file(p7_bc_count).json_writer(report_file)



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