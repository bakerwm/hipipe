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
import logging




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


def key_len(d):
    """Length of keys in dict
    skip null, '' 
    return len
    """
    k = list(d.keys())
    k = [len(i) for i in k if len(i) > 0 and not i.upper() == 'NULL']
    n = list(set(k))
    if len(n) == 1:
        return n[0]
    else:
        return None


def p7_parser(fn):
    """Parse the p7-index
    return dict
    """
    assert os.path.isfile(fn)
    p7_dict = {}
    with open(fn, 'rt') as fi:
        for line in fi:
            p = line.strip().split('\t') # <p7> <barcode> <name>
            if not len(p) == 3:
                sys.exit('require <p7> <barcode> <name> in sample list')
            p[0] = p[0].upper()
            if p[0] in p7_dict:
                raise ValueError('p7 index duplecates detected')
            p7_dict[p[0]] = p[2]
    return p7_dict


def bc_parser(fn):
    """Parse the barcode
    return dict
    """
    assert os.path.isfile(fn)
    bc_dict = {}
    with open(fn, 'rt') as fi:
        for line in fi:
            p = line.strip().split('\t') # <p7> <barcode> <name>
            if not len(p) == 3:
                sys.exit('require <p7> <barcode> <name> in sample list')
            p[1] = p[1].upper()
            if p[1] in bc_dict:
                raise ValueError('p7 index duplecates detected')
            bc_dict[p[1]] = p[2]
    return bc_dict



def seq_index_finder(s, cut=False, bc_len=0, bc_n_left=0, bc_n_right=0, mode=0):
    """Extract the barcode, p7-index from sequence
    : arg s, a list contins fastq elements
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
    #s_description, s_seq, _, s_quality = s
    s_description, s_seq, _, s_quality = s
    # extract index
    p7_query = s_description.split(':')[-1]
    bc_query = None
    # extract barcode
    if bc_len > 0:
        _x = (bc_n_left + bc_len)
        _y = (bc_n_left + bc_len + bc_n_right)
        bc_query = s_seq[bc_n_left:_x]
        bc_random = str(s_seq[:bc_n_left]) + str(s_seq[_x:_y])
        if cut is True:
            _name = s_description.split(' ')
            _name.insert(1, bc_random)
            s_description = ' '.join(_name)
            s_seq = s_seq[_y:]
            s_quality = s_quality[_y:]

    return [(s_description, s_seq, '+', s_quality), p7_query, bc_query]


def seq_validator(q, s_dict, mm=0):
    """Validate the barcode or P7-index in dictionary
    compare 'binary' mode
    check number of mismatch
    if null in dict, return 0
    # bug:
    !!!multiple hits!!!!
    """
    assert isinstance(q, str)
    assert isinstance(s_dict, dict)
    assert isinstance(mm, int)
    # match in dict
    # n = [x for x in list(s_dict.keys()) if
    #     str_mismatch(q, x) <= mm]
    # create tag
    for x in list(s_dict.keys()):
        m = str_mismatch(q, x)
        if m <= mm:
            tag = x
            break
        else:
            tag = 'undemx'

    return tag


def p7_demx_pe(fn1, fn2, smp_list, path_out, n_mismatch=0, bc_in_read=1, 
    n_left=3, n_right=2, cut=False):
    """Demultiplex PE reads, index in NAME
    """
    p7_dict = p7_parser(smp_list)
    p7_len = key_len(p7_dict)
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



def bc_demx_pe(fn1, fn2, smp_list, path_out, n_mismatch=0, bc_in_read=1, 
    n_left=3, n_right=2, cut=False):
    """Demultiplex PE reads
    barcode located in 5-prime end of read1/2
    with_bc: read1 | read2
    """
    bc_dict = bc_parser(smp_list)
    bc_len = key_len(bc_dict)
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
        bc_writer[bc] = [open(os.path.join(path_out, 
                                   bc_dict[bc] + r1_suffix), 'wt'),
                         open(os.path.join(path_out, 
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


# def idx_reader(x):
#     """Read index file
#     1. P7-index
#     2. inline barcode
#     Two steps
    
#     File format:

#     <p7-index>  <bc>    <sample_name>

#     """

#     p7_dict = {}
#     bc_dict = {}
#     dd_dict = {} # check unique


#     idx = []
#     with open(x, 'wt') as fi:
#         for line in fi:
#             # idx.append(line.strip())
#             ps = line.strip().split('\t')
#             if not len(ps) == 3:
#                 raise Exception('unknown file format - %s' % line)
#             p7, bc, name = ps
#             if p7 in p7_dict:
#                 bc_dict[p7] = p7
#                 p7_dict[p7] = '\t'.join([p7, bc, p7])
#             else:
#                 p7_dict[p7] = line.strip()

#     # iterate 2nd
#     for i in idx:
        
        



    # split p7, bc






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
    # print('BBB: n_mismatch: %s' % n_mismatch)
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