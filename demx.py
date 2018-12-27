#!/usr/bin/env python


__author__ = 'Ming Wang'
__email__ = 'wangm08@hotmail.com'
__date__ = '2018-12-25'
__version__ = '0.3'


import os, sys
import gzip, binascii, json


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


    def json_writer(self, to=None, sorted=True):
        """Write dict to file in json format"""
        fn = self.fn

        if to is None:
            to = self._tmp()

        if isinstance(fn, Json_file):
            fn = fn.fn
        elif isinstance(fn, dict):
            fn = fn
        with open(to, 'wt') as ff:
            json.dump(fn, ff, indent=4, sort_keys=sorted)


class Demx(object):
    """Demultiplexing Illumina SE and PE reads according to the 
    P7-index and inline barcode

    The library was constructed according to Illumina TruSeq protocol and 
    P7-index was saved in fastq 'name field'

    The inline barcode was designed for customized library, barcode was 
    located in the beginning of read1 or read2 which was consist with 
    randomer and specific barcode, [randomer] + [barcode] + [randomer]

    The sample_info should be in the following format with 3 columns:
    <p7-index> <bc-index> <name>

    Example:

    fq1 = 'data/raw_data/demo_01m_1.fq.gz'
    fq2 = 'data/raw_data/demo_01m_2.fq.gz'
    sample_info = 'data/sample_info/samplesheet.txt'
    path_out = 'results/demo'

    Demx(fq1=fq1, fq2=fq2, sample_info=sample_info, path_out=path_out,
        n_mismatch=2).p7_demx_pe()

    Demx(fq1=fq1, fq2=fq2, sample_info=sample_info, path_out=path_out,
        n_mismatch=2).p7_demx_se()

    """

    def __init__(self, fq1, sample_info, path_out, **kwargs):
        self.fq1 = fq1
        self.sample_info = sample_info
        self.path_out = path_out
        self.kwargs = kwargs
        if not os.path.exists(path_out):
            os.makedirs(path_out)
        # wrap args
        self.args = self._args_init()


    def _args_init(self):
        args = self.kwargs
        args['fq1'] = self.fq1
        args['sample_info'] = self.sample_info
        args['path_out'] = self.path_out
        args['fq2'] = args.get('fq2', None)
        args['n_mismatch'] = args.get('n_mismatch', 0)
        args['bc_in_read'] = args.get('bc_in_read', 1)
        args['n_left'] = args.get('n_left', 0)
        args['n_right'] = args.get('n_right', 0)
        args['cut'] = args.get('cut', False)
        return args


    def key_len(self, d):
        """Return the length of keys in dict
        skip null, '' 
        """
        assert isinstance(d, dict)
        k = [len(i) for i in list(d.keys()) if len(i) > 0 and not i.upper() == 'NULL']
        n = list(set(k))
        if len(n) == 1:
            return n[0]
        else:
            return None


    def str_mismatch(self, q, s):
        """Return the number of mismatches between two strings
        if subject is 'null', return 0
        """
        assert isinstance(q, str)
        assert isinstance(s, str)
        m = sum(c1!=c2 for c1,c2 in zip(q, s))
        if s.lower() == 'null':
            m = 0
        return m            


    def sampleinfo_spliter(self):
        """Split the origin sample_info file by the P7-barcode, 
        File format:
        <p7-index>  <bc>    <sample_name>

        split index file according to the p7-index

        1. p7-index only
        2. bc-1
        3. bc-2 ...

        """
        x = self.sample_info

        # read file
        idx = []
        p7_dict = {}
        p7_bc_dict = {}

        with open(x, 'rt') as fi:
            for line in fi:
                idx.append(line.strip())
                if ' ' in line:
                    sep = ' '
                elif '\t' in line:
                    sep = '\t'
                else:
                    raise Exception('unknown file format, expect <TAB> separated - %s' % line)
                ps = line.strip().split(sep)
                if not len(ps) == 3:
                    raise Exception('unknown file format, expect 3 columns- %s' % line)
                p7, bc, name = ps
                if p7 in p7_dict:
                    p7_bc_dict[p7] = p7
                    p7_dict[p7] = '\t'.join([p7, bc, p7])
                else:
                    p7_dict[p7] = line.strip()

        # save p7 index
        p7_list = list(p7_dict.values())
            
        # save bc index
        bc_list = []
        for p7_bc in p7_bc_dict:
            tmp_list = []
            for i in idx:
                p7, bc, name = i.split('\t')
                if p7 == p7_bc:
                    tmp_list.append(i)
            bc_list.append(tmp_list)
                    
        return [p7_list, bc_list]


    def index_to_dict(self, x, col_key=1, col_val=2):
        """Save index file to dict
        choose the column of <key> and <value>

        input file:
        <p7>[TAB]<bc>[TAB]<name>

        """
        assert isinstance(col_key, int)
        assert isinstance(col_val, int)

        # python starts with 0
        col_key -= 1
        col_val -= 1

        dd = {}
        with open(x, 'rt') as fi:
            for line in fi:
                if ' ' in line:
                    sep = ' '
                elif '\t' in line:
                    sep = '\t'
                else:
                    raise Exception('unknown file format, expect <TAB> separated file')
                ps = line.strip().split(sep)
                #
                if not len(ps) == 3:
                    raise Exception('unknown file format, expected <p7>[TAB]<barcode>[TAB]<name> - %s' % line)
                ps[col_key] = ps[col_key].upper()
                if ps[col_key] == 'NULL':
                    continue
                if ps[col_key] in dd:
                    raise Exception('Duplicates detected in column: %s' % str(col_key + 1))
                dd[ps[col_key]] = ps[col_val]
        return dd


    def p7_extractor(self, x, p7_len=0):
        """Extract the p7-index from fastq file
        eg:
        p7 in the tail of name: 
        @ST-E00310:586:HJH7JCCXY:7:1101:27275:1555 1:N:0:NCAACAAT

        the first 6-nt was the p7-index
        """
        assert isinstance(x, list)
        assert isinstance(p7_len, int)
        x_name, x_seq, _, x_qual = x
        p7_query = x_name.split(':')[-1]
        if not p7_len in list(range(4, 8)):
            raise Exception('Illegal <p7_len>, should be [4, 8] in length: %s' % p7_len)
        p7_query = p7_query[:p7_len]
        return [(x_name, x_seq, '+', x_qual), p7_query]


    def barcode_extractor(self, x, cut=False, bc_len=0, bc_n_left=0, bc_n_right=0):
        """Extract the barcode, p7-index from sequence
        : arg x, a list contins fastq elements [name, seq, +, qual]
        : arg cut, cut barcode from fastq sequence
        : arg bc_len, the length of barcode
        : arg bc_n_left, number of bases at the left of barcode
        : arg bc_n_right, number of bases at the right of barcode
        """
        assert isinstance(x, list)
        x_name, x_seq, _, x_qual = x
        if bc_len > 0:
            _x = (bc_n_left + bc_len)
            _y = (bc_n_left + bc_len + bc_n_right)
            bc_query = x_seq[bc_n_left:_x]
            bc_random = str(x_seq[:bc_n_left]) + str(x_seq[_x:_y])
            if cut is True:
                _name = s_description.split(' ')
                _name.insert(1, bc_random)
                x_name = ' '.join(_name)
                x_se = x_seq[_y:]
                x_qual = x_qual[_y:]
        else:
            raise Exception('Illegal <bc_len>, should be positive integer: %s' % bc_len)
        return [(x_name, x_seq, '+', x_qual), bc_query]


    def str_validator(self, x, x_dict, mm=0):
        """Validate the string in dictionary keys
        : arg x, a str, indicate the barcode sequence, or p7 index sequence
        : arg x_dict, a dict, the known
        : arg mm, maximum mismatch allowed

        return tag:
        x, the input string
        'undemx' 
        # bug:
        !!!multiple hits!!!!
        """
        assert isinstance(x, str)
        assert isinstance(x_dict, dict)
        assert isinstance(mm, int)

        # create tag
        for i in  list(x_dict.keys()):
            m = self.str_mismatch(x, i)
            if m <= mm:
                tag = i
                break
            else:
                tag = 'undemx'

        return tag


    def p7_demx_pe(self):
        """Demultiplex PE reads by P7 index
        : arg fq1, read1 of PE reads
        : arg fq2, read2 of PE reads
        : arg sample_info, the file saving P7-index information
        : arg path_out, the directory to save demultiplexed fastq files
        : arg mm, maximum mismatch allowed
        """
        args = self.args
        fq1 = args['fq1']
        fq2 = args['fq2']
        sample_info = self.sample_info
        path_out = self.path_out
        mm = args['n_mismatch']

        p7_dict = self.index_to_dict(sample_info, col_key=1, col_val=3)
        p7_len = self.key_len(p7_dict)
        if p7_len is None:
            raise ValueError('Empty sampleinfo or the length of P7-index are not consistent: %s' % sample_info)
        p7_dict['undemx'] = 'undemx' # undemx file
        # p7_dict['multi'] = 'multi' # barcode match multi hits

        ####################
        ## prepare read12 ##
        ####################
        r1_suffix = '_1.fq.gz'
        r2_suffix = '_2.fq.gz'

        ####################
        ## prepare writer ##
        ####################
        p7_count = {}
        p7_writer = {}
        for idx in p7_dict:
            p7_count[idx] = {}
            p7_count[idx]['count'] = 0
            p7_count[idx]['name'] = p7_dict[idx]
            idx_r1_name = p7_dict[idx] + r1_suffix
            idx_r2_name = p7_dict[idx] + r2_suffix
            p7_writer[idx] = [gzip.open(os.path.join(path_out, idx_r1_name), 'wt'),
                              gzip.open(os.path.join(path_out, idx_r2_name), 'wt'),]

        ####################
        ## prepare reader ##
        ####################
        fq_reader1 = gzip.open if is_gz(fq1) else open
        fq_reader2 = gzip.open if is_gz(fq2) else open
        with fq_reader1(fq1, 'rt') as f1, fq_reader2(fq2, 'rt') as f2:
            while True:
                try:
                    s1 = [next(f1), next(f1), next(f1), next(f1)]
                    s2 = [next(f2), next(f2), next(f2), next(f2)]
                    s1 = [i.strip() for i in s1]
                    s2 = [i.strip() for i in s2]
                    # check read names from read1
                    if not s1[0].split(' ')[0] == s2[0].split(' ')[0]:
                        raise ValueError('name of read1/2 are not consistent: \
                            %s \n %s' % (s1[0], s2[0]))
                    # extract P7 sequence
                    s1_new, p7_query = self.p7_extractor(s1, p7_len=p7_len)
                    tag = self.str_validator(p7_query, p7_dict, mm=mm)
                    # write fastq
                    p7_writer[tag][0].write('\n'.join(s1_new) + '\n')
                    p7_writer[tag][1].write('\n'.join(s2) + '\n') # not processed
                    p7_count[tag]['count'] += 1
                except StopIteration:
                    break
        # close writers
        for idx in p7_writer: 
            p7_writer[idx][0].close() # close writers
            p7_writer[idx][1].close() # close writers
        
        ## save log
        lib = os.path.join(path_out, 'run_args.lib')
        args['sub-command'] = 'p7_demx_pe'
        Json_file(args).json_writer(to=lib)

        report_file = os.path.join(path_out, 'report_demx.json')
        Json_file(p7_count).json_writer(report_file)


    def p7_demx_se(self):
        """Demultiplex PE reads by P7 index
        : arg fq1, read1 of PE reads
        : arg sample_info, the file saving P7-index information
        : arg path_out, the directory to save demultiplexed fastq files
        : arg mm, maximum mismatch allowed
        """
        args = self.args
        fq1 = args['fq1']
        sample_info = self.sample_info
        path_out = self.path_out
        mm = args['n_mismatch']

        p7_dict = self.index_to_dict(sample_info, col_key=1, col_val=3)
        p7_len = self.key_len(p7_dict)
        if p7_len is None:
            raise ValueError('Empty sampleinfo or the length of P7-index are not consistent: %s' % sample_info)
        p7_dict['undemx'] = 'undemx' # undemx file
        # p7_dict['multi'] = 'multi' # barcode match multi hits

        ####################
        ## prepare writer ##
        ####################
        p7_count = {}
        p7_writer = {}
        for idx in p7_dict:
            p7_count[idx] = {}
            p7_count[idx]['count'] = 0
            p7_count[idx]['name'] = p7_dict[idx]
            idx_name = p7_dict[idx] + '.fq.gz'
            p7_writer[idx] = open(os.path.join(path_out, idx_name), 'wt')

        ####################
        ## prepare reader ##
        ####################
        fq_reader = gzip.open if is_gz(fq1) else open
        with fq_reader(fq1, 'rt') as f1:
            while True:
                try:
                    s1 = [next(f1), next(f1), next(f1), next(f1)]
                    s1 = [i.strip() for i in s1]
                    # extract P7 sequence
                    s1_new, p7_query = self.p7_extractor(s1, p7_len=p7_len)
                    tag = self.str_validator(p7_query, p7_dict, mm=mm)
                    p7_writer[tag].write('\n'.join(s1_new) + '\n')
                    p7_count[tag]['count'] += 1
                except StopIteration:
                    break
        # close writers
        for idx in p7_writer: 
            p7_writer[idx].close() # close writers

        ## save log
        lib = os.path.join(path_out, 'run_args.lib')
        args['sub-command'] = 'p7_demx_se'
        Json_file(args).json_writer(to=lib)

        report_file = os.path.join(path_out, 'report_demx.json')
        Json_file(p7_count).json_writer(report_file)


    def bc_demx_pe(self):
        """Demultiplex PE reads
        : arg fq1, read1 of PE reads
        : arg fq2, read2 of PE reads
        : arg sample_info, the file saving P7-index information
        : arg path_out, the directory to save demultiplexed fastq files
        : arg mm, maximum mismatch allowed
        : arg bc_in_read, barcode in read1 or read2, default: 1
        : arg n_left, number of bases in the left of barcode
        : arg n_right, number of bases in the right of barcode
        : arg cut, cut the barcode and ramdomer from sequence
        """
        args = self.args
        fq1 = args['fq1']
        fq2 = args['fq2']
        sample_info = self.sample_info
        path_out = self.path_out
        mm = args['n_mismatch']
        bc_in_read = args['bc_in_read']
        n_left = args['n_left']
        n_right = args['n_right']
        cut = args['cut']

        bc_dict = self.index_to_dict(sample_info, col_key=2, col_val=3)
        bc_len = self.key_len(bc_dict)
        if bc_len is None:
            raise Exception('Empty sampleinfo or length of barcodes are not consistent: %s' % sample_info)
        bc_dict['undemx'] = 'undemx' # add undemx
        # bc_dict['multi'] = 'multi' # barcode match multi hits

        ####################
        ## prepare read12 ##
        ####################
        r1_suffix = '_1.fq.gz'
        r2_suffix = '_2.fq.gz'
        if bc_in_read == 2:
            (fq1, fq2) = (fq2, fq1) # switch read 1/2
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
            bc_r1_name = bc_dict[bc] + r1_suffix
            bc_r2_name = bc_dict[bc] + r2_suffix
            bc_writer[bc] = [gzip.open(os.path.join(path_out, bc_r1_name), 'wt'),
                             gzip.open(os.path.join(path_out, bc_r2_name), 'wt'),]

        ####################
        ## prepare reader ##
        ####################
        fq_reader1 = gzip.open if is_gz(fq1) else open
        fq_reader2 = gzip.open if is_gz(fq2) else open
        with fq_reader1(fq1, 'rt') as f1, fq_reader2(fq2, 'rt') as f2:
            while True:
                try:
                    s1 = [next(f1), next(f1), next(f1), next(f1)]
                    s2 = [next(f2), next(f2), next(f2), next(f2)]
                    s1 = [i.strip() for i in s1]
                    s2 = [i.strip() for i in s2]
                    # check read names from read1
                    if not s1[0].split(' ')[0] == s2[0].split(' ')[0]:
                        raise ValueError('name of read1/2 are not consistent: \
                            %s \n %s' % (s1[0], s2[0]))
                    # default bc in s1
                    s1_new, bc_query = self.barcode_extractor(s1, cut=False, 
                        bc_len=bc_len, bc_n_left=n_left, bc_n_right=n_right)
                    tag = self.str_validator(bc_query, bc_dict, mm=mm)
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

        ## save log
        lib = os.path.join(path_out, 'run_args.lib')
        args['sub-command'] = 'bc_demx_pe'
        Json_file(args).json_writer(to=lib)

        report_file = os.path.join(path_out, 'report_demx.json')
        Json_file(bc_count).json_writer(report_file)


# fq1 = 'data/raw_data/demo_01m_1.fq.gz'
# fq2 = 'data/raw_data/demo_01m_2.fq.gz'
# sample_info = 'data/sample_info/samplesheet.txt'
# path_out = 'results/demo'

# Demx(fq1=fq1, fq2=fq2, sample_info=sample_info, path_out=path_out,
#     n_mismatch=2).p7_demx_pe()

# Demx(fq1=fq1, fq2=fq2, sample_info=sample_info, path_out=path_out,
#     n_mismatch=2).p7_demx_se()

