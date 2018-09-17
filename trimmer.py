#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trimming reads
## Required
1. 3' adapter, (default: TruSeq, optional)
2. low-quality (q=20) at 3' end
3. trim-n
4. remove sequences do not contain adapters (--rm-untrim)
## Optional
1. trim N-bases at either ends
2. trim sliding window of 3'-adapter
3. trim adapter by multiple times (--times=N)
4. remove PCR duplicates (--rm-PCR-dup)
"""

import os
import sys
import re
import shlex
import subprocess
import logging

import goldclip
from utils_parser import *
from helper import *

class Trimmer(object):
    """Trim fastq files by 3' adapter, further quality control
    ## Required
    1. 3' adapter, (default: TruSeq, optional)
    2. low-quality (q=20) at 3' end
    3. trim-n
    4. remove sequences do not contain adapters (--rm-untrim)
    ## Optional
    1. trim N-bases at either ends
    2. trim sliding window of 3'-adapter
    3. trim adapter by multiple times (--times=N)
    4. remove PCR duplicates (--rm-PCR-dup)

    ## Basic trimming ##

    ## cutadapt
    input: fq, adapter3, path_out, len_min, qual_min, err_rate, overlap, 
           adapter_sliding, double_trim, multi_cores, rm_untrim,
           cut_before_trim, overwrite
    output: [Trimmer class], [args], fq, untrim_fq, log

    ## Further trimming ##
    
    ## Trim N-bases after cutadapt
    input: [Trimmer class], fq, cut_after_trim
    output: [Trimmer class], fq, log

    ## Remove PCR-duplicates
    input: [Trimmer class], fq, rm_dup
    output: [Trimmer class], fq, log

    """

    def __init__(self, fq, adapter3, path_out=None, adapter5=None, len_min=15,
                 adapter_sliding=False, trim_times=1, double_trim=False, 
                 qual_min=20, err_rate=0.1, overlap=3, multi_cores=1, read12=1, 
                 rm_untrim=False, rm_dup=False, cut_before_trim='0', 
                 cut_after_trim='0', overwrite=False):
        """Process parameters"""
        self.fq = fq
        self.adapter3 = adapter3
        self.path_out = path_out
        self.adapter5 = adapter5
        self.len_min = len_min
        self.adapter_sliding = adapter_sliding
        self.trim_times = trim_times
        self.double_trim = double_trim
        self.qual_min =qual_min
        self.err_rate = err_rate
        self.overlap = overlap
        self.multi_cores = multi_cores
        self.read12 = read12
        self.rm_untrim = rm_untrim
        self.rm_dup = rm_dup
        self.cut_before_trim = cut_before_trim
        self.cut_after_trim = cut_after_trim
        self.overwrite = overwrite
        ## wrapper parameters
        self.args = self.get_args()
        # fq_prefix, fq_clean, fq_log, fq_untrim = self.init_dir()
        # if rm_dup is True and self.cut_after_trim == '0':
        #     self.fq_clean = self.run()
        #     self.fq_clean = self.rm_duplicate()
        #     self.fq_clean = self.trim_ends()
        # elif rm_dup is True:
        #     self.fq_clean = self.run()
        #     self.fq_clean = self.rm_duplicate()
        # elif not self.cut_after_trim == '0':
        #     self.fq_clean = self.run()
        #     self.fq_clean = self.trim_ends()
        # else:
        #     self.fq_clean = self.run()
        #     # self.fq_clean = fq_clean
        


    def get_args(self):
        """Parse parameters"""
        args = {'fq' : self.fq, 
                'adapter3' : self.adapter3, 
                'path_out' : self.path_out, 
                'adapter5' : self.adapter5,
                'len_min' : self.len_min,
                'adapter_sliding' : self.adapter_sliding,
                'trim_times' : self.trim_times,
                'double_trim' : self.double_trim,
                'qual_min' : self.qual_min,
                'err_rate' : self.err_rate,
                'overlap' : self.overlap,
                'multi_cores' : self.multi_cores,
                'read12' : self.read12,
                'rm_untrim' : self.rm_untrim,
                'rm_dup' : self.rm_dup,
                'cut_before_trim' : self.cut_before_trim,
                'cut_after_trim' : self.cut_after_trim,
                'overwrite' : self.overwrite
                }
        # check
        assert isinstance(args['fq'], str)
        assert os.path.exists(args['fq'])
        assert isinstance(args['adapter3'], str)
        assert is_path(args['path_out'])
        # assert isinstance(args['adatper5'], str)
        assert isinstance(args['len_min'], int)
        assert isinstance(args['adapter_sliding'], bool)
        assert isinstance(args['double_trim'], bool)
        assert isinstance(args['qual_min'], int)
        assert isinstance(args['err_rate'], float)
        assert isinstance(args['overlap'], int)
        assert isinstance(args['multi_cores'], int)
        assert isinstance(args['read12'], int)
        assert isinstance(args['rm_untrim'], bool)
        assert isinstance(args['rm_dup'], bool)
        assert isinstance(args['cut_before_trim'], str)
        assert isinstance(args['cut_after_trim'], str)
        assert isinstance(args['overwrite'], bool)
        return args



    def init_dir(self, fq_in=None):
        """Prepare directory, filename for each file"""
        args = self.args

        if fq_in is None:
            fq_in = args['fq']
        fq_prefix = file_prefix(fq_in)[0]
        # fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
        fq_type = seq_type(fq_in)

        fq_clean = os.path.join(args['path_out'], '%s.clean.fastq' % fq_prefix)
        fq_log = os.path.join(args['path_out'], '%s.cutadapt.log' % fq_prefix)
        fq_untrim = os.path.join(args['path_out'], '%s.untrim.fastq' % fq_prefix)

        return [fq_prefix, fq_clean, fq_log, fq_untrim]



    def ad_chopper(self, ad=None, step=2, window=15):
        """Chop adapter by given length and step, to create a series of adapters"""
        args = self.args
        if ad is None:
            ad = args['adapter3']
        assert isinstance(ad, str)
        assert isinstance(step, int)
        assert isinstance(window, int)
        p = []
        if len(ad) < window:
            p.append(ad)
        else:
            for i in range(int(len(ad) / step)):
                a = i * step
                b = a + window
                if b > len(ad):
                    continue
                p.append(ad[a:b])
        return p



    def wrap_args(self):
        """Wrapper args for cutadapt"""
        args = self.args

        fq_prefix, fq_clean, fq_log, fq_untrim = self.init_dir()

        cutadapt_arg = 'cutadapt -m %s -q %s --overlap=%s --error-rate=%s \
                        --trim-n --max-n=0.1 --times=%s' % (args['len_min'],
                        args['qual_min'], args['overlap'], args['err_rate'],
                        args['trim_times'])
 
        # adapter
        ad3_list = [args['adapter3'], ]
        if args['adapter_sliding'] is True:
            ad3_list = self.ad_chopper(args['adapter3'])
        arg_ad3 = ' '.join(['-a %s' % i for i in ad3_list])

        # un_trim
        if args['rm_untrim'] is True:
            arg_untrim = '--untrimmed-output=%s --cores=1' % fq_untrim
        else:
            arg_untrim = '--cores=%s' % args['multi_cores']

        # cut before
        if args['cut_before_trim'] == '0':
            arg_cut = ''
        else:
            arg_cut = '--cut %s' % args['cut_before_trim']

        arg_line = ' '.join([cutadapt_arg, arg_ad3, arg_untrim, arg_cut, args['fq']])

        return arg_line



    def cutadapt_cut(self, s, cut_para=True):
        """
        recognize para: cut for cutadapt
        eg: cut=6, cut=-3, cut=6,-3
        """
        if ',' in s:
            n = s.split(',')
            if len(n) > 2:
                raise ValueError('illegal ad_cut: %s' % s)
            else:
                c1, c2 = (int(n[0]), int(n[1]))
                if c1 < 0 or c2 > 0:
                    raise ValueError('illegal ad_cut: %s' % s)
            if cut_para:
                c_para = '--cut %s --cut %s' % (c1, c2)
            else:
                c_para = [c1, c2]
        else:
            if cut_para:
                c_para = '--cut %s' % int(s)
            else:
                c_para = [int(s), ]
        return c_para




    def _is_non_empty(self, fn):
        """file if not empyt"""
        if os.path.getsize(fn) > 0:
            return True
        else:
            return False


    def _is_file(self, fn):
        """file exists"""
        if os.path.isfile(fn):
            return True
        else:
            return False


    def _tmp(self):
        """Create a temp file"""
        args = self.args
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.fq',
                                            dir=args['path_out'],
                                            delete=False)
        return tmpfn.name



    def saveas(self, _out=None):
        """Save output fastq file"""
        if _out is None:
            _out = self._tmp()

        self.fq_clean
        assert os.path.exists(self.fq_clean)
        os.rename(self.fq_clean, _out)
        return _out



    def trim_se(self):
        """Trimming reads using cutadapt"""
        args = self.args
        fq_prefix, fq_clean, fq_log, fq_untrim = self.init_dir()
        logging.info('trimming adapter: %s' % fq_prefix)

        arg_line = self.wrap_args()
        with open(fq_clean, 'wt') as ff, open(fq_log, 'wt') as fg:
            p1 = subprocess.run(shlex.split(arg_line), stdout=ff, stderr=fg)
        Cutadapt_log(fq_log).saveas()
        return fq_clean



    def rm_duplicate(self, fq_in=None):
        """Remove possible PCR duplicates"""
        pkg_dir = os.path.split(goldclip.__file__)[0]
        fa2fq = os.path.join(pkg_dir, 'bin', 'fasta_to_fastq.pl')
        args = self.args
        if fq_in is None:
            fq_in = args['fq']
        fq_prefix, fq_clean, fq_log, fq_untrim = self.init_dir(fq_in)
        logging.info('removing PCR duplicates')

        if args['rm_dup'] is True:
            fq_clean = os.path.splitext(fq_in)[0] + '.nodup.fastq'
            if os.path.exists(fq_clean) and args['overwrite'] is False:
                logging.info('file exists, skip rm_dup: %s' % fq_prefix)
            else:
                # create new clean fastq file
                with open(fq_clean, 'wt') as ff:
                    c1 = 'fastx_collapser -Q33 -i %s' % fq_in
                    c2 = 'perl %s -' % fa2fq
                    p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE)
                    p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, stdout=ff)
                    px = p2.communicate()
        else:
            logging.info('Checkout if --rm_dup True, skipped this round')
        return fq_clean



    def trim_ends(self, fq_in=None):
        """Trim N-bases at either tail of fastq after cutadapt"""
        args = self.args
        if fq_in is None:
            fq_in = args['fq']
        fq_prefix, fq_clean, fq_log, fq_untrim = self.init_dir(fq_in)
        logging.info('trimming ends')

        if args['cut_after_trim'] == '0':
            logging.info('Checkout if --cut-after-trim 0, skipped this round')
        else:
            # move clean fq to temp
            fq_clean = os.path.splitext(fq_in)[0] + '.cut.fastq'
            trim_args = self.cutadapt_cut(args['cut_after_trim'], False)
            if len(trim_args) == 2:
                trim_5, trim_3 = trim_args
            elif len(trim_args) == 1:
                trim_5 = 0 if(int(trim_args[0]) < 0) else trim_args[0]
                trim_3 = 0 if(int(trim_args[0]) > 0) else trim_args[0]
            else:
                trim_5 = trim_3 = 0
            # run trimming
            if os.path.exists(fq_clean) and args['overwrite'] is False:
                logging.info('file exists, skip trimming ends: %s' % fq_prefix)
            else:
                with open(fq_in, 'rt') as fi, open(fq_clean, 'wt') as ff:
                    while True:
                        try:
                            fq_id, fq_seq, fq_plus, fq_qual = [next(fi).strip(),
                                                               next(fi).strip(),
                                                               next(fi).strip(),
                                                               next(fi).strip()]
                            if len(fq_seq) < args['len_min'] + abs(trim_5) + abs(trim_3): 
                                continue # skip short reads
                            fq_seq = fq_seq[trim_5:trim_3] if(trim_3 < 0) else fq_seq[trim_5:]
                            fq_qual = fq_qual[trim_5:trim_3] if(trim_3 < 0) else fq_qual[trim_5:]
                            ff.write('\n'.join([fq_id, fq_seq, fq_plus, fq_qual]) + '\n')
                        except StopIteration:
                            break
        return fq_clean



    def run(self):
        """Trimming all"""
        args = self.args
        fq_clean = self.trim_se()
        if not args['cut_after_trim'] == '0' and args['rm_dup'] is True:
            fq1 = self.rm_duplicate(fq_clean)
            fq2 = self.trim_ends(fq1)
            os.remove(fq_clean)
            os.remove(fq1)
            fq_clean = fq2
        elif not args['cut_after_trim'] == '0':
            fq1 = self.trim_ends(fq_clean)
            os.remove(fq_clean)
            fq_clean = fq1
        elif args['rm_dup'] is True:
            fq1 = self.rm_duplicate(fq_clean)
            os.remove(fq_clean)
            fq_clean = fq1
        else:
            pass
            # raise ValueError('unknown --rm-dup %s, --cut-after-trim %s' % (
            #                   args['rm_dup'], args['cut_after_trim']))
        return fq_clean

