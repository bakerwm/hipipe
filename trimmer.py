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


2018-12-13 to-do

1. cutadapt log parser for PE trimming

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
from arguments import args_init


# class Trimmer(object):
#     """Trim fastq files by 3' adapter, further quality control
#     ## Required
#     1. 3' adapter, (default: TruSeq, optional)
#     2. low-quality (q=20) at 3' end
#     3. trim-n
#     4. remove sequences do not contain adapters (--rm-untrim)
#     ## Optional
#     1. trim N-bases at either ends
#     2. trim sliding window of 3'-adapter
#     3. trim adapter by multiple times (--times=N)
#     4. remove PCR duplicates (--rm-PCR-dup)

#     ## Basic trimming ##

#     ## cutadapt
#     input: fq, adapter3, path_out, len_min, qual_min, err_rate, overlap, 
#            adapter_sliding, double_trim, multi_cores, rm_untrim,
#            cut_before_trim, overwrite
#     output: [Trimmer class], [args], fq, untrim_fq, log

#     ## Further trimming ##
    
#     ## Trim N-bases after cutadapt
#     input: [Trimmer class], fq, cut_after_trim
#     output: [Trimmer class], fq, log

#     ## Remove PCR-duplicates
#     input: [Trimmer class], fq, rm_dup
#     output: [Trimmer class], fq, log

#     """


#     def __init__(self, fq1, adapter3, path_out=None, len_min=15, **kwargs):
#         """Parsing the parameters for reads trimming
#         support both SE and PE reads
#         """
#         ## SE options
#         self.fq1 = fq1
#         self.adapter3 = adapter3
#         self.path_out = path_out
#         self.len_min = len_min
#         self.adapter5 = kwargs.get('adapter5', None)
#         self.read12 = kwargs.get('read12', 1)

#         ## Paired-end options
#         self.fq2 = kwargs.get('fq2', None)
#         self.AD3 = kwargs.get('AD3', None)
#         self.AD5 = kwargs.get('AD5', None)

#         ## global options
#         self.qual_min = kwargs.get('qual_min', 20)
#         self.err_rate = kwargs.get('err_rate', 0.1)
#         self.overlap = kwargs.get('overlap', 3)
#         self.rm_untrim = kwargs.get('rm_untrim', False)
#         self.threads = kwargs.get('threads', 1)
#         self.overwrite = kwargs.get('overwrite', False)
#         self.keep_name = kwargs.get('keep_name', True) # do not change name

#         ## extra options
#         self.adapter_sliding = kwargs.get('adapter_sliding', False)
#         self.trim_times = kwargs.get('trim_times', 1)
#         self.double_trim = kwargs.get('double_trim', False)
#         self.rm_dup = kwargs.get('rm_dup', False)
#         self.cut_before_trim = kwargs.get('cut_before_trim', 0)
#         self.cut_after_trim = kwargs.get('cut_after_trim', 0)
#         self.trim_to_length = kwargs.get('trim_to_length', 0)

#         ## validate options
#         is_path(self.path_out)
#         if self.trim_to_length > 0 and self.trim_to_length < self.len_min:
#             raise ValueError('[fatal] --trim-to-length [%s] shorter than -m [%s]') % (self.trim_to_length, self.len_min)



#     def init_dir(self, fq_in=None):
#         """Prepare directory, filename for each file"""
#         # args = self.args
#         if fq_in is None:
#             fq_in = self.fq1

#         # SE mode
#         if self.fq2 is None:
#             fq_prefix = file_prefix(fq_in)[0]
#             # fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
#             fq_type = seq_type(fq_in)
#             if self.keep_name:
#                 fq_clean = os.path.join(self.path_out, '%s.fastq' % fq_prefix)
#             else:
#                 fq_clean = os.path.join(self.path_out, '%s.clean.fastq' % fq_prefix)
#             fq_log = os.path.join(self.path_out, '%s.cutadapt.log' % fq_prefix)
#             fq_untrim = os.path.join(self.path_out, '%s.untrim.fastq' % fq_prefix)
#             return [fq_prefix, fq_clean, fq_log, fq_untrim]
#         # PE mode
#         elif os.path.exists(self.fq2):
#             fq1_name = file_prefix(self.fq1)[0]
#             fq2_name = file_prefix(self.fq2)[0]
#             if self.keep_name:
#                 fq1_clean = os.path.join(self.path_out, '%s.fastq' % fq1_name)
#                 fq2_clean = os.path.join(self.path_out, '%s.fastq' % fq2_name)
#             else:
#                 fq1_clean = os.path.join(self.path_out, '%s.clean.fastq' % fq1_name)
#                 fq2_clean = os.path.join(self.path_out, '%s.clean.fastq' % fq2_name)
#             fq_log = os.path.join(self.path_out, '%s.cutadapt.log' % fq1_name)
#             fq1_untrim = os.path.join(self.path_out, '%s.untrim.fastq' % fq1_name)
#             fq2_untrim = os.path.join(self.path_out, '%s.untrim.fastq' % fq2_name)
#             return [fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim]
#         else:
#             raise ValueError('checkout -fq2 option')


#     def ad_chopper(self, ad=None, step=2, window=15):
#         """Chop adapter by given length and step, to create a series of adapters"""
#         args = self.args
#         if ad is None:
#             ad = args['adapter3']
#         assert isinstance(ad, str)
#         assert isinstance(step, int)
#         assert isinstance(window, int)
#         p = []
#         if len(ad) < window:
#             p.append(ad)
#         else:
#             for i in range(int(len(ad) / step)):
#                 a = i * step
#                 b = a + window
#                 if b > len(ad):
#                     continue
#                 p.append(ad[a:b])
#         return p


#     def get_cutadapt_cmd(self):
#         """Parse arguments, create cutadapt command line"""
#         ## command
#         cut_exe = which('cutadapt')

#         ## determine the minimum length of reads
#         ## consider cut-after-trim
#         trim_args = self.cutadapt_cut(self.cut_after_trim, False) # return the numbers
#         if len(trim_args) == 2:
#             trim_5, trim_3 = trim_args[:2]
#         elif len(trim_args) == 1:
#             trim_5 = 0 if(int(trim_args[0]) < 0) else trim_args[0]
#             trim_3 = 0 if(int(trim_args[0]) > 0) else trim_args[0]
#         else:
#             trim_5 = trim_3 = 0
#         len_min = self.len_min + abs(trim_5) + abs(trim_3)
#         ## basic command
#         cut_args = '%s -m %s -q %s --overlap=%s --error-rate=%s --trim-n \
#             --max-n=0.1 --times=%s' % (cut_exe, len_min, self.qual_min, 
#                 self.overlap, self.err_rate, self.trim_times)

#         ## adapter3
#         ad3_list = [self.adapter3, ]
#         if self.adapter_sliding:
#             ad3_list = self.ad_chopper(self.adapter3)
#         ad3_arg = ' '.join(['-a %s' % i for i in ad3_list])

#         ## SE mode
#         if self.fq2 is None:
#             fq_prefix, fq_clean, fq_log, fq_untrim = self.init_dir()

#             ## untrim
#             if self.rm_untrim:
#                 untrim_arg = '--untrimmed-output=%s --cores=1 -o %s' % (fq_untrim, fq_clean)
#             else:
#                 untrim_arg = '--cores=%s -o %s' % (self.threads, fq_clean)

#             ## cut before
#             if self.cut_before_trim > '0':
#                 cut_arg ='--cut %s'  % self.cut_before_trim
#             else:
#                 cut_arg = ''

#             arg_cmd = ' '.join([cut_args, ad3_arg, untrim_arg, cut_arg, self.fq1])

#         ## PE mode
#         elif os.path.exists(self.fq2):
#             fq_prefix, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim = self.init_dir()

#             AD3_list = [self.AD3, ]
#             if self.adapter_sliding:
#                 AD3_list = self.ad_chopper(self.AD3)
#             AD3_arg = ' '.join(['-A %s' % i for i in AD3_list])

#             ## untrim
#             if self.rm_untrim:
#                 untrim_arg = '--untrimmed-output=%s --untrimmed-paired-output=%s \
#                     --cores=1 -o %s --paired-output=%s' % (fq1_untrim, fq2_untrim,
#                         fq1_clean, fq2_clean)
#             else:
#                 untrim_arg = '--cores=%s -o %s --paired-output=%s' % (self.threads,
#                         fq1_clean, fq2_clean)

#             ## cut before
#             if self.cut_before_trim > '0':
#                 cut_arg ='--cut %s'  % self.cut_before_trim
#             else:
#                 cut_arg = ''

#             ## merge
#             arg_cmd = ' '.join([cut_args, ad3_arg, AD3_arg, untrim_arg, cut_arg, self.fq1, self.fq2])
#         else:
#             raise ValueError('--fq2 value illegal')

#         return arg_cmd


#     def cutadapt_cut(self, s, cut_para=True):
#         """Cut N-bases, parse --cut argument
#         recognize para: cut for cutadapt
#         eg: cut=6, cut=-3, cut=6,-3
#         """
#         if ',' in s:
#             n = s.split(',')
#             if len(n) > 2:
#                 raise ValueError('illegal ad_cut: %s' % s)
#             else:
#                 c1, c2 = (int(n[0]), int(n[1]))
#                 if c1 < 0 or c2 > 0:
#                     raise ValueError('illegal ad_cut: %s' % s)
#             if cut_para:
#                 c_para = '--cut %s --cut %s' % (c1, c2)
#             else:
#                 c_para = [c1, c2]
#         else:
#             if cut_para:
#                 c_para = '--cut %s' % int(s)
#             else:
#                 c_para = [int(s), ]
#         return c_para


#     def _is_non_empty(self, fn):
#         """file if not empyt"""
#         if os.path.getsize(fn) > 0:
#             return True
#         else:
#             return False


#     def _is_file(self, fn):
#         """file exists"""
#         if os.path.isfile(fn):
#             return True
#         else:
#             return False


#     def _tmp(self):
#         """Create a temp file"""
#         # args = self.args
#         tmpfn = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.fastq',
#                                             dir=self.path_out,
#                                             delete=False)
#         return tmpfn.name


#     def saveas(self, _out=None):
#         """Save output fastq file"""
#         ## SE mode
#         if self.fq2 is None:
#             fq_clean = self.init_dir()[1]

#             if _out is None:
#                 _out = fq_clean

#             if os.path.exists(fq_clean):
#                 os.rename(fq_clean, _out)
#             else:
#                 raise ValueError('need run cutadapt first')

#         elif os.path.exists(self.fq2):
#             # [fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim]
#             fq1_clean, fq2_clean = self.init_dir()[1:3]
#             if _out is None:
#                 _out = [fq1_clean, fq2_clean]

#             if os.path.exists(fq1_clean) and os.path.exists(fq2_clean):
#                 os.rename(fq1_clean, _out[0])
#                 os.rename(fq2_clean, _out[1])
#             else:
#                 raise ValueError('need run cutadapt first')

#         else:
#             raise ValueError('--fq2 value illegal')

#         return _out


#     def trim_se(self):
#         """Trimming reads using cutadapt"""
#         fq_prefix, fq_clean, fq_log, fq_untrim = self.init_dir()
#         logging.info('trimming SE reads: %s' % fq_prefix)

#         arg_cmd = self.get_cutadapt_cmd()
#         if os.path.exists(fq_clean) and self.overwrite is False:
#             logging.info('file exists, cutadapt skipped: %s' % fq_prefix)
#         else:
#             with open(fq_log, 'wt') as fo:
#                 p1 = subprocess.run(shlex.split(arg_cmd), stdout=fo, stderr=fo)
#             Cutadapt_log(fq_log).saveas()
#         return fq_clean


#     def trim_pe(self):
#         """Trimming Paired-end reads using cutadapt
#         example: 
#         cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq

#         -A  see -a for read1
#         -G  see -g for read1
#         -B  see -b for read1
#         -u  see --cut for read1

#         -p  write the second read to FILE
#         --untrimmed-paired-output  write the second read in a pair to this FILE when no
#             adapter was found in the first read. together with --untrimmed-output
#         --too-short-paired-output  write the second read in a pair to this FILE if pair
#             is too short, together with --too-short-output
#         --too-long-paired-output    write the second read in a pair to this FILE if pair
#             is too long, together with --too-long-output
#         """
#         # [fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim]
#         fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim = self.init_dir()
#         logging.info('trimming PE reads: %s' % fq1_name)

#         arg_cmd = self.get_cutadapt_cmd()
#         if os.path.exists(fq1_clean) and os.path.exists(fq2_clean) and self.overwrite is False:
#             logging.info('file exists, cutadapt skipped: %s' % fq1_name)
#         else:
#             with open(fq_log, 'wt') as fo:
#                 p1 = subprocess.run(shlex.split(arg_cmd), stdout=fo, stderr=fo)
#             # Cutadapt_log(fq_log).saveas()
#         return [fq1_clean, fq2_clean]


#     def rm_duplicate(self, fq_in=None):
#         """Remove possible PCR duplicates
#         collapse reads before trimming
#         ****
#         only support SE reads (not PE reads)
#         """
#         pkg_dir = os.path.split(goldclip.__file__)[0]
#         fa2fq = os.path.join(pkg_dir, 'bin', 'fasta_to_fastq.pl')

#         if fq_in is None:
#             fq_in = self.fq1
#         fq_prefix, fq_clean, fq_log, fq_untrim = self.init_dir(fq_in)
#         logging.info('removing PCR duplicates')

#         if self.rm_dup:
#             fq_clean = os.path.splitext(fq_in)[0] + '.nodup.fastq'
#             if os.path.exists(fq_clean) and self.overwrite is False:
#                 logging.info('file exists, skip rm_dup: %s' % fq_prefix)
#             else:
#                 # create new clean fastq file
#                 with open(fq_clean, 'wt') as ff:
#                     c1 = 'fastx_collapser -Q33 -i %s' % fq_in
#                     c2 = 'perl %s -' % fa2fq
#                     p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE)
#                     p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, stdout=ff)
#                     px = p2.communicate()
#         else:
#             logging.info('Checkout if --rm_dup True, skipped this round')
#         return fq_clean


#     def trim_ends(self, fq_in=None, fq_out=None, rm_too_short=True):
#         """Trim N-bases at either tail of fastq after cutadapt
#         Using python to process fastq files
#         support SE reads only
#         """
#         # if self.fq2 is None:
#         #     fq_prefix, fq_clean, fq_log, fq_untrim = self.init_dir(fq_in)
#         # else:
#         #     fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim = self.init_dir()

#         if fq_in is None or fq_out is None:
#             raise ValueError('fq_in and fq_out are required')
#             # fq_in = self.fq1
#         logging.info('trimming ends')

#         if self.cut_after_trim == '0':
#             logging.info('Checkout --cut-after-trim, skipped...')
#         else:
#             # move clean fq to temp
#             # fq_clean = os.path.splitext(fq_in)[0] + '.cut.fastq'
#             trim_args = self.cutadapt_cut(self.cut_after_trim, False)
#             if len(trim_args) == 2:
#                 trim_5, trim_3 = trim_args[:2]
#             elif len(trim_args) == 1:
#                 trim_5 = 0 if(int(trim_args[0]) < 0) else trim_args[0]
#                 trim_3 = 0 if(int(trim_args[0]) > 0) else trim_args[0]
#             else:
#                 trim_5 = trim_3 = 0
#             # run trimming
#             # if os.path.exists(fq_clean) and args['overwrite'] is False:
#             #     logging.info('file exists, skip trimming ends: %s' % fq_prefix)
#             # else:
#             # fq_tmp = fq_clean + '.tmp'
#             if fq_in == fq_out:
#                 fq_tmp = fq_in + '.cut_ends.tmp'
#                 os.rename(fq_in, fq_tmp)
#                 fq_in = fq_tmp
#             with open(fq_in, 'rt') as fi, open(fq_out, 'wt') as ff:
#                 while True:
#                     try:
#                         fq_id, fq_seq, fq_plus, fq_qual = [next(fi).strip(),
#                                                            next(fi).strip(),
#                                                            next(fi).strip(),
#                                                            next(fi).strip()]
#                         if len(fq_seq) < self.len_min + abs(trim_5) + abs(trim_3):
#                             if rm_too_short:
#                                 continue # skip short reads
#                             else:
#                                 fq_seq = 'A'
#                                 fq_qual = 'J' # quality
#                         else:
#                             fq_seq = fq_seq[trim_5:trim_3] if(trim_3 < 0) else fq_seq[trim_5:]
#                             fq_qual = fq_qual[trim_5:trim_3] if(trim_3 < 0) else fq_qual[trim_5:]
#                         # trim to length
#                         if self.trim_to_length > self.len_min:
#                             n_right = min([self.trim_to_length, len(fq_seq)])
#                             fq_seq = fq_seq[:n_right]
#                             fq_qual = fq_qual[:n_right]
#                         ff.write('\n'.join([fq_id, fq_seq, fq_plus, fq_qual]) + '\n')
#                     except StopIteration:
#                         break

#         return fq_out


#     def trim_ends2(self, fq_in=None):
#         """Trim reads from right, save reads with maximum length"""
#         if fq_in is None:
#             fq_in = self.fq1
#         fq_prefix, fq_clean, fq_log, fq_untrim = self.init_dir(fq_in)

#         if self.trim_to_length > 0:
#             with open(fq_in, 'rt') as fi, open(fq_clean, 'wt') as ff:
#                 while True:
#                     try:
#                         fq_id, fq_seq, fq_plus, fq_qual = [next(fi).strip(),
#                                                            next(fi).strip(),
#                                                            next(fi).strip(),
#                                                            next(fi).strip()]
#                         n_right = min([self.trim_to_length, len(fq_seq)])
#                         fq_seq = fq_seq[:n_right]
#                         fq_qual = fq_qual[:n_right]
#                         ff.write('\n'.join([fq_id, fq_seq, fq_plus, fq_qual]) + '\n')
#                     except StopIteration:
#                         break
#         return fq_clean


#     def run(self):
#         """Run trimming"""
#         if self.fq2 is None:
#             fq_clean = self.trim_se()
#             if self.cut_after_trim == '0' and self.rm_dup:
#                 fq_tmp1 = self.rm_duplicate(fq_clean)
#                 fq_tmp2 = self.trim_ends(fq_tmp1, fq_tmp1)
#                 os.remove(fq_clean)
#                 os.remove(fq_tmp1)
#                 fq_clean = fq_tmp2
#             elif not self.cut_after_trim == '0':
#                 fq_tmp1 = self.trim_ends(fq_clean, fq_clean)
#                 # os.remove(fq_clean)
#                 fq_clean = fq_tmp1
#             elif self.rm_dup:
#                 fq_tmp1 = self.rm_duplicate(fq_clean)
#                 if self.trim_to_length > self.len_min:
#                     fq_tmp2 = self.trim_ends2(fq_tmp1)
#                     os.remove(fq_clean)
#                     os.remove(fq_tmp1)
#                     fq_clean = fq_tmp2
#                 else:
#                     os.remove(fq_tmp1)
#                     pass
#             else:
#                 if self.trim_to_length > self.len_min:
#                     fq1 = self.trim_ends2(fq_clean)
#                     # os.remove(fq_clean)
#                     fq_clean = fq1
#             return fq_clean
#         elif os.path.exists(self.fq2):
#             [fq1_clean, fq2_clean] = self.trim_pe()
#             if not self.cut_after_trim == '0':
#                 fq_tmp1 = self.trim_ends(fq1_clean, fq1_clean)
#                 fq_tmp2 = self.trim_ends(fq2_clean, fq2_clean)
#             return [fq1_clean, fq2_clean]

        
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
    input: fq, adapter3, path_out, len_min, qual_min, error_rate, overlap, 
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

    ## Collapse
    input: fq
    output: fq

    """

    def __init__(self, fq1, adapter3, path_out=None, len_min=15, **kwargs):
        """Parsing the parameters for reads trimming
        support both SE and PE reads
        """
        assert isinstance(fq1, str) # process one file
        args1 = args_init(kwargs) # default parameters
        args2 = {
            'fq1': fq1,
            'adapter3': adapter3,
            'path_out': path_out,
            'len_min': len_min}
        args = {**args1, **args2}

        ## validate options
        assert is_path(path_out)
        if args['trim_to_length'] > 0 and args['trim_to_length'] < len_min:
            raise ValueError('[fatal] --trim-to-length [%s] shorter than -m [%s]') % (args['trim_to_length'], len_min)

        self.kwargs = args # global 


    def trim_init(self, fq_in=None):
        """Prepare directory for trimming, filename for each file"""
        args = self.kwargs.copy()
        if fq_in is None:
            fq_in = args['fq1']

        # SE mode
        if args['fq2'] is None:
            fq_prefix = file_prefix(fq_in)[0]
            fq_type = seq_type(fq_in)
            if args['keep_name']:
                fq_clean = os.path.join(args['path_out'], '%s.fastq' % fq_prefix)
            else:
                fq_clean = os.path.join(args['path_out'], '%s.clean.fastq' % fq_prefix)
            fq_log = os.path.join(args['path_out'], '%s.cutadapt.log' % fq_prefix)
            fq_untrim = os.path.join(args['path_out'], '%s.untrim.fastq' % fq_prefix)
            return [fq_prefix, fq_clean, fq_log, fq_untrim]

        # PE mode
        elif os.path.exists(args['fq2']):
            fq1_name = file_prefix(args['fq1'])[0]
            fq2_name = file_prefix(args['fq2'])[0]
            if args['keep_name']:
                fq1_clean = os.path.join(args['path_out'], '%s.fastq' % fq1_name)
                fq2_clean = os.path.join(args['path_out'], '%s.fastq' % fq2_name)
            else:
                fq1_clean = os.path.join(args['path_out'], '%s.clean.fastq' % fq1_name)
                fq2_clean = os.path.join(args['path_out'], '%s.clean.fastq' % fq2_name)
            fq_log = os.path.join(args['path_out'], '%s.cutadapt.log' % fq1_name)
            fq1_untrim = os.path.join(args['path_out'], '%s.untrim.fastq' % fq1_name)
            fq2_untrim = os.path.join(args['path_out'], '%s.untrim.fastq' % fq2_name)
            return [fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim]
        else:
            raise ValueError('checkout -fq2 option')


    def adapter_chopper(self, ad=None, step=2, window=15):
        """Chop adapter by given length and step, to create a series of adapters"""
        args = self.kwargs.copy()
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


    def cut_parser(self, s, cut_para=True):
        """Cut N-bases, before adapter trimming
        parse the argument
        if cut_para is True,
        return the cut numbers in --cut {arg} version

        else
        return the numbers        
        """
        if ',' in s:
            n = s.split(',')
            if len(n) > 2:
                raise ValueError('illegal argument: %s' % s)
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


    def get_cutadapt_cmd(self):
        """Parse arguments, create cutadapt command line"""
        args = self.kwargs.copy()

        ## command
        cut_exe = which('cutadapt')

        ## determine the minimum length of reads
        ## consider cut-after-trim
        trim_args = self.cut_parser(args['cut_after_trim'], False) # return the numbers
        if len(trim_args) == 2:
            trim_5, trim_3 = trim_args[:2]
        elif len(trim_args) == 1:
            trim_5 = 0 if(int(trim_args[0]) < 0) else trim_args[0]
            trim_3 = 0 if(int(trim_args[0]) > 0) else trim_args[0]
        else:
            trim_5 = trim_3 = 0
        len_min = args['len_min'] + abs(trim_5) + abs(trim_3) # minimum lenght

        ## basic command
        cut_basic = '%s -m %s \
                    -q %s \
                    --overlap=%s \
                    --error-rate=%s \
                    --trim-n \
                    --max-n=0.1 \
                    --times=%s' % (
                    cut_exe, 
                    len_min, 
                    args['qual_min'], 
                    args['overlap'], 
                    args['error_rate'], 
                    args['trim_times'])

        ## adapter3
        ad3_list = [args['adapter3'], ]
        ## add short version of adapter3
        if args['adapter_sliding']:
            ad3_list = self.adapter_chopper(args['adapter3'])
        ad3_arg = ' '.join(['-a %s' % i for i in ad3_list])

        ## SE mode
        if args['fq2'] is None:
            fq_prefix, fq_clean, fq_log, fq_untrim = self.trim_init()

            ## untrim
            if args['rm_untrim']:
                untrim_arg = '--untrimmed-output=%s --cores=1 -o %s' % (fq_untrim, fq_clean)
            else:
                untrim_arg = '--cores=%s -o %s' % (args['threads'], fq_clean)

            ## cut before
            if args['cut_before_trim'] > '0':
                cut_arg ='--cut %s'  % args['cut_before_trim']
            else:
                cut_arg = ''

            arg_cmd = ' '.join([cut_basic, ad3_arg, untrim_arg, cut_arg, args['fq1']])

        ## PE mode
        elif os.path.exists(args['fq2']):
            fq_prefix, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim = self.trim_init()

            AD3_list = [args['AD3'], ]
            if args['adapter_sliding']:
                AD3_list = self.adapter_chopper(args['AD3'])
            AD3_arg = ' '.join(['-A %s' % i for i in AD3_list])

            ## untrim
            if args['rm_untrim']:
                untrim_arg = '--untrimmed-output=%s --untrimmed-paired-output=%s \
                    --cores=1 -o %s --paired-output=%s' % (fq1_untrim, fq2_untrim,
                        fq1_clean, fq2_clean)
            else:
                untrim_arg = '--cores=%s -o %s --paired-output=%s' % (args['threads'],
                        fq1_clean, fq2_clean)

            ## cut before
            if args['cut_before_trim'] > '0':
                cut_arg ='--cut %s'  % args['cut_before_trim']
            else:
                cut_arg = ''

            ## merge
            arg_cmd = ' '.join([cut_basic, ad3_arg, AD3_arg, untrim_arg, cut_arg, args['fq1'], args['fq2']])
        else:
            raise ValueError('--fq2 value illegal')

        return arg_cmd


    def trim_se(self):
        """Trimming SE reads using cutadapt"""
        args = self.kwargs.copy()
        fq_prefix, fq_clean, fq_log, fq_untrim = self.trim_init()
        logging.info('trimming SE reads: %s' % fq_prefix)

        arg_cmd = self.get_cutadapt_cmd()
        if os.path.exists(fq_clean) and args['overwrite'] is False:
            logging.info('file exists, cutadapt skipped: %s' % fq_prefix)
        else:
            with open(fq_log, 'wt') as fo:
                p1 = subprocess.run(shlex.split(arg_cmd), stdout=fo, stderr=fo)
            Cutadapt_log(fq_log).saveas()
        return fq_clean


    def trim_pe(self):
        """Trimming Paired-end reads using cutadapt
        example: 
        cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq

        -A  see -a for read1
        -G  see -g for read1
        -B  see -b for read1
        -u  see --cut for read1

        -p  write the second read to FILE
        --untrimmed-paired-output  write the second read in a pair to this FILE when no
            adapter was found in the first read. together with --untrimmed-output
        --too-short-paired-output  write the second read in a pair to this FILE if pair
            is too short, together with --too-short-output
        --too-long-paired-output    write the second read in a pair to this FILE if pair
            is too long, together with --too-long-output
        """
        args = self.kwargs.copy()
        fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim = self.trim_init()
        logging.info('trimming PE reads: %s' % fq1_name)

        arg_cmd = self.get_cutadapt_cmd()
        if os.path.exists(fq1_clean) and os.path.exists(fq2_clean) and args['overwrite'] is False:
            logging.info('file exists, cutadapt skipped: %s' % fq1_name)
        else:
            with open(fq_log, 'wt') as fo:
                p1 = subprocess.run(shlex.split(arg_cmd), stdout=fo, stderr=fo)
            # Cutadapt_log(fq_log).saveas()
        return [fq1_clean, fq2_clean]


    def rm_duplicate(self, fq_in=None, fq_out=None):
        """Remove possible PCR duplicates, using fastx_collapse
        collapse reads before trimming
        ****
        only support SE reads (not PE reads)
        """
        logging.info('removing PCR duplicates')
        args = self.kwargs.copy()
        path_rmdup = os.path.join(args['path_out'], 'rm_PCR_dup')
        assert is_path(path_rmdup)
        pkg_dir = os.path.split(goldclip.__file__)[0]
        fa2fq = os.path.join(pkg_dir, 'bin', 'fasta_to_fastq.pl')

        if fq_in is None:
            fq_in = args['fq1']
        if fq_out is None:
            fq_out = os.path.join(path_rmdup, os.path.basename(fq_in))
        fq_prefix, fq_clean, fq_log, fq_untrim = self.trim_init(fq_in)

        fastx_collapser = which('fastx_collapser')

        if os.path.exists(fq_out) and args['overwrite'] is False:
            logging.info('file exists, skip rm_dup: %s' % fq_prefix)
        else:
            # create new clean fastq file
            with open(fq_out, 'wt') as ff:
                c1 = '%s -Q33 -i %s' % (fastx_collapser, fq_in)
                c2 = 'perl %s -' % fa2fq
                p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE)
                p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, stdout=ff)
                px = p2.communicate()
        return fq_out


    def cut_ends(self, fq_in=None, fq_out=None, rm_too_short=True):
        """Cut N-bases at either tail of fastq after cutadapt
        Using python to process fastq files
        support SE reads only
        """
        logging.info('Trimming ends of fastq file')
        args = self.kwargs.copy()
        path_cut_ends = os.path.join(args['path_out'], 'cut_ends')
        assert is_path(path_cut_ends)

        if fq_in is None:
            raise ValueError('either fq_in or fq_out are required')
        if fq_out is None:
            fq_out = os.path.join(path_cut_ends, os.path.basename(fq_in))

        args_cut_ends = self.cut_parser(args['cut_after_trim'], False)

        if len(args_cut_ends) == 2:
            trim_5, trim_3 = args_cut_ends
        elif len(args_cut_ends) == 1:
            trim_5 = 0 if(int(args_cut_ends[0]) < 0) else args_cut_ends[0]
            trim_3 = 0 if(int(args_cut_ends[0]) > 0) else args_cut_ends[0]
        else:
            trim_5 = trim_3 = 0

        with open(fq_in, 'rt') as fi, open(fq_out, 'wt') as ff:
            while True:
                try:
                    fq_id, fq_seq, fq_plus, fq_qual = [next(fi).strip(),
                                                       next(fi).strip(),
                                                       next(fi).strip(),
                                                       next(fi).strip()]
                    if len(fq_seq) < args['len_min'] + abs(trim_5) + abs(trim_3):
                        if rm_too_short:
                            continue # skip short reads
                        else:
                            fq_seq = 'A'
                            fq_qual = 'J' # quality
                    else:
                        fq_seq = fq_seq[trim_5:trim_3] if(trim_3 < 0) else fq_seq[trim_5:]
                        fq_qual = fq_qual[trim_5:trim_3] if(trim_3 < 0) else fq_qual[trim_5:]
                    # trim to length
                    if args['trim_to_length'] > args['len_min']:
                        n_right = min([args['trim_to_length'], len(fq_seq)])
                        fq_seq = fq_seq[:n_right]
                        fq_qual = fq_qual[:n_right]
                    ff.write('\n'.join([fq_id, fq_seq, fq_plus, fq_qual]) + '\n')
                except StopIteration:
                    break

        return fq_out


    def cut_ends2(self, fq_in=None, fq_out=None):
        """Cut reads from 3' end to specific length, 
        save reads with maximum length
        """
        logging.info('Trimming reads to specific length')
        args = self.kwargs.copy()
        path_cut_ends = os.path.join(args['path_out'], 'cut_ends2')
        assert is_path(path_cut_ends)

        if fq_in is None:
            fq_in = args['fq1']
        if fq_out is None:
            fq_out = os.path.join(path_cut_ends, os.path.basename(fq_in))
        fq_prefix, fq_clean, fq_log, fq_untrim = self.trim_init(fq_in)

        with open(fq_in, 'rt') as fi, open(fq_out, 'wt') as ff:
            while True:
                try:
                    fq_id, fq_seq, fq_plus, fq_qual = [next(fi).strip(),
                                                       next(fi).strip(),
                                                       next(fi).strip(),
                                                       next(fi).strip()]
                    n_right = min([args['trim_to_length'], len(fq_seq)])
                    fq_seq = fq_seq[:n_right]
                    fq_qual = fq_qual[:n_right]
                    ff.write('\n'.join([fq_id, fq_seq, fq_plus, fq_qual]) + '\n')
                except StopIteration:
                    break
        return fq_out


    def run_se(self):
        """Run cutadapt for SE reads
        1. cut N-bases before trimming (optional)
        2. trim 3' adapter
        3. remove PCR duplicates (optional)
        4. cut N-bases after trimming (optional)
        5. cut reads into specific length (optional)
        """
        logging.info('Trimming reads: SE mode')
        args = self.kwargs.copy()

        fq_prefix, fq_clean, fq_log, fq_untrim = self.trim_init()
        fq_in = args['fq1']

        ## 1. cut N-bases before trim
        ## 2. trim 3' adapter
        fq_out = self.trim_se() # fq_out == fq_clean

        ## 3. remove PCR dup
        if args['rm_dup']:
            fq_tmp1 = fq_in + '.before_rmdup.tmp'
            os.rename(fq_clean, fq_tmp1)
            fq_clean = self.rm_duplicate(fq_in=fq_tmp1, fq_out=fq_clean)
            os.remove(fq_tmp1)

        ## 4. cut N-bases after trim
        if not args['cut_after_trim'] == '0':
            fq_tmp2 = fq_in + '.before_cut_after_trim.tmp'
            os.rename(fq_clean, fq_tmp2)
            fq_clean = self.cut_ends(fq_in=fq_tmp2, fq_out=fq_clean)
            os.remove(fq_tmp2)

        ## 5. cut to specific length
        if args['trim_to_length'] > args['len_min']:
            fq_tmp3 = fq_in + '.before_trim_to_length.tmp'
            os.rename(fq_clean, fq_tmp3)
            fq_clean = self.cut_ends2(fq_in=fq_tmp3, fq_out=fq_clean)
            os.remove(fq_tmp3)

        return fq_clean


    def run_pe(self):
        """Run cutadapt for PE reads
        1. cut N-bases before trimming (optional)
        2. trim 3' adapter
        2. cut N-bases after trimming (optional)
        """
        logging.info('Trimming reads: PE mode')
        args = self.kwargs.copy()
        fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim = self.trim_init()
        fq_in = args['fq1']

        ## 1. cut N-bases before trim
        ## 2. trim 3' adapter
        fq1_clean, fq2_clean = self.trim_pe() # fq_out == fq_clean

        ## 3. cut N-bases after trim
        if not args['cut_after_trim'] == '0':
            fq1_tmp1 = fq1_clean + '.before_cut_after_trim.tmp'
            fq2_tmp2 = fq2_clean + '.before_cut_after_trim.tmp'
            os.rename(fq1_clean, fq1_tmp1)
            os.rename(fq2_clean, fq2_tmp1)
            fq1_clean = self.cut_ends(fq_in=fq1_tmp1, fq_out=fq1_clean)
            fq2_clean = self.cut_ends(fq_in=fq1_tmp1, fq_out=fq2_clean)
            os.remove(fq1_tmp1)
            os.remove(fq2_tmp1)

        ## 4. cut to specific length
        if args['trim_to_length'] > args['len_min']:
            fq1_tmp2 = fq1_clean + '.before_trim_to_length.tmp'
            fq2_tmp2 = fq2_clean + '.before_trim_to_legnth,tmp'
            os.rename(fq1_clean, fq1_tmp2)
            os.rename(fq2_clean, fq2_tmp2)
            fq1_clean = self.cut_ends2(fq_in=fq1_tmp2, fq_out=fq1_clean)
            fq2_clean = self.cut_ends2(fq_in=fq2_tmp2, fq_out=fq2_clean)

        return [fq1_clean, fq2_clean]


    def run(self):
        """Run trimming for all"""
        args = self.kwargs.copy()
        fq_prefix = file_prefix(args['fq1'])[0]

        ## save arguments
        args_file = os.path.join(args['path_out'], fq_prefix + '.arguments.txt')
        args_pickle = os.path.join(args['path_out'], fq_prefix + '.arguments.pickle')
        if args_checker(args, args_pickle) and args['overwrite'] is False:
            logging.info('files exists, arguments not changed, trimming skipped - %s' % fq_prefix)
            return True #!!! fastq files
        else:
            args_logger(args, args_file, True) # update arguments.txt

        ## SE mode
        if args['fq2'] is None:
            fq_return = self.run_se()
            fq1_clean = fq_return

        ## PE mode
        elif os.path.exists(args['fq2']):
            fq_return = self.run_pe()
            fq1_clean = fq_return[0]

        else:
            raise Exception('Illegal --fq2 argument: %s' % args['fq2'])

        ## save read count
        fq_prefix = file_prefix(args['fq1'])[0]
        fq_count_file = os.path.join(args['path_out'], fq_prefix + '.clean_reads.txt')
        if not os.path.exists(fq_count_file) or args['overwrite'] is True:
            fq_count = int(file_row_counter(fq1_clean) / 4)
            with open(fq_count_file, 'wt') as fo:
                fo.write(str(fq_count) + '\n')

        return fq_return



## EOF
