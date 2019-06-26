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


2019-04-22

add modules

Trimmer, 
Cutter,
DupRemover,
LogParser,


1. Trim adapters (ad-3, ad-5, ...)
2. cut N-bases after trim
3. rm duplicates
4. cut N-bases after rm-dup

"""

import os
import sys
import re
import shlex
import logging
import subprocess
from utils_parser import *
from helper import *
from arguments import args_init


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)


class Fastx(object):
    """Collection of tools to manipulate fastx file
    1. trimmer: cutadapt [trimmomatic, ...]
    2. collapse: fastx_collapse [fastx_toolkit]
    3. collapse: seqkit rmdup [faster, discard read numbers]
    4. fq2fa: fastq_to_fasta [fastx_toolkit]
    5. fa2fq: [fasta_to_fastq.pl]
    6. revcomp: [seqkit seq]
    7. sample: [seqkit sample -n]
    ...
    """
    def __init__(self, fx_input, path_out, **kwargs):

        ## init args
        overwrite = kwargs.get('overwrite', False)
        output_fasta = kwargs.get('output_fasta', False)
        gzip = kwargs.get('gzip', True)

        ## sequence type
        fx_type = seq_type(fx_input) # fastq, fasta, None
        if fx_type is None:
            raise Exception('unknown input file type. expect [fastq, fasta]')

        ## output filename
        fx_prefix = file_prefix(fx_input)[0]
        if output_fasta:
            output_suffix = 'fasta'
        else:
            output_suffix = fx_type
        fx_output = os.path.join(path_out, fx_prefix + '.' + output_suffix)

        if gzip:
            fx_output += '.gz'

        ## check file exist
        is_path(path_out)
        output_exists = False # default
        if os.path.exists(fx_output) and overwrite is False:
            log.info('output file exists. try overwrite=True ?')
            output_exists = True
            sys.exit('update arguments, re-run again...')

        ## open file
        if is_gz(fx_input):
            cmd_open = 'zcat ' + fx_input
        else:
            cmd_open = 'cat ' + fx_input

        self.input = fx_input
        self.path_out = path_out
        self.output = fx_output
        self.fx_type = fx_type
        self.gzip = gzip
        self.overwrite = overwrite
        self.output_fasta = output_fasta
        self.output_exists = output_exists
        self.cmd_open = cmd_open
        self.args = kwargs


    def cut(self, cut_left=0, cut_right=0, cut_to_length=0, 
        len_min=20, cut_from_right=True, discard_tooshort=True):
        """Cut N-bases from fastq or fasta file
        cut_left: 5' end of read
        cut_right: 3' end of read
        cut_to_length: cut from right (default)

        1. cut first N bases
        2. cut last N bases
        3. cut from right, keep specific length
        4. specific length (seqkit subseq, bioawk, ...)

        Using fastx_toolkit
        http://hannonlab.cshl.edu/fastx_toolkit/
        """
        ## external program
        cut_exe = which('fastx_trimmer')
        seqkit_exe = which('seqkit')
        if cut_exe is None:
            raise Exception('command not found in $PATH: fastx_trimmer')
        if seqkit_exe is None:
            raise Exception('command not found in $PATH: seqkit')

        log.info('trimming ends: cut5=%d, cut3=%d, cut_to_length=%d' % (cut_left, cut_right, cut_to_length))

        ## fastq file
        if self.fx_type == 'fastq':
            cut_exe += ' -Q33'

        if cut_to_length > 0:
            if cut_to_length < len_min:
                raise Exception('arguments conflic, cut_to_length=%s short than len_min=%s' % (cut_to_length, len_min))
            # cut right
            if cut_from_right:
                cmd = '{} -l {}'.format(
                    cut_exe,
                    cut_to_length)
                if discard_tooshort:
                    cmd += ' | {} -m {}'.format(
                        seqkit_exe,
                        len_min)
            else:
                # cut left
                # rev-comp, then cut-right, rev-comp
                revcomp = which('fastx_reverse_complement')
                if revcomp is None:
                    raise Exception('command not found in $PATH: fastx_reverse_complement')
                cmd = '{} | ' # revcomp
                cmd += '{} -l {} | ' # cut 
                cmd += '{} '
                cmd = cmd.foramt(
                    revcomp,
                    cut_exe,
                    cut_to_length,
                    revcomp)
                if discard_tooshort:
                    cmd += ' | {} -m {}'.format(
                        seqkit_exe,
                        len_min)

        elif cut_left > 0 and cut_right > 0:
            cmd = '{} -f {} | '
            cmd += '{} -t {} '
            cmd = cmd.format(
                cut_exe,
                cut_left + 1,
                cut_exe,
                cut_right)
            if discard_tooshort:
                cmd += '-m {}'.format(len_min)

        elif cut_left > 0:
            cmd = '{} -f {} '.format(
                cut_exe,
                cut_left + 1)
            if discard_tooshort:
                cmd += ' | {} -m {}'.format(
                    seqkit_exe,
                    len_min)

        elif cut_right > 0:
            cmd = '{} -t {} '
            if discard_tooshort:
                cmd += ' -m {} '.format(len_min)
            cmd = cmd.format(
                cut_exe,
                cut_right)

        else:
            pass

        # run program
        cmd = '{} | {} > {}'.format(
            self.cmd_open,
            cmd,
            self.output)

        run_shell_cmd(cmd)

        return self.output


    def rmdup(self):
        """Remove PCR duplicates (sequence) 
        fastx_collapser
        fasta_to_fastq
        """
        fastx_collapser = which('fastx_collapser')
        if fastx_collapser is None:
            raise Exception('command not found in $PATH: fastx_collapser')

        log.info('removing PCR dup')

        # cmd = fastx_collapser
        if self.fx_type == 'fastq':
            fastx_collapser += ' -Q33'

        # require
        # 1. collapser
        cmd = '{} | {} > {}'.format(
            self.cmd_open,
            fastx_collapser,
            self.output)
        run_shell_cmd(cmd)

        if not self.output_fasta:
            # fa2fq
            out_tmp = self.output + '.tmp'
            os.rename(self.output, out_tmp)
            self.fa2fq(out_tmp, self.output)

        return self.output


    def revcomp(self):
        """Reverse complement sequences"""
        revcomp = which('fastx_reverse_complement')
        if revcomp is None:
            raise Exception('command not found in $PATH: fastx_reverse_complement')

        if self.fx_type == 'fastq':
            revcomp += ' -Q33'

        cmd = '{} | {} > {}'.format(
            self.cmd_open, 
            revcomp,
            self.output)

        run_shell_cmd(cmd)

        return self.output


    def sample(self, n=1000, p=0.01):
        """Extract sample size of fastx file
        using seqkit tool

        $ seqkit sample -n 1000

        $ seqkit sample -p 0.01

        """
        assert isinstance(n, int)
        assert isinstance(p, float)

        seqkit_exe = which('seqkit')
        if seqkit_exe is None:
            raise Exception('command not found in $PATH: seqkit')

        log.info('extracting sub-sample: n=%d' % n)

        ## warning
        if n > 1000000 or p > 0.1:
            log.warning('choose a smaller number.')

        cmd = '{} | {} sample -n {} > {}'.format(
            self.cmd_open,
            seqkit_exe,
            n,
            self.output)
        run_shell_cmd(cmd)

        return self.output


    def fa2fq(self, fa, fq):
        """Convert fasta to fastq, quality=J
        fasta

        >name
        sequence


        fastq

        @name
        sequence
        +
        quality
        """
        with xopen(fa, 'rb') as fi, xopen(fq, 'wb') as fo:
            name = None
            for line in fi:
                if line.startswith('>'):
                    if isinstance(name, str):
                        fo.write('\n'.join(['@' + name, sequence, '+', seq_quality]) + '\n')
                    name = line.strip()
                    sequence = ""
                    seq_quality = ""

                else:
                    sequence += line.strip()
                    seq_quality = 'J' * len(sequence)
            # the last record
            fo.write('\n'.join(['@' + name, sequence, '+', seq_quality]) + '\n')

        return fq


    def fq2fa(self):
        seqkit_exe = which('seqkit')
        if seqkit_exe is None:
            raise Exception('command not found in $PATH: seqkit')

        cmd = '{} | {} fq2fa > {}'.format(
            self.cmd_open,
            seqkit_exe,
            self.output)
        run_shell_cmd(cmd)

        return self.output


class Cutadapt(object):
    """Cut adapters using cutadapt
    1. cut 3-ad
    2. cut 5-ad
    3. save log

    both SE and PE reads

    """
    def __init__(self, fq1, adapter3, path_out=None, len_min=15, fq2=None, **kwargs):
        """Parsing the parameters for reads trimming
        support both SE and PE reads
        """
        assert isinstance(fq1, str) # process one file
        args = args_init(kwargs) # default parameters
        # print(args)
        args['fq1'] = fq1
        args['fq2'] = fq2
        args['adapter3'] = adapter3
        args['path_out'] = path_out
        args['len_min'] = len_min

        ## validate options
        assert is_path(path_out)
        if args['cut_to_length'] > 0 and args['cut_to_length'] < len_min:
            raise Exception('illegal length, --trim-to-length %s --len-min %s' % (args['cut_to_length'], len_min))

        self.kwargs = args # global


    def cutadapt_init(self, fq1_in=None, fq2_in=None):
        """Prepare directory for trimming, filename for each file
        fq
        fq_clean
        fq_log
        fq_untrim
        fq_too_short
        fq_too_long
        """
        args = self.kwargs.copy()
        if fq1_in is None:
            fq1_in = args['fq1']
        if fq2_in is None:
            fq2_in = args['fq2']
        
        ## fq1
        fq1_prefix = file_prefix(fq1_in)[0]
        fq1_type = seq_type(fq1_in)
        ## file type
        if fq1_type is None:
            raise Exception('unknown filetype, only support: [fasta, fastq]')
        ## file name
        fq1_clean_prefix = os.path.join(args['path_out'], fq1_prefix)
        if args['keep_name']:
            fq1_clean = fq1_clean_prefix + '.' + fq1_type
        else:
            fq1_clean = fq1_clean_prefix + '.clean.' + fq1_type
        ## other files
        fq1_log = fq1_clean_prefix + '.cutadapt.log'
        fq1_untrim = fq1_clean_prefix + '.untrim.' + fq1_type
        fq1_too_short = fq1_clean_prefix + '.tooshort.' + fq1_type
        fq1_too_long = fq1_clean_prefix + '.toolong.' + fq1_type
        
        # SE mode
        if args['fq2'] is None:
            return [fq1_prefix, fq1_clean, fq1_log, 
                fq1_untrim, fq1_too_short, fq1_too_long]

        # PE mode
        elif os.path.exists(args['fq2']):
            fq2_prefix = file_prefix(fq2_in)[0]
            fq2_type = seq_type(fq2_in)
            ## file type
            if fq2_type is None:
                raise Exception('unknown filetype, only support: [fasta, fastq]')
            ## file name
            fq2_clean_prefix = os.path.join(args['path_out'], fq2_prefix)
            if args['keep_name']:
                fq2_clean = fq2_clean_prefix + '.' + fq2_type
            else:
                fq2_clean = fq2_clean_prefix + '.clean.' + fq2_type
            ## other files
            fq2_untrim = fq2_clean_prefix + '.untrim.' + fq2_type
            fq2_too_short = fq2_clean_prefix + '.tooshort.' + fq2_type
            fq2_too_long = fq2_clean_prefix + '.toolong.' + fq2_type
            return [fq1_prefix, fq1_clean, fq2_clean, fq1_log, 
                fq1_untrim, fq2_untrim, fq1_too_short, fq2_too_short, 
                fq1_too_long, fq2_too_long]
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


    def get_cutadapt_cmd(self):
        """Parse arguments, create cutadapt command line
        SE: 
        cutadapt -a ADAPTER_FWD -A ADAPTER_REV reads.fastq 1> clean.fq 2> log

        --minimum-length=LENGTH
        --max-n=COUNT
        --discard-untrimmed
        --discard-trimmed
        --length=LENGTH
        --too-short-output=FILE
        --too-long-output=FILE
        --untrimmed-output=FILE

        PE:
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

        ## command
        cutadapt = which('cutadapt')

        ## basic
        cmd_basic = '%s -m %s -q %s --overlap=%s --error-rate=%s \
            --trim-n --max-n=0.1 --times=%s' % (
            cutadapt,
            args['len_min'],
            args['qual_min'],
            args['overlap'],
            args['error_rate'],
            args['trim_times'])

        ## option-1: adapter
        ## 3-adapters
        ad3_list = args['adapter3']

        ## short version of adapter
        if args['adapter_sliding']:
            i_list = []
            for i in ad3_list:
                i_list.extend(self.adapter_chopper(i))
            # ad3_list = self.adapter_chopper(args['adapter3'])
            # ad3_list.extend(ad3_ext)
            ad3_list = i_list
        args_ad3 = ' '.join(['-a ' + i for i in ad3_list])

        ## default: SE mode
        if args['fq2'] is None:
            fq_prefix, fq_clean, fq_log, fq_untrim, fq_too_short, fq_too_long = self.cutadapt_init()

            ## option-2: filter
            ## 2.1 untrim
            if args['save_untrim']:
                args_untrim = '--untrimmed-output=%s --cores=1' % fq_untrim
            elif args['rm_untrim']:
                args_untrim = '--discard-untrimmed --cores=%s' % args['threads']
            else:
                args_untrim = '--cores=%s' % args['threads']

            ## 2.2 too-short
            if args['save_too_short']:
                args_too_short = '--too-short-output=%s' % fq_too_short
            else:
                args_too_short = ''

            ## 2.3 too-long
            if args['save_too_long']:
                args_too_long = '--too-long-output=%s' % fq_too_long
            else:
                args_too_long = ''

            ## 2.4 cut-to-length
            if args['cut_to_length'] >= args['len_min']:
                args_cut_to_length = '--length=%s' % args['cut_to_length']
            else:
                args_cut_to_length = ''

            ## SE command
            cmd = ' '.join([cmd_basic, args_ad3, args_untrim, args_too_short, 
                args_too_long, args_cut_to_length, args['fq1']])

        ## PE command
        elif args['fq2']:
            tmp = self.cutadapt_init()
            ## too-long in one-line
            fq1_prefix, fq1_clean, fq2_clean, fq1_log = tmp[:4]
            fq1_untrim, fq2_untrim, fq1_too_short, fq2_too_short = tmp[4:8]
            fq1_too_long, fq2_too_long = tmp[8:]

            ## option-1: adapter
            AD3_list = args['AD3']
            ## short version of adapter
            if args['adapter_sliding']:
                i_list = []
                for i in AD3_list:
                    i_list.extend(self.adapter_chopper(i))
                AD3_list = i_list
            args_AD3 = ' '.join(['-A ' + i for i in AD3_list])
            # args_ad3 = ' '.join(['-a ' + i for i in ad3_list])

            ## option-2: filter
            ## 2.1 untrim 
            args_untrim = '-o %s --paired-output=%s' % (fq1_clean, fq2_clean)
            if args['save_untrim']:
                args_untrim += ' --untrimmed-output=%s \
                    --untrimmed-paired-output=%s --cores=1' % (fq1_untrim, fq2_untrim)
            elif args['rm_untrim']:
                args_untrim += ' --discard-untrimmed --cores=%s' % args['threads']
            else:
                args_untrim += ' --cores=%s' % args['threads']

            ## 2.2 too-short
            if args['save_too_short']:
                args_too_short = '--too-short-output=%s \
                    --too-short-paired-output=%s' % (
                        fq1_too_short,
                        fq2_too_short)
            else:
                args_too_short = ''

            ## 2.3 too-long
            if args['save_too_long']:
                args_too_long = '--too-long-output=%s \
                    --too-long-paired-output=%s' % (
                        fq1_too_long,
                        fq2_too_long)
            else:
                args_too_long = ''

            ## 2.4 cut-to-length
            if args['cut_to_length'] >= args['len_min']:
                args_cut_to_length = '--length=%s' % args['cut_to_length']
            else:
                args_cut_to_length = ''                    

            cmd = ' '.join([cmd_basic, args_ad3, args_AD3, args_untrim, 
                args_too_short, args_too_long, args_cut_to_length, 
                args['fq1'], args['fq2']])
        # others
        else:
            raise Exception('illegal --fq2')

        return cmd


    def trim_se(self):
        """Trimming SE reads using cutadapt
        cutadapt -a ADAPTER_FWD -A ADAPTER_REV reads.fastq 1> clean.fq 2> log
        """
        args = self.kwargs.copy()
        fq_prefix, fq_clean, fq_log, fq_untrim, fq_too_short, fq_too_long = self.cutadapt_init()
        log.info('trimming SE reads: %s' % fq_prefix)

        if args['gzip']:
            fq_clean += '.gz'

        ## run program
        cmd = self.get_cutadapt_cmd()
        if os.path.exists(fq_clean) and args['overwrite'] is False:
            log.info('file exists, cutadapt skipped: %s' % fq_prefix)
        else:
            if args['gzip']:
                cmd += ' 2> {} | gzip -nc > {}'
            else:
                cmd += ' 2> {} 1> {}'
            cmd = cmd.format(fq_log, fq_clean)
            run_shell_cmd(cmd)
            Cutadapt_log(fq_log).saveas()

        return fq_clean


    def trim_pe(self):
        """Trimming Paired-end reads using cutadapt
        example: 
        cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq
        """
        args = self.kwargs.copy()
        fq1_prefix, fq1_clean, fq2_clean, fq1_log = self.cutadapt_init()[:4]
        log.info('trimming PE reads: %s' % fq1_prefix)

        cmd = self.get_cutadapt_cmd()
        ## check exists
        fq1_check = fq1_clean
        fq2_check = fq2_clean
        if args['gzip']:
            fq1_check += '.gz'
            fq2_check += '.gz'
        ## compress fastq file after cutadapt
        if os.path.exists(fq1_check) and os.path.exists(fq2_check) and args['overwrite'] is False:
            log.info('file exists, cutadapt skipped: %s' % fq1_prefix)
        else:
            cmd += ' 1> {}'.format(fq1_log)
            run_shell_cmd(cmd)
            if args['gzip']:
                fq1_clean = gzip_file2(fq1_clean, rm=True)
                fq2_clean = gzip_file2(fq2_clean, rm=True)
            # Cutadapt_log(fq_log).saveas()
        return [fq1_clean, fq2_clean]


    def run(self):
        """Run trimming for all"""
        args = self.kwargs.copy()

        args_init = self.cutadapt_init()
        fq_prefix, fq1_clean = args_init[:2]

        ## save arguments
        args_file = os.path.join(args['path_out'], fq_prefix + '.arguments.txt')
        args_pickle = os.path.join(args['path_out'], fq_prefix + '.arguments.pickle')
        if args_checker(args, args_pickle) and args['overwrite'] is False:
            log.info('arguments not changed, trimming skipped - %s' % fq_prefix)
        
        if args['fq2'] is None:
            return_fq = self.trim_se() # fq1
            fq1_clean = return_fq
        else:
            return_fq = self.trim_pe() # fq1 fq2
            fq1_clean = return_fq[0]

        ## save read count
        fq_count_file = os.path.join(args['path_out'], fq_prefix + '.clean_reads.txt')
        if not os.path.exists(fq_count_file) or args['overwrite'] is True:
            fq_count = int(file_row_counter(fq1_clean) / 4)
            with open(fq_count_file, 'wt') as fo:
                fo.write(str(fq_count) + '\n')

        return return_fq


class Trimmer(object):
    """Processing fastq files
    
    General
      a. low quality bases at 3' 
      b. trim-n

    1. CLIP reads

      a. trim adapter (sliding)
      b. cut inline-barcode
      c. cut N-bases (5', 3' ends)
      e. collapse (remove duplicates)
      e. cut random barcode

    2. RNAseq reads

      a. trim adapter
      b. cut N-bases

    3. ChIPseq reads
      a. trim adapter
      b. cut N-bases (5', 3' ends)

    4. smRNAseq reads
      a. trim adapter
      b. discard untrimmed reads
      c. cut N-bases

    5. ATACseq reads
      a. trim adapter (sliding)

    """

    def __init__(self, fq1, adapter3, path_out=None, len_min=15, **kwargs):
        """Parsing the parameters for reads trimming
        support both SE and PE reads
        """
        assert isinstance(fq1, str) # process one file
        args = args_init(kwargs) # default parameters
        # print(args)
        args['fq1'] = fq1
        args['adapter3'] = adapter3
        args['path_out'] = path_out
        args['len_min'] = len_min

        ## validate options
        assert is_path(path_out)
        if args['cut_to_length'] > 0 and args['cut_to_length'] < len_min:
            raise Exception('illegal length, --cut-to-length %s --len-min %s' % (args['cut_to_length'], len_min))

        self.kwargs = args # global


    def cut_parser(self, s, return_arguments=False):
        """Cut N-bases, before adapter trimming
        parse the argument
        if cut_para is True,
        return the cut numbers in --cut {arg} version

        else
        return the numbers        
        """
        ## match
        p = re.compile('^\-?\d+$|^\d+\,\-\d+$') # 5,  -4,   5,-4
        if not p.match(s):
            raise Exception('illegal argument: %s, supported [7 | 7,-7 | -7]' % s)

        ## length 1, 2
        ## return int
        n = s.strip().split(',')
        if len(n) == 2:
            cut5, cut3 = list(map(int, n))
            cut3 = abs(cut3)
        elif '-' in n[0]: # len() == 1, minus
            cut5 = 0
            cut3 = abs(int(n[0]))
        else:
            cut5 = int(n[0])
            cut3 = 0
        ##
        return_code = [cut5, cut3]

        if return_arguments:
            return_code = ['--cut %d' % i for i in return_code]

        return return_code


    def trim_se(self):
        """Define the processing steps
        a. trim adapter
        b. cut N bases
        c. collapse
        d. cut N bases
        e. cut to length
        """
        args = self.kwargs.copy()

        ## output file 
        fq1_prefix = file_prefix(args['fq1'])[0]
        fx_type = seq_type(args['fq1'])
        return_output_name = fq1_prefix + '.' + fx_type
        if args['gzip']:
            return_output_name += '.gz'
        return_output = os.path.join(args['path_out'], return_output_name)

        ## save arguments
        args_file = os.path.join(args['path_out'], fq1_prefix + '.arguments.txt')
        args_pickle = os.path.join(args['path_out'], fq1_prefix + '.arguments.pickle')
        if args_checker(args, args_pickle) and os.path.exists(return_output) and args['overwrite'] is False:
            log.info('file exists, arguments not changed, trimming skipped : %s' % fq1_prefix)
            return return_output
        else:
            args_logger(args, args_file, True) # update arguments.txt

        ## determine the len_min 
        ## cut_after_trim
        ## cut_after_rmdup
        cut1, cut2 = self.cut_parser(args['cut_after_trim'])
        cut3, cut4 = self.cut_parser(args['cut_after_rmdup'])

        ## update arguments for Cutadapt
        path_out = args.pop('path_out')
        fq1 = args.pop('fq1')
        adapter3 = args.pop('adapter3')
        len_min = args.pop('len_min')
        len_min_cutadapt = sum([len_min, cut1, cut2, cut3, cut4])
        ## ignore cut-to-length
        cut_to_length = args.pop('cut_to_length') # save for further analysis
        args['cut_to_length'] = 0 # update

        ## 0. input file
        count_input = fx_counter(fq1)

        ## 1. trim adapter
        path_cutadapt = os.path.join(path_out, '1_trim_adapter')
        if args['not_trim_adapter']:
            return_trim = fq1
            count_trim = count_input
        else:
            return_trim = Cutadapt(
                fq1, 
                adapter3, 
                path_cutadapt, 
                len_min_cutadapt, **args).run()
            count_trim = fx_counter(return_trim)

        ## 2. cut N bases
        path_cut1 = os.path.join(path_out, '2_cut')
        return_cut1 = os.path.join(path_cut1, os.path.basename(return_trim))
        if args['cut_after_trim'] == '0' or is_empty_file(return_trim):
            return_cut1 = return_trim
            count_cut1 = count_trim
        else:
            cut5, cut3 = self.cut_parser(args['cut_after_trim'])
            Fastx(return_trim, path_cut1, **args).cut(cut_left=cut5, 
                cut_right=cut3, len_min=len_min)
            count_cut1 = fx_counter(return_cut1)

        ## 3. collapse
        path_rmdup = os.path.join(path_out, '3_rmdup')
        return_rmdup = os.path.join(path_rmdup, os.path.basename(return_trim))
        if not args['rmdup'] or is_empty_file(return_cut1):
            return_rmdup = return_cut1
            count_rmdup = count_cut1
        else:
            Fastx(return_cut1, path_rmdup, **args).rmdup()
            count_rmdup = fx_counter(return_rmdup)
        
        ## 4. cut Nbases
        path_cut2 = os.path.join(path_out, '4_cut')
        return_cut2 = os.path.join(path_cut2, os.path.basename(return_trim))
        if args['cut_after_rmdup'] == '0' or is_empty_file(return_rmdup):
            return_cut2 = return_rmdup
            count_cut2 = count_rmdup
        else:
            cut5, cut3 = self.cut_parser(args['cut_after_rmdup'])
            Fastx(return_rmdup, path_cut2, **args).cut(cut_left=cut5, 
                cut_right=cut3, len_min=len_min)
            count_cut2 = fx_counter(return_cut2)

        ## 5. cut to length
        path_cut3 = os.path.join(path_out, '5_cut_to_length')
        return_cut3 = os.path.join(path_cut3, os.path.basename(return_trim))
        if cut_to_length < len_min or is_empty_file(return_cut2):
            return_cut3 = return_cut2
            count_cut3 = count_cut2
        else:
            Fastx(return_cut2, path_cut3, **args).cut(cut_to_length=cut_to_length,
                len_min=len_min)
            count_cut3 = fx_counter(return_cut3)
        
        ## 6. organize output
        ## save results
        # return_output = os.path.join(path_out, os.path.basename(return_trim))
        os.rename(return_cut3, return_output)

        ## remove temp files
        if not args['keep_temp_files']:
            log.info('[remving temp files]')
            rm_file([return_trim, return_cut1, return_rmdup, return_cut2, return_cut3])

        ## 6. summarise
        trim_stat_file = os.path.join(path_out, fq1_prefix + '.trim.stat')
        n_trim = count_input - count_trim
        n_cut1 = count_trim - count_cut1
        n_dup = count_cut1 - count_rmdup
        n_cut2 = count_rmdup - count_cut2
        n_cut3 = count_cut2 - count_cut3

        with open(trim_stat_file, 'wt') as fo:
            fo.write('%13s : %10d : %6.2f%%\n' % ('input', count_input, 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('trim', n_trim, n_trim / count_input * 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('cut_end1', n_cut1, n_cut1 / count_input * 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('duplicate', n_dup, n_dup / count_input * 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('cut_end2', n_cut2, n_cut2 / count_input * 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('cut_to_length', n_cut3, n_cut3 / count_input * 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('output', count_cut3, count_cut3 / count_input * 100))

        ## return
        return return_output


    def trim_pe(self):
        """Define the processing steps
        a. trim adapter
        b. cut N bases
        c. collapse
        d. cut N bases
        e. cut to length
        """
        args = self.kwargs.copy()

        ## output file 
        fq1_prefix = file_prefix(args['fq1'])[0]
        fq2_prefix = file_prefix(args['fq2'])[0]
        fx_type = seq_type(args['fq1'])
        return_output_name1 = fq1_prefix + '.' + fx_type
        return_output_name2 = fq2_prefix + '.' + fx_type
        if args['gzip']:
            return_output_name1 += '.gz'
            return_output_name2 += '.gz'

        # return_output = [os.path.join(args['path_out'], i) for i in return_output_name]
        return_output = [
            os.path.join(args['path_out'], return_output_name1),
            os.path.join(args['path_out'], return_output_name2)]

        ## save arguments
        args_file = os.path.join(args['path_out'], fq1_prefix + '.arguments.txt')
        args_pickle = os.path.join(args['path_out'], fq1_prefix + '.arguments.pickle')
        if args_checker(args, args_pickle) and os.path.exists(return_output) and args['overwrite'] is False:
            log.info('file exists, arguments not changed, trimming skipped : %s' % fq1_prefix)
            return return_output
        else:
            args_logger(args, args_file, True) # update arguments.txt

        ## determine the len_min 
        ## cut_after_trim
        ## cut_after_rmdup
        cut1, cut2 = self.cut_parser(args['cut_after_trim'])
        cut3, cut4 = self.cut_parser(args['cut_after_rmdup'])

        ## update arguments for Cutadapt
        path_out = args.pop('path_out')
        fq1 = args.pop('fq1')
        fq2 = args.pop('fq2')
        adapter3 = args.pop('adapter3')
        len_min = args.pop('len_min')
        len_min_cutadapt = sum([len_min, cut1, cut2, cut3, cut4]) #!!!!
        ## ignore cut-to-length
        cut_to_length = args.pop('cut_to_length') # save for further analysis
        args['cut_to_length'] = 0 # update

        ## 0. input file
        count_input = fx_counter(fq1)

        ## 1. trim adapter
        path_cutadapt = os.path.join(path_out, '1_trim_adapter')
        if args['not_trim_adapter']:
            return_trim = [fq1, fq2]
            count_trim = count_input
        else:
            return_trim = Cutadapt(
                fq1, 
                adapter3, 
                path_cutadapt, 
                len_min_cutadapt, 
                fq2=fq2, **args).run()
            count_trim = fx_counter(return_trim[0])

        ## 2. cut N bases
        path_cut1 = os.path.join(path_out, '2_cut')
        return_cut1 = [os.path.join(path_cut1, os.path.basename(i)) for i in return_trim]
        if args['cut_after_trim'] == '0' or is_empty_file(return_trim[0]):
            return_cut1 = return_trim
            count_cut1 = count_trim
        else:
            cut5, cut3 = self.cut_parser(args['cut_after_trim'])
            [Fastx(i, path_cut1, **args).cut(cut_left=cut5, cut_right=cut3, 
                len_min=len_min, discard_tooshort=False) for i in return_trim]
            count_cut1 = fx_counter(return_cut1[0])

        ## 3. collapse
        path_rmdup = os.path.join(path_out, '3_rmdup')
        return_rmdup = [os.path.join(path_rmdup, os.path.basename(i)) for i in return_cut1]
        if not args['rmdup'] or is_empty_file(return_cut1[0]):
            return_rmdup = return_cut1
            count_rmdup = count_cut1
        else:
            [Fastx(i, path_rmdup, **args).rmdup() for i in return_cut1]
            count_rmdup = fx_counter(return_rmdup[0])

        ## 4. cut Nbases
        path_cut2 = os.path.join(path_out, '4_cut')
        return_cut2 = [os.path.join(path_cut2, os.path.basename(i)) for i in return_rmdup]
        if args['cut_after_rmdup'] == '0' or is_empty_file(return_rmdup[0]):
            return_cut2 = return_rmdup
            count_cut2 = count_rmdup
        else:
            cut5, cut3 = self.cut_parser(args['cut_after_rmdup'])
            [Fastx(i, path_cut2, **args).cut(cut_left=cut5, cut_right=cut3, 
                len_min=len_min, discard_tooshort=False) for i in return_rmdup]
            count_cut2 = fx_counter(return_cut2[0])

        ## 5. cut to length
        path_cut3 = os.path.join(path_out, '5_cut_to_length')
        return_cut3 = [os.path.join(path_cut3, os.path.basename(i)) for i in return_cut2]
        if cut_to_length < len_min or is_empty_file(return_cut2[0]):
            return_cut3 = return_cut2
            count_cut3 = count_cut2
        else:
            [Fastx(i, path_cut3, **args).cut(cut_to_length=cut_to_length,
                len_min=len_min, discard_tooshort=False) for i in return_cut2]
            count_cut3 = fx_counter(return_cut3)
        
        ## 6. organize output
        ## save results
        os.rename(return_cut3[0], return_output[0])
        os.rename(return_cut3[1], return_output[1])

        ## remove temp files
        if not args['keep_temp_files']:
            log.info('[remving temp files]')
            rm_file(return_trim + return_cut1 + return_rmdup + return_cut2 + return_cut3)

        ## 6. summarise
        trim_stat_file = os.path.join(path_out, fq1_prefix + '.trim.stat')
        n_trim = count_input - count_trim
        n_cut1 = count_trim - count_cut1
        n_dup = count_cut1 - count_rmdup
        n_cut2 = count_rmdup - count_cut2
        n_cut3 = count_cut2 - count_cut3

        with open(trim_stat_file, 'wt') as fo:
            fo.write('%13s : %10d : %6.2f%%\n' % ('input', count_input, 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('trim', n_trim, n_trim / count_input * 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('cut_end1', n_cut1, n_cut1 / count_input * 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('duplicate', n_dup, n_dup / count_input * 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('cut_end2', n_cut2, n_cut2 / count_input * 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('cut_to_length', n_cut3, n_cut3 / count_input * 100))
            fo.write('%13s : %10d : %6.2f%%\n' % ('output', count_cut3, count_cut3 / count_input * 100))

        ## return
        return return_output


    def trimmer(self):
        """Choose se or pe"""
        args = self.kwargs.copy()

        ## SE
        if args['fq2'] is None:
            return_output = self.trim_se()
        else:
            return_output = self.trim_pe()

        return return_output

