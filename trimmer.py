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
import subprocess
import logging
from utils_parser import *
from helper import *
from arguments import args_init


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
            logging.info('output file exists. try overwrite=True ?')
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
        len_min=20, cut_from_right=True):
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

        logging.info('trimming ends: cut5=%d, cut3=%d, cut_to_length=%d' % (cut_left, cut_right, cut_to_length))

        ## option-0
        ## filter by length
        cmdx = '%s seq -m %s' % (seqkit_exe, len_min)

        ## option-1
        ## cut left
        ## fastx_trimmer -f <cut_left> -m <len_min>
        cut_left += 1
        cmd1 = '%s -f %s' % (cut_exe, cut_left)

        ## option-2
        ## cut right
        cmd2 = '%s -t %s -m %s' % (cut_exe, cut_right, len_min)

        ## option-3
        ## cut to length
        cmd3 = '%s -l %s' % (cut_exe, cut_to_length)

        ## fastq file
        if self.fx_type == 'fastq':
            cmd1 += ' -Q33'
            cmd2 += ' -Q33'
            cmd3 += ' -Q33'

        ## combine options
        if cut_to_length > 0:
            ## len_min
            if cut_to_length < len_min:
                raise Exception('arguments conflic, cut_to_length=%s short than len_min=%s' % (cut_to_length, len_min))
            ## cut from right
            if cut_from_right:
                with xopen(self.output, 'w') as fo:
                    p0 = subprocess.Popen(self.cmd_open, shell=True, stdout=subprocess.PIPE)
                    p1 = subprocess.Popen(cmd3, shell=True, stdin=p0.stdout, stdout=subprocess.PIPE)
                    p2 = subprocess.Popen(cmdx, shell=True, stdin=p1.stdout, stdout=subprocess.PIPE)
                    for line in p2.stdout:
                        fo.write(line)
                    ret_code = p2.communicate()
            else:
                revcomp = which('fastx_reverse_complement')
                if revcomp is None:
                    raise Exception('command not found in $PATH: fastx_reverse_complement')
                with xopen(self.output, 'w') as fo:
                    p0 = subprocess.Popen(self.cmd_open, shell=True, stdout=subprocess.PIPE)
                    p1 = subprocess.Popen(revcomp, shell=True, stdin=p0.stdout, stdout=subprocess.PIPE)
                    p2 = subprocess.Popen(cmd3, shell=True, stdin=p1.stdout, stdout=subprocess.PIPE)
                    p3 = subprocess.Popen(cmdx, shell=True, stdin=p2.stdout, stdout=subprocess.PIPE)
                    p4 = subprocess.Popen(revcomp, shell=True, stdin=p3.stdout, stdout=subprocess.PIPE)
                    for line in p4.stdout:
                        fo.write(line)
                    ret_code = p4.communicate()

        elif cut_left > 0 and cut_right > 0:
            ## left + right
            with xopen(self.output, 'w') as fo:
                p0 = subprocess.Popen(self.cmd_open, shell=True, stdout=subprocess.PIPE)
                p1 = subprocess.Popen(cmd1, shell=True, stdin=p0.stdout, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(cmd2, shell=True, stdin=p1.stdout, stdout=subprocess.PIPE)
                for line in p2.stdout:
                    fo.write(line)
                ret_code = p2.communicate()

        elif cut_left > 0:
            ## left
            with xopen(self.output, 'w') as fo:
                p0 = subprocess.Popen(self.cmd_open, shell=True, stdout=subprocess.PIPE)
                p1 = subprocess.Popen(cmd1, shell=True, stdin=p0.stdout, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(cmdx, shell=True, stdin=p1.stdout, stdout=subprocess.PIPE)
                # write to file
                for line in p2.stdout:
                    fo.write(line)
                ret_code = p2.communicate() # wait exteral code finish
                # see https://www.jianshu.com/p/11d3a0a9c7d1

        elif cut_right > 0:
            ## right
            with xopen(self.output, 'w') as fo:
                p0 = subprocess.Popen(self.cmd_open, shell=True, stdout=subprocess.PIPE)
                p1 = subprocess.Popen(cmd2, shell=True, stdin=fi, stdout=subprocess.PIPE)
                for line in p1.stdout:
                    fo.write(line)
                p1.communicate()
        else:
            ## no
            pass

        return self.output


    def rmdup(self):
        """Remove PCR duplicates (sequence) 
        fastx_collapser
        fasta_to_fastq
        """
        fastx_collapser = which('fastx_collapser')
        if fastx_collapser is None:
            raise Exception('command not found in $PATH: fastx_collapser')

        logging.info('removing PCR dup')

        cmd = fastx_collapser
        if self.fx_type == 'fastq':
            cmd += ' -Q33'

        if self.output_fasta:
            # require
            # 1. collapser
            with xopen(self.output, 'wt') as fo:
                p0 = subprocess.Popen(self.cmd_open, shell=True, stdout=subprocess.PIPE)
                p1 = subprocess.Popen(cmd, shell=True, stdin=p0.stdout, stdout=subprocess.PIPE)
                for line in p1.stdout:
                    fo.write(line)
                re_code = p1.communicate()
        else:
            # require
            # 1. collapser
            # 2. fa2fq
            with xopen(self.output, 'wt') as fo:
                p0 = subprocess.Popen(self.cmd_open, shell=True, stdout=subprocess.PIPE)
                p1 = subprocess.Popen(cmd, shell=True, stdin=p0.stdout, stdout=subprocess.PIPE)
                name = None
                for line in p1.stdout:
                    line = line.decode() # bytes to str
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
                re_code = p1.communicate()

        return self.output


    def revcomp(self):
        """Reverse complement sequences"""
        revcomp = which('fastx_reverse_complement')
        if revcomp is None:
            raise Exception('command not found in $PATH: fastx_reverse_complement')

        cmd = revcomp
        if self.fx_type == 'fastq':
            cmd += ' -Q33'

        with xopen(self.input, 'rb') as fi, xopen(self.output, 'w') as fo:
            p = subprocess.Popen(cmd, shell=True, stdin=fi, stdout=subprocess.PIPE)
            for line in p.stdout:
                fo.write(line)
            p.communicate()

        return self.output


    def sample(self, n=1000, p=0.01):
        """Extract sample size of fastx file
        using seqkit tool

        $ seqkit sample -n 1000

        $ seqkit sample -p 0.01

        """
        assert isinstance(n, int)
        assert isinstance(p, float)

        seqkit = which('seqkit')
        if seqkit is None:
            raise Exception('command not found in $PATH: seqkit')

        logging.info('extracting sub-sample: n=%d' % n)

        ## warning
        if n > 1000000 or p > 0.1:
            logging.warning('choose a smaller number.')

        cmd = '%s sample -n %s' % (seqkit, n)

        with xopen(self.input, 'rb') as fi, xopen(self.output, 'w') as fo:
            p = subprocess.Popen(cmd, shell=True, stdin=fi, stdout=subprocess.PIPE)
            for line in p.stdout:
                fo.write(line)
            p.communicate()

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
        pass


class Cutadapt(object):
    """Cut adapters using cutadapt
    1. cut 3-ad
    2. cut 5-ad
    3. save log

    both SE and PE reads

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
            fq2_type = seq_type(fq_in)
            ## file type
            if fq2_type is None:
                raise Exception('unknown filetype, only support: [fasta, fastq]')
            ## file name
            fa2_clean_prefix = os.path.join(args['path_out'], fq2_prefix)
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
            AD3_list = [args['AD3']]
            ## short version of adapter
            if args['adapter_sliding']:
                i_list = []
                for i in AD3_list:
                    i_list.extend(self.adapter_chopper(i))
                AD3_list = i_list

            args_AD3 = ' '.join(['-A ' + i for i in AD3_list])

            ## option-2: filter
            ## 2.1 untrim 
            args_untrim = '-o %s --paired-output=%s' % (fq1_clean, fq2_clean)
            if args['save_untrim']:
                args_untrim += ' --untrimmed-output=%s \
                    --untrimmed-paired-output=%s --cores=1' % (fq1_untrim, fq2_untrim)
            elif args['rm_untrim']:
                args_untrim += ' --discard-trimmed --cores=%s' % args['threads']
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
        logging.info('trimming SE reads: %s' % fq_prefix)

        cmd = self.get_cutadapt_cmd()
        ## gzip file
        if args['gzip']:
            fq_clean += '.gz'
        ## run program
        if os.path.exists(fq_clean) and args['overwrite'] is False:
            logging.info('file exists, cutadapt skipped: %s' % fq_prefix)
        else:
            with xopen(fq_clean, 'wb') as fo, open(fq_log, 'wt') as flog:
                p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=flog)
                for line in p.stdout:
                    fo.write(line)
                p.communicate()
            Cutadapt_log(fq_log).saveas()
        return fq_clean


    def trim_pe(self):
        """Trimming Paired-end reads using cutadapt
        example: 
        cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq
        """
        args = self.kwargs.copy()
        fq1_prefix, fq1_clean, fq2_clean, fq1_log = self.cutadapt_init()[:4]
        logging.info('trimming PE reads: %s' % fq1_prefix)

        cmd = self.get_cutadapt_cmd()
        ## check exists
        fq1_check = fq1_clean
        fq2_check = fq2_clean
        if args['gzip']:
            fq1_check += '.gz'
            fq2_check += '.gz'
        ## compress fastq file after cutadapt
        if os.path.exists(fq1_check) and os.path.exists(fq2_check) and args['overwrite'] is False:
            logging.info('file exists, cutadapt skipped: %s' % fq1_prefix)
        else:
            with open(fq1_log, 'wt') as fo:
                p1 = subprocess.run(shlex.split(arg_cmd), stdout=fo, stderr=fo)
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
            logging.info('arguments not changed, trimming skipped - %s' % fq_prefix)
            ## SE mode
            if args['fq2'] is None:
                return_fq = fq1_clean
                if args['gzip']:
                    return_fq += '.gz'
            ## PE mode
            else:
                fq2_clean = args_init[2]
                return_fq = [fq1_clean, fq2_clean]
                if args['gzip']:
                    return_fq = [fq1_clean + '.gz', fq2_clean + '.gz']
        else:
            ## run program
            args_logger(args, args_file, overwrite=True) # update arguments.txt
            ## SE mode
            if args['fq2'] is None:
                fq1_clean = self.trim_se() # fq1
                return_fq = fq1_clean
                ## write to gzip file using gzip.open()
                # ## gzip compress files
                # if args['gzip']:
                #     fq1_clean = gzip_file(return_fq, rm=True)
                #     return_fq = fq1_clean
            else:
                return_fq = self.trim_pe() # fq1, fq2
                ## gzip compress files
                if args['gzip']:
                    fq1_clean = gzip_file(return_fq[0], rm=True)
                    fq2_clean = gzip_file(return_fq[1], rm=True)
                    return_fq = [fq1_clean, fq2_clean]

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


    def trimmer(self):
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
            logging.info('file exists, arguments not changed, trimming skipped : %s' % fq1_prefix)
            return return_output
        else:
            args_logger(args, args_file, True) # update arguments.txt

        ## update arguments for Cutadapt
        path_out = args.pop('path_out')
        fq1 = args.pop('fq1')
        adapter3 = args.pop('adapter3')
        len_min = args.pop('len_min')
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
            return_trim = Cutadapt(fq1, adapter3, path_cutadapt, len_min, **args).run()
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
            logging.info('[remving temp files]')
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
#     input: fq, adapter3, path_out, len_min, qual_min, error_rate, overlap, 
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

#     ## Collapse
#     input: fq
#     output: fq

#     """

#     def __init__(self, fq1, adapter3, path_out=None, len_min=15, **kwargs):
#         """Parsing the parameters for reads trimming
#         support both SE and PE reads
#         """
#         assert isinstance(fq1, str) # process one file
#         args1 = args_init(kwargs) # default parameters
#         args2 = {
#             'fq1': fq1,
#             'adapter3': adapter3,
#             'path_out': path_out,
#             'len_min': len_min}
#         args = {**args1, **args2}

#         ## validate options
#         assert is_path(path_out)
#         if args['cut_to_length'] > 0 and args['cut_to_length'] < len_min:
#             raise ValueError('[fatal] --trim-to-length [%s] shorter than -m [%s]') % (args['cut_to_length'], len_min)

#         self.kwargs = args # global 


#     def trim_init(self, fq_in=None):
#         """Prepare directory for trimming, filename for each file"""
#         args = self.kwargs.copy()
#         if fq_in is None:
#             fq_in = args['fq1']

#         # SE mode
#         if args['fq2'] is None:
#             fq_prefix = file_prefix(fq_in)[0]
#             fq_type = seq_type(fq_in)
#             if args['keep_name']:
#                 fq_clean = os.path.join(args['path_out'], '%s.fastq' % fq_prefix)
#             else:
#                 fq_clean = os.path.join(args['path_out'], '%s.clean.fastq' % fq_prefix)
#             fq_log = os.path.join(args['path_out'], '%s.cutadapt.log' % fq_prefix)
#             fq_untrim = os.path.join(args['path_out'], '%s.untrim.fastq' % fq_prefix)
#             return [fq_prefix, fq_clean, fq_log, fq_untrim]

#         # PE mode
#         elif os.path.exists(args['fq2']):
#             fq1_name = file_prefix(args['fq1'])[0]
#             fq2_name = file_prefix(args['fq2'])[0]
#             if args['keep_name']:
#                 fq1_clean = os.path.join(args['path_out'], '%s.fastq' % fq1_name)
#                 fq2_clean = os.path.join(args['path_out'], '%s.fastq' % fq2_name)
#             else:
#                 fq1_clean = os.path.join(args['path_out'], '%s.clean.fastq' % fq1_name)
#                 fq2_clean = os.path.join(args['path_out'], '%s.clean.fastq' % fq2_name)
#             fq_log = os.path.join(args['path_out'], '%s.cutadapt.log' % fq1_name)
#             fq1_untrim = os.path.join(args['path_out'], '%s.untrim.fastq' % fq1_name)
#             fq2_untrim = os.path.join(args['path_out'], '%s.untrim.fastq' % fq2_name)
#             return [fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim]
#         else:
#             raise ValueError('checkout -fq2 option')


#     def adapter_chopper(self, ad=None, step=2, window=15):
#         """Chop adapter by given length and step, to create a series of adapters"""
#         args = self.kwargs.copy()
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


#     def cut_parser(self, s, cut_para=True):
#         """Cut N-bases, before adapter trimming
#         parse the argument
#         if cut_para is True,
#         return the cut numbers in --cut {arg} version

#         else
#         return the numbers        
#         """
#         if ',' in s:
#             n = s.split(',')
#             if len(n) > 2:
#                 raise ValueError('illegal argument: %s' % s)
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


#     def get_cutadapt_cmd(self):
#         """Parse arguments, create cutadapt command line"""
#         args = self.kwargs.copy()

#         ## command
#         cut_exe = which('cutadapt')

#         ## determine the minimum length of reads
#         ## consider cut-after-trim
#         trim_args = self.cut_parser(args['cut_after_trim'], False) # return the numbers
#         if len(trim_args) == 2:
#             trim_5, trim_3 = trim_args[:2]
#         elif len(trim_args) == 1:
#             trim_5 = 0 if(int(trim_args[0]) < 0) else trim_args[0]
#             trim_3 = 0 if(int(trim_args[0]) > 0) else trim_args[0]
#         else:
#             trim_5 = trim_3 = 0
#         len_min = args['len_min'] + abs(trim_5) + abs(trim_3) # minimum lenght

#         ## basic command
#         cut_basic = '%s -m %s \
#                     -q %s \
#                     --overlap=%s \
#                     --error-rate=%s \
#                     --trim-n \
#                     --max-n=0.1 \
#                     --times=%s' % (
#                     cut_exe, 
#                     len_min, 
#                     args['qual_min'], 
#                     args['overlap'], 
#                     args['error_rate'], 
#                     args['trim_times'])

#         ## adapter3
#         ad3_list = args['adapter3'] #[args['adapter3'], ]
#         ## add short version of adapter3
#         if args['adapter_sliding']:
#             ad3_list = self.adapter_chopper(args['adapter3'])
#         ad3_arg = ' '.join(['-a %s' % i for i in ad3_list])

#         ## SE mode
#         if args['fq2'] is None:
#             fq_prefix, fq_clean, fq_log, fq_untrim = self.trim_init()

#             ## untrim
#             if args['rm_untrim']:
#                 untrim_arg = '--untrimmed-output=%s --cores=1 -o %s' % (fq_untrim, fq_clean)
#             else:
#                 untrim_arg = '--cores=%s -o %s' % (args['threads'], fq_clean)

#             ## cut before
#             if args['cut_before_trim'] > '0':
#                 cut_arg ='--cut %s'  % args['cut_before_trim']
#             else:
#                 cut_arg = ''

#             arg_cmd = ' '.join([cut_basic, ad3_arg, untrim_arg, cut_arg, args['fq1']])

#         ## PE mode
#         elif os.path.exists(args['fq2']):
#             fq_prefix, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim = self.trim_init()

#             AD3_list = [args['AD3'], ]
#             if args['adapter_sliding']:
#                 AD3_list = self.adapter_chopper(args['AD3'])
#             AD3_arg = ' '.join(['-A %s' % i for i in AD3_list])

#             ## untrim
#             if args['rm_untrim']:
#                 untrim_arg = '--untrimmed-output=%s --untrimmed-paired-output=%s \
#                     --cores=1 -o %s --paired-output=%s' % (fq1_untrim, fq2_untrim,
#                         fq1_clean, fq2_clean)
#             else:
#                 untrim_arg = '--cores=%s -o %s --paired-output=%s' % (args['threads'],
#                         fq1_clean, fq2_clean)

#             ## cut before
#             if args['cut_before_trim'] > '0':
#                 cut_arg ='--cut %s'  % args['cut_before_trim']
#             else:
#                 cut_arg = ''

#             ## merge
#             arg_cmd = ' '.join([cut_basic, ad3_arg, AD3_arg, untrim_arg, cut_arg, args['fq1'], args['fq2']])
#         else:
#             raise ValueError('--fq2 value illegal')

#         return arg_cmd


#     def trim_se(self):
#         """Trimming SE reads using cutadapt"""
#         args = self.kwargs.copy()
#         fq_prefix, fq_clean, fq_log, fq_untrim = self.trim_init()
#         logging.info('trimming SE reads: %s' % fq_prefix)

#         arg_cmd = self.get_cutadapt_cmd()
#         if args['gzipped']:
#             fq_clean = fq_clean + '.gz'

#         if os.path.exists(fq_clean) and args['overwrite'] is False:
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
#         args = self.kwargs.copy()
#         fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim = self.trim_init()
#         logging.info('trimming PE reads: %s' % fq1_name)

#         arg_cmd = self.get_cutadapt_cmd()
#         if os.path.exists(fq1_clean) and os.path.exists(fq2_clean) and args['overwrite'] is False:
#             logging.info('file exists, cutadapt skipped: %s' % fq1_name)
#         else:
#             with open(fq_log, 'wt') as fo:
#                 p1 = subprocess.run(shlex.split(arg_cmd), stdout=fo, stderr=fo)
#             # Cutadapt_log(fq_log).saveas()
#         return [fq1_clean, fq2_clean]


#     def rm_duplicate(self, fq_in=None, fq_out=None):
#         """Remove possible PCR duplicates, using fastx_collapse
#         collapse reads before trimming
#         ****
#         only support SE reads (not PE reads)
#         """
#         logging.info('removing PCR duplicates')
#         args = self.kwargs.copy()
#         path_rmdup = os.path.join(args['path_out'], 'rm_PCR_dup')
#         assert is_path(path_rmdup)
#         basedir = os.path.dirname(os.path.realpath(__file__))
#         fa2fq = os.path.join(basedir, 'fasta_to_fastq.pl')

#         if not os.path.exists(fa2fq):
#             raise Exception('file not exists - %s' % fa2fq)

#         if fq_in is None:
#             fq_in = args['fq1']
#         if fq_out is None:
#             fq_out = os.path.join(path_rmdup, os.path.basename(fq_in))
#         fq_prefix, fq_clean, fq_log, fq_untrim = self.trim_init(fq_in)

#         fastx_collapser = which('fastx_collapser')

#         if os.path.exists(fq_out) and args['overwrite'] is False:
#             logging.info('file exists, skip rm_dup: %s' % fq_prefix)
#         else:
#             # create new clean fastq file
#             with open(fq_out, 'wt') as ff:
#                 c1 = '%s -Q33 -i %s' % (fastx_collapser, fq_in)
#                 c2 = 'perl %s -' % fa2fq
#                 p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE)
#                 p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, stdout=ff)
#                 px = p2.communicate()
#         return fq_out


#     def cut_ends(self, fq_in=None, fq_out=None, rm_too_short=True):
#         """Cut N-bases at either tail of fastq after cutadapt
#         Using python to process fastq files
#         support SE reads only
#         """
#         logging.info('Trimming ends of fastq file')
#         args = self.kwargs.copy()
#         path_cut_ends = os.path.join(args['path_out'], 'cut_ends')
#         assert is_path(path_cut_ends)

#         if fq_in is None:
#             raise ValueError('either fq_in or fq_out are required')
#         if fq_out is None:
#             fq_out = os.path.join(path_cut_ends, os.path.basename(fq_in))

#         args_cut_ends = self.cut_parser(args['cut_after_trim'], False)

#         if len(args_cut_ends) == 2:
#             trim_5, trim_3 = args_cut_ends
#         elif len(args_cut_ends) == 1:
#             trim_5 = 0 if(int(args_cut_ends[0]) < 0) else args_cut_ends[0]
#             trim_3 = 0 if(int(args_cut_ends[0]) > 0) else args_cut_ends[0]
#         else:
#             trim_5 = trim_3 = 0

#         with open(fq_in, 'rt') as fi, open(fq_out, 'wt') as ff:
#             while True:
#                 try:
#                     fq_id, fq_seq, fq_plus, fq_qual = [next(fi).strip(),
#                                                        next(fi).strip(),
#                                                        next(fi).strip(),
#                                                        next(fi).strip()]
#                     if len(fq_seq) < args['len_min'] + abs(trim_5) + abs(trim_3):
#                         if rm_too_short:
#                             continue # skip short reads
#                         else:
#                             fq_seq = 'A'
#                             fq_qual = 'J' # quality
#                     else:
#                         fq_seq = fq_seq[trim_5:trim_3] if(trim_3 < 0) else fq_seq[trim_5:]
#                         fq_qual = fq_qual[trim_5:trim_3] if(trim_3 < 0) else fq_qual[trim_5:]
#                     # trim to length
#                     if args['cut_to_length'] > args['len_min']:
#                         n_right = min([args['cut_to_length'], len(fq_seq)])
#                         fq_seq = fq_seq[:n_right]
#                         fq_qual = fq_qual[:n_right]
#                     ff.write('\n'.join([fq_id, fq_seq, fq_plus, fq_qual]) + '\n')
#                 except StopIteration:
#                     break

#         return fq_out


#     def cut_ends2(self, fq_in=None, fq_out=None):
#         """Cut reads from 3' end to specific length, 
#         save reads with maximum length
#         """
#         logging.info('Trimming reads to specific length')
#         args = self.kwargs.copy()
#         path_cut_ends = os.path.join(args['path_out'], 'cut_ends2')
#         assert is_path(path_cut_ends)

#         if fq_in is None:
#             fq_in = args['fq1']
#         if fq_out is None:
#             fq_out = os.path.join(path_cut_ends, os.path.basename(fq_in))
#         fq_prefix, fq_clean, fq_log, fq_untrim = self.trim_init(fq_in)

#         with open(fq_in, 'rt') as fi, open(fq_out, 'wt') as ff:
#             while True:
#                 try:
#                     fq_id, fq_seq, fq_plus, fq_qual = [next(fi).strip(),
#                                                        next(fi).strip(),
#                                                        next(fi).strip(),
#                                                        next(fi).strip()]
#                     n_right = min([args['cut_to_length'], len(fq_seq)])
#                     fq_seq = fq_seq[:n_right]
#                     fq_qual = fq_qual[:n_right]
#                     ff.write('\n'.join([fq_id, fq_seq, fq_plus, fq_qual]) + '\n')
#                 except StopIteration:
#                     break
#         return fq_out


#     def run_se(self):
#         """Run cutadapt for SE reads
#         1. cut N-bases before trimming (optional)
#         2. trim 3' adapter
#         3. remove PCR duplicates (optional)
#         4. cut N-bases after trimming (optional)
#         5. cut reads into specific length (optional)
#         """
#         logging.info('Trimming reads: SE mode')
#         args = self.kwargs.copy()

#         fq_prefix, fq_clean, fq_log, fq_untrim = self.trim_init()
#         fq_in = args['fq1']

#         ## 1. cut N-bases before trim
#         ## 2. trim 3' adapter
#         fq_out = self.trim_se() # fq_out == fq_clean

#         ## 3. remove PCR dup
#         if args['rm_dup']:
#             fq_tmp1 = fq_in + '.before_rmdup.tmp'
#             os.rename(fq_clean, fq_tmp1)
#             fq_clean = self.rm_duplicate(fq_in=fq_tmp1, fq_out=fq_clean)
#             os.remove(fq_tmp1)

#         ## 4. cut N-bases after trim
#         if not args['cut_after_trim'] == '0':
#             fq_tmp2 = fq_in + '.before_cut_after_trim.tmp'
#             os.rename(fq_clean, fq_tmp2)
#             fq_clean = self.cut_ends(fq_in=fq_tmp2, fq_out=fq_clean)
#             os.remove(fq_tmp2)

#         ## 5. cut to specific length
#         if args['cut_to_length'] > args['len_min']:
#             fq_tmp3 = fq_in + '.before_cut_to_length.tmp'
#             os.rename(fq_clean, fq_tmp3)
#             fq_clean = self.cut_ends2(fq_in=fq_tmp3, fq_out=fq_clean)
#             os.remove(fq_tmp3)

#         return fq_clean


#     def run_pe(self):
#         """Run cutadapt for PE reads
#         1. cut N-bases before trimming (optional)
#         2. trim 3' adapter
#         2. cut N-bases after trimming (optional)
#         """
#         logging.info('Trimming reads: PE mode')
#         args = self.kwargs.copy()
#         fq1_name, fq1_clean, fq2_clean, fq_log, fq1_untrim, fq2_untrim = self.trim_init()
#         fq_in = args['fq1']

#         ## 1. cut N-bases before trim
#         ## 2. trim 3' adapter
#         fq1_clean, fq2_clean = self.trim_pe() # fq_out == fq_clean

#         ## 3. cut N-bases after trim
#         if not args['cut_after_trim'] == '0':
#             fq1_tmp1 = fq1_clean + '.before_cut_after_trim.tmp'
#             fq2_tmp1 = fq2_clean + '.before_cut_after_trim.tmp'
#             os.rename(fq1_clean, fq1_tmp1)
#             os.rename(fq2_clean, fq2_tmp1)
#             fq1_clean = self.cut_ends(fq_in=fq1_tmp1, fq_out=fq1_clean)
#             fq2_clean = self.cut_ends(fq_in=fq1_tmp1, fq_out=fq2_clean)
#             os.remove(fq1_tmp1)
#             os.remove(fq2_tmp1)

#         ## 4. cut to specific length
#         if args['cut_to_length'] > args['len_min']:
#             fq1_tmp2 = fq1_clean + '.before_cut_to_length.tmp'
#             fq2_tmp2 = fq2_clean + '.before_cut_to_length,tmp'
#             os.rename(fq1_clean, fq1_tmp2)
#             os.rename(fq2_clean, fq2_tmp2)
#             fq1_clean = self.cut_ends2(fq_in=fq1_tmp2, fq_out=fq1_clean)
#             fq2_clean = self.cut_ends2(fq_in=fq2_tmp2, fq_out=fq2_clean)

#         return [fq1_clean, fq2_clean]


#     def run(self):
#         """Run trimming for all"""
#         args = self.kwargs.copy()
#         fq_prefix = file_prefix(args['fq1'])[0]

#         ## save arguments
#         args_file = os.path.join(args['path_out'], fq_prefix + '.arguments.txt')
#         args_pickle = os.path.join(args['path_out'], fq_prefix + '.arguments.pickle')
#         if args_checker(args, args_pickle) and args['overwrite'] is False:
#             logging.info('files exists, arguments not changed, trimming skipped - %s' % fq_prefix)
#             return True #!!! fastq files
#         else:
#             args_logger(args, args_file, True) # update arguments.txt

#         ## SE mode
#         if args['fq2'] is None:
#             fq_return = self.run_se()
#             fq_return = gzip_file(fq_return, rm=True)
#             fq1_clean = fq_return

#         ## PE mode
#         elif os.path.exists(args['fq2']):
#             fq_return = self.run_pe()
#             fq_return1 = gzip_file(fq_return[0], rm=True)
#             fq_return2 = gzip_file(fq_return[1], rm=True)
#             fq1_clean = fq_return1

#         else:
#             raise Exception('Illegal --fq2 argument: %s' % args['fq2'])

#         ## save read count
#         fq_prefix = file_prefix(args['fq1'])[0]
#         fq_count_file = os.path.join(args['path_out'], fq_prefix + '.clean_reads.txt')
#         if not os.path.exists(fq_count_file) or args['overwrite'] is True:
#             fq_count = int(file_row_counter(fq1_clean) / 4)
#             with open(fq_count_file, 'wt') as fo:
#                 fo.write(str(fq_count) + '\n')

#         return fq_return


## EOF
