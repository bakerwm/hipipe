#!/usr/bin/env python3

import os
import logging
import fnmatch
import subprocess
import sys

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)


class QC_reporter(object):
    """Make html report for fastqc statistics file
    input fastq file or directory contains fastq file
    output output directory
    template (optional) customized template
    """

    def __init__(self, input, output, template=None):
        if isinstance(input, list):
            self.input = input
        else:
            self.input = [input, ]
        self.output = output
        self.template = template


    def findfiles(self, which, where='.'):
        """Returns list of filenames from `where` path matched by 'which'
        shell pattern. Matching is case-insensitive.
        # findfiles('*.ogg')
        """    
        # TODO: recursive param with walk() filtering
        return [os.path.join(where, f) for f in os.listdir(where) if fnmatch.fnmatch(f, which)]


    def get_fq_files(self, x):
        """
        if x is directory, list all fastq files within x directory
        if x is a file, return the file
        fastq file format:
        *.fq, *.fastq, *.fq.gz, *.fastq.gz
        """
        fq_files = []
        if os.path.isdir(x):
            fq_files = self.findfiles("*.f[astq]*", x)
        elif os.path.isfile(x):
            # f1 = [x]
            fq_files = [x]
        else:
            # f1 = None
            fq_files = None

        return fq_files


    def run_shell_cmd(self, cmd): 
        """Run shell command"""
        try:
            p = subprocess.run(cmd, shell=True,
                stdout=subprocess.PIPE,
                # stderr=subprocess.STDOUT,
                universal_newlines=True,
                preexec_fn=os.setsid)
            logging.info('run_shell_cmd: CMD={}'.format(cmd))
        except:
             raise Exception('Killed CMD={}\nSTDOUT={}'.format(
                 cmd, ret))


    def check_fastqc(self, path, fqs, overwrite=False):
        """Check the fastqc output files exists in path
        input: <prefix>.fq.gz
        output: <prefix>_fastqc.html
                <prefix>_fastqc.html
        """
        fq_tmp = []

        for fq in sorted(fqs):
            prefix = os.path.basename(os.path.splitext(fq)[0])
            if fq.endswith('.gz'):
                prefix = os.path.splitext(prefix)[0]
            f1 = os.path.join(path, prefix + '_fastqc.html')
            f2 = os.path.join(path, prefix + '_fastqc.zip')
            if os.path.exists(f1) and os.path.exists(f2) and overwrite is False:
                logging.info('file exists, fastqc skipped - %s' % prefix)
                continue
            else:
                logging.info('run fastqc - %s' % prefix)
                fq_tmp.append(fq)
        return fq_tmp


    def run(self):
        """Run fastqc and fastqc report"""
        fq_files = []
        # list fastq files
        for x in self.input:
            q1 = self.get_fq_files(x)
            fq_files.extend(q1)

        if not os.path.exists(self.output):
            os.makedirs(self.output)

        # check fastqc output
        logging.info('fastqc')
        fq_files = self.check_fastqc(self.output, fq_files, overwrite=False)
        
        # run fastqc
        fq_list = ' '.join(fq_files)
        cmd1 = 'fastqc -o %s %s' % (self.output, fq_list)
        if len(fq_files) > 0:
            self.run_shell_cmd(cmd1)

        # fastqc report
        logging.info('fastqc-report')
        main_script = os.path.realpath(__file__)
        home = os.path.dirname(main_script)
        report_r = os.path.join(home, 'qc_report.R')
        report_html = os.path.join(self.output, 'fastqc_report.html')
        ## for update template
        if self.template is None:
            self.template = ''
        elif os.path.exists(self.template):
            self.template = os.path.abspath(self.template)
        else:
            self.template = ''

        cmd2 = 'Rscript %s %s %s' % (report_r, self.output, self.template)

        if os.path.exists(report_html):
            logging.info('file exists, fastqc-report skipped')
        else:
            self.run_shell_cmd(cmd2)

        if not os.path.exists(report_html):
            logging.error('failed, generating html file')

        # finish
        logging.info('output - %s' % report_html)


class Alignment_reporter(object):
    """Make html report for alignment statistics file
    input map_sta file or directory
    output html report
    """

    def __init__(self, input, output, 
        suffix='*.mapping_stat.csv',
        stat_list_name='mapping_stat.list', 
        out_html_name='alignment_report.html',
        template=None):
        if isinstance(input, list):
            self.input = input
        else:
            self.input = [input, ]
        self.output = output
        self.template = template
        self.suffix = suffix
        self.stat_list = os.path.join(output, stat_list_name)
        self.out_html = os.path.join(output, out_html_name)


    def findfiles(self, which, where='.'):
        """Returns list of filenames from `where` path matched by 'which'
        shell pattern. Matching is case-insensitive.
        # findfiles('*.ogg')
        """    
        # TODO: recursive param with walk() filtering
        return [os.path.join(where, f) for f in os.listdir(where) if fnmatch.fnmatch(f, which)]


    def get_stat_files(self):
        """Search the stat files, or return the directory
        match the filename: *mapping_stat.csv
        """
        fs = []

        for x in self.input:
            if os.path.isdir(x):
                tmp1 = self.findfiles(self.suffix, x)
                fs.extend(tmp1)
            elif os.path.isfile(x):
                fs = [x]
            else:
                continue

        # convert to absolute path
        fs = [os.path.abspath(f) for f in fs if not f is None]

        return fs


    def save_stat_files(self, x, t=None):
        """Write list of files to file t"""
        if isinstance(x, str):
            x_list = [x, ]
        elif isinstance(x, list):
            x_list = x
        else:
            raise Exception('character or list supported, unknown types detected')

        if t is None:
            list_file = self.stat_list
        else:
            list_file = t

        with open(list_file, 'wt') as fo:
            fo.write('\n'.join(x_list) + '\n')

        return list_file


    def run_shell_cmd(self, cmd): 
        """Run shell command"""
        try:
            p = subprocess.run(cmd, shell=True,
                stdout=subprocess.PIPE,
                #stderr=subprocess.STDOUT,
                universal_newlines=True,
                preexec_fn=os.setsid)
            logging.info('run_shell_cmd: CMD={}'.format(cmd))
        except:
             raise Exception('Killed CMD={}\nSTDOUT={}'.format(
                 cmd, ret))


    def run(self):
        """Run fastqc and fastqc report"""
        ## input
        input_files = self.get_stat_files()
        if len(input_files) == 0:
            raise Exception(self.suffix + ' file not detected in ' + self.output)

        ## check template
        if self.template is None:
            self.template = ' '
        elif os.path.exists(self.template):
            self.template = os.path.abspath(self.template)
        else:
            self.template = ' '

        ## output
        if not os.path.exists(self.output):
            os.makedirs(self.output)
        stat_list = self.save_stat_files(input_files)
        report_html = self.out_html

        ## run
        logging.info('alignment-report')
        main_script = os.path.realpath(__file__)
        home = os.path.dirname(main_script)
        align_stat_r = os.path.join(home, 'alignment_stat.R')
        cmd1 = '%s %s %s %s %s' % ('Rscript', align_stat_r, stat_list, 
            self.output, self.template)

        if os.path.exists(report_html):
            logging.info('file exists, alignment-report skipped')
        else:
            self.run_shell_cmd(cmd1)

        if not os.path.exists(report_html):
            logging.error('failed generating html file')

        # finish
        logging.info('saving results in - %s' % report_html)


