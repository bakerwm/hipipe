#!/usr/bin/env python3
"""
Create fastqc report for multiple fastq files
"""

__author__ = 'Ming Wang'
__email__ = 'wangm08@hotmail.com'
__date__ = '2018-12-25'
__version__ = '0.2'

import os
import re
import logging
import argparse
import fnmatch
import subprocess

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

def get_args():
    """
    Mapping SE read or one of PE reads to reference genome
    using bowtie, STAR, ... (universal)
    """
    parser = argparse.ArgumentParser(prog='aligner', 
                                     description='mapping reads')
    parser.add_argument('-i', nargs='+', required=True, metavar='INPUT', 
        help='fastq files or directories contain fastq files')
    parser.add_argument('-o', default=None, 
        metavar='OUTPUT',  help='The directory to save results, default, \
        current working directory.')
    args = parser.parse_args()
    return args


def findfiles(which, where='.'):
    """Returns list of filenames from `where` path matched by 'which'
    shell pattern. Matching is case-insensitive.
    # findfiles('*.ogg')
    """    
    # TODO: recursive param with walk() filtering
    return [os.path.join(where, f) for f in os.listdir(where) if fnmatch.fnmatch(f, which)]


def get_fq_files(x):
    """
    if x is directory, list all fastq files within x directory
    if x is a file, return the file
    fastq file format:
    *.fq, *.fastq, *.fq.gz, *.fastq.gz
    """
    if os.path.isdir(x):
        f1 = findfiles("*.f[astq]*[.gz]*", "data/clean_data/")
    elif os.path.isfile(x):
        f1 = [x]
    else:
        f1 = None

    return f1


def run_shell_cmd(cmd): 
    """Run shell command"""
    try:
        p = subprocess.run(cmd, shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            preexec_fn=os.setsid)
        logging.info('run_shell_cmd: CMD={}'.format(cmd))
    except:
         raise Exception('Killed CMD={}\nSTDOUT={}'.format(
             cmd, ret))


def check_fastqc(path, fqs, overwrite=False):
    """Check the fastqc output files exists in path
    input: <prefix>.fq.gz
    output: <prefix>_fastqc.html
            <prefix>_fastqc.html
    """
    fq_tmp = []

    for fq in fqs:
        prefix = os.path.basename(os.path.splitext(fq)[0])
        if fq.endswith('.gz'):
            prefix = os.path.splitext(prefix)[0]
        f1 = os.path.join(path, prefix + '_fastqc.html')
        f2 = os.path.join(path, prefix + '_fastqc.zip')
        if os.path.exists(f1) and os.path.exists(f2) and overwrite is False:
            logging.info('file exists, fastqc skipped - %s' % prefix)
            continue
        else:
            fq_tmp.append(fq)
    return fq_tmp


def main():
    """Run fastqc and fastqc report"""
    args = get_args()

    fq_files = []
    # list fastq files
    for x in args.i:
        q1 = get_fq_files(x)
        fq_files.extend(q1)

    if not os.path.exists(args.o):
        os.makedirs(args.o)

    # check fastqc output
    fq_files = check_fastqc(args.o, fq_files, overwrite=False)
    
    # run fastqc
    logging.info('running FastQC')
    fq_list = ' '.join(fq_files)
    cmd1 = 'fastqc -o %s %s' % (args.o, fq_list)
    if len(fq_files) > 0:
        run_shell_cmd(cmd1)

    # fastqc report
    logging.info('running fastqc-report')
    main_script = os.path.realpath(__file__)
    home = os.path.dirname(main_script)
    report_r = os.path.join(home, 'qc_report.R')
    report_html = os.path.join(args.o, 'report.html')
    cmd2 = 'Rscript %s %s' % (report_r, args.o)
    if not os.path.exists(report_html):
        run_shell_cmd(cmd2)

    # finish
    logging.info('saving results in - %s' % report_html)

if __name__ == '__main__':
    main()


## EOF