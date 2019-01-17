#!/usr/bin/env python3

import os
import re
import logging
import argparse
import fnmatch
import subprocess

class Alignment_reporter(object):
    """Make html report for alignment statistics file
    input map_sta file or directory
    output html report
    """

    def __init__(self, input, output, template=None):
        self.input = input
        self.output = output
        self.template = template


    def findfiles(self, which, where='.'):
        """Returns list of filenames from `where` path matched by 'which'
        shell pattern. Matching is case-insensitive.
        # findfiles('*.ogg')
        """    
        # TODO: recursive param with walk() filtering
        return [os.path.join(where, f) for f in os.listdir(where) if fnmatch.fnmatch(f, which)]


    def get_stat_files(self, x):
        """Search the stat files, or return the directory"""
        stat_files = []
        if os.path.isdir(x):
            tmp1 = self.findfiles("*.mapping_stat.csv", x)
            stat_files.extend(tmp1)
        elif os.path.isfile(x):
            stat_files = [x]
        else:
            stat_files = [None]

        # convert to absolute path
        stat_files = [os.path.abspath(f) for f in stat_files if not f is None]

        return stat_files


    def run_shell_cmd(self, cmd): 
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


    def run(self):
        """Run fastqc and fastqc report"""
        ## input
        input_files = []
        for i in self.input:
            fi = self.get_stat_files(i)
            if not fi is None:
                input_files.extend(fi)

        if len(input_files) == 0:
            raise Exception('*mapping_stat.csv file not detected')

        ## output
        if not os.path.exists(self.output):
            os.makedirs(self.output)

        ## stat_list to file
        list_file = os.path.join(self.output, 'mapping_stat.list')
        with open(list_file, 'wt') as fo:
            fo.write('\n'.join(input_files) + '\n')

        ## run
        logging.info('running alignment-report')
        main_script = os.path.realpath(__file__)
        home = os.path.dirname(main_script)
        align_stat_r = os.path.join(home, 'alignment_stat.R')
        align_stat_html = os.path.join(self.output, 'alignment_report.html')
        if self.template is None:
            self.template = ' '
        elif os.path.exists(self.template):
            self.template = os.path.abspath(self.template)
        else:
            self.template = ' '

        ## run command
        cmd1 = '%s %s %s %s %s' % ('Rscript', align_stat_r, list_file, self.output, self.template)

        if not os.path.exists(align_stat_html):
            self.run_shell_cmd(cmd1)

        # finish
        logging.info('saving results in - %s' % align_stat_html)

