#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MACS2 

calpeak
bdgcmp

"""

import os
import sys
import glob
import tempfile
import platform
import shlex
import subprocess
import logging
from helper import *



class Venv(object):
    """Check virtualenv created by virtualenv"""

    def __init__(self, venv):
        self.venv = venv


    # current, run
    def _base_prefix(self):
        if hasattr(sys, 'real_prefix'):
            return sys.real_prefix
        elif hasattr(sys, 'base_prefix'):
            return sys.base_prefix
        else:
            return None


    ## virtualenv 
    def _venv_status(self):
        """Determine if inside virtualenv"""
        return (hasattr(sys, 'real_prefix') or 
            (hasattr(sys, 'base_prefix') and sys.base_prefix !=  sys.prefix))


    ## enter env
    def _venv_switcher(self, venv):
        """Switch to specific virtualenv, venv"""
        venv_bin = os.path.join(venv, 'bin', 'activate_this.py')
        if sys.version_info[0:2] ==  (2, 7):
            execfile(venv_bin, dict(__file__ = venv_bin)) # python2      
        elif sys.version_info[0:2] >=  (3, 0):
            # execfile(venv_bin, dict(__file__ = venv_bin)) # python2
            exec(open(venv_bin).read(), {'__file__': venv_bin}) #python3
        else:
            logging.error('unknown version of python: ' + sys.version)


    def _venv_validator(self, into_venv=True):
        """Validate the virtualenv, and switch to it
        check virtualenv 
        if not, go into
        into_venv, into env or out env
        """
        venv_in = os.path.expanduser(self.venv)
        
        if into_venv: # switch to env
            if self._venv_status() and venv_in ==  self._base_prefix():
                logging.info('already in venv: %s' % venv_in)
            elif os.path.exists(venv_in):
                self._venv_switcher(venv_in)
            else:
                logging.error('virtualenv not exists: %s' % self.venv)
        else: # exit env
            if self._venv_status() and venv_in ==  sys.prefix:
                os.system('deactivate')
            else:
                logging.info('not in venv: ' + venv_in)


class Macs2(object):
    """Run macs2 for BAM files
    1. macs2 callpeak -f BAM -t {IP.bam} -c {input.bam} -g {gsize} --outdir {out_dir} 
        -n {prefix} -B --SPMR {--broad} {--keep-dup auto}
    2. macs2 bdgcmp -t {prefix}_treat_pileup.bdg -c {prefix}_control_lambda.bdg -o {prefix}.ppois.bdg -m ppois
    3. macs2 bdgcmp -t {prefix}_treat_pileup.bdg -c {prefix}_control_lambda.bdg -o {prefix}.FE.bdg -m FE
    4. macs2 bdgcmp -t {prefix}_treat_pileup.bdg -c {prefix}_control_lambda.bdg -o {prefix}.logLR.bdg -m logLR -p 0.00001
    5. 
    """

    def __init__(self, ip, control, genome, output, prefix=None, venv=None, overwrite=False):
        """Parse the parameters
        venv, the virtualenv created for macs2, running in Python2
        """
        self.ip = ip
        self.control = control
        self.genome = genome
        self.output = output
        self.overwrite = overwrite
        if prefix is None:
            prefix = file_prefix(ip)[0]
            # prefix = os.path.splitext(os.path.basename(ip))[0]
        self.prefix = prefix
        if venv is None:
            venv = '/home/wangming/envs/piPipes/'
            # function to create virtualenv
            # pass
        self.venv = venv
        is_path(self.output)


    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.out',
                                            delete=False)
        return tmpfn.name


    def get_gsize(self):
        """Return the genome size of the genome"""
        gsize_file = Genome(self.genome).get_fasize()
        gsize = 0
        with open(gsize_file, 'rt') as fi:
            for a in fi:
                c, n = a.strip().split('\t')
                gsize += int(n)
        return gsize


    def python2_run(self, cmd, log=None):
        """Run commands in python2 ENV"""
        ## log file
        if log is None:
            log = self._tmp()

        ## switch to python2,
        Venv(self.venv)._venv_validator(True)
        # if platform.python_version() > '3.0':
        #     raise ValueError('macs2 needs to be run in Python 2')
        # cmd = 'macs2 --help'
        # print(cmd)
        with open(log, 'wt') as fo:
            subprocess.call(shlex.split(cmd), stdout=fo, stderr=fo)


    def callpeak(self):
        """Call peaks using MACS"""
        genome_size = self.get_gsize()
        log = os.path.join(self.output, self.prefix + '.macs2.callpeak.out')
        macs2_cmd = "macs2 callpeak -f BAM -t %s -c %s -g %s --outdir %s -n %s --broad \
            --keep-dup auto -B --SPMR" % (self.ip, self.control, genome_size,
            self.output, self.prefix)
        # output file
        macs2_out = os.path.join(self.output, self.prefix + '_peaks.xls')
        if os.path.exists(macs2_out) and self.overwrite is False:
            logging.info('file exists, skip macs2 callpeak')
        else:
            logging.info('run macs2 callpeak')
            self.python2_run(macs2_cmd, log)



    def bdgcmp(self, opt='ppois'):
        """Options for -m:
        {ppois,qpois,subtract,logFE,FE,logLR,slogLR,max}
        """
        # use output of callpeak
        ip_bdg = os.path.join(self.output, self.prefix + '_treat_pileup.bdg')
        input_bdg = os.path.join(self.output, self.prefix + '_control_lambda.bdg')

        if not os.path.exists(ip_bdg) or not os.path.exists(input_bdg):
            raise ValueError('*.bdg file not found, need to run .callpeak() first')
        if not opt in ['ppois', 'qpois', 'subtract', 'logFE', 'FE', 'logLR', 'slogLR', 'max']:
            raise ValueError('unknown option: opt=%s' % opt)

        if opt == 'logLR':
            opt_ext = '-p 0.00001'
        else:
            opt_ext = ''

        out_bdg = os.path.join(self.output, self.prefix + '.' + opt + '.bdg')
        log = os.path.join(self.output, self.prefix + '.macs2.bdgcmp.out')
        c = "macs2 bdgcmp -t %s -c %s -o %s -m %s %s" % (ip_bdg, input_bdg, out_bdg, opt, opt_ext)
        
        if os.path.exists(out_bdg) and self.overwrite is False:
            logging.info('file exists, skip macs2 bdgcmp')
        else:
            self.python2_run(c)

        # sort output *.bdg
        c1 = 'sort -k1,1 -k2,2n -o %s %s' % (out_bdg, out_bdg)

        # cnvert *.bdg to *.bigWig
        gsize_file = Genome(self.genome).get_fasize()
        out_bw  = os.path.join(self.output, self.prefix + '.' + opt + '.bigWig')
        c2 = 'bedGraphToBigWig %s %s %s' % (out_bdg, gsize_file, out_bw)
        if os.path.exists(out_bw) and self.overwrite is False:
            logging.info('file exists, skip bg2bw')
        else:
            subprocess.run(shlex.split(c1)) # sort bdg
            subprocess.run(shlex.split(c2)) # convert bdg to bw


    def bdgpeakcall(self):
        pass


    def bdgopt(self):
        pass


    def cmbreps(self):
        pass


    def bdgdiff(self):
        pass


    def filterdup(self):
        pass


    def predicted(self):
        pass


    def pileup(self):
        pass


    def randsample(self):
        pass


    def refinepeak(self):
        pass


    def broadpeak_annotation(self):
        """Annotate broadpeak using HOMER annotatePeaks.pl script
        BED foramt
        """
        anno_exe = which('annotatePeaks.pl')
        if not os.path.exists(anno_exe):
            logging.error('command not exists, skip annotation - %s' % anno_exe)
            return None

        # broad peak file
        broadpeak = os.path.join(self.output, self.prefix + '_peaks.broadPeak')
        if not os.path.exists(broadpeak):
            raise ValueError('file not found, need to run .callpeak() first - %s' % broadpeak)

        # run
        anno_peak = os.path.join(self.output, self.prefix + '_peaks.broadPeak.annotation')
        anno_log = os.path.join(self.output, self.prefix + '_peaks.broadPeak.annotation.log')
        cmd = 'perl %s %s %s' % (anno_exe, broadpeak, self.genome)
        with open(anno_peak, 'wt') as fo, open(anno_log, 'wt') as fe:
            subprocess.run(shlex.split(cmd), stdout=fo, stderr=fe)
        return anno_peak


    def get_effect_size(self):
        """Extract the effective depth of macs2 files
        parse the file: output/*_peaks.xls
        tags after filtering in treatment
        tags in treatment
        tags after filtering in control
        tags in control
        # tag size is determined as 100 bps
        # total tags in treatment: 7978071
        # tags after filtering in treatment: 2384854
        # maximum duplicate tags at the same position in treatment = 1
        # Redundant rate in treatment: 0.70
        # total tags in control: 10555283
        # tags after filtering in control: 6639591
        # maximum duplicate tags at the same position in control = 1
        # Redundant rate in control: 0.37
        # d = 122
        # alternative fragment length(s) may be 122 bps
        """

        # search the xls file
        f = glob.glob(os.path.join(self.output, '*_peaks.xls'))
        if not os.path.exists(f[0]):
            raise ValueError('file missing in macs2 callpeak output: %s' % f)

        # top
        topN = 100
        counter = 0
        dep = {}
        # ip_depth = ip_scale = input_depth = input_scale = 0
        with open(f[0], 'rt') as fi:
            for line in fi:
                if not line.startswith('#'): 
                    continue
                if counter > 100: # nrows
                    break # stop
                num = line.strip().split()[-1]
                if 'tags after filtering in treatment' in line:
                    dep['ip_depth'] = num
                if 'tags in treatment' in line:
                    s = 1e6 / int(num)
                    dep['ip_scale'] = '%.6f' % s
                if 'tags after filtering in control' in line:
                    dep['input_depth'] = num
                if 'tags in control' in line:
                    s = 1e6 / int(num)
                    dep['input_scale'] = '%.6f' % s
                counter += 1

        return dep




# ctl = '/home/wangming/work/yu_2018/projects/20181106_YY25_zk05_chipseq/results/piPipies_output/97FL_Nxf2MT_rep1/genome_mapping_unique_only/97FL_Nxf2MT_H3K9me3_rep1_1.clean.dm3.Input.sorted.bam'
# tre = '/home/wangming/work/yu_2018/projects/20181106_YY25_zk05_chipseq/results/piPipies_output/97FL_Nxf2MT_rep1/genome_mapping_unique_only/97FL_Nxf2MT_H3K9me3_rep1_1.clean.dm3.IP.sorted.bam'
# out = '/home/wangming/work/yu_2018/projects/20181106_YY25_zk05_chipseq/results/piPipies_output/97FL_Nxf2MT_rep1/aaaaaa'
# prefix = 'aaaaaa'

# p = Macs2(ip=tre, control=ctl, genome='dm3', output=out, prefix=None)
# # p.callpeak()
# # p.bdgcmp()
# d = p.get_effect_size()
# print(d)
# # p.bdgcmp(opt='ppois')


