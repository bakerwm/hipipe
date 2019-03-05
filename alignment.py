#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapping fastq to reference genome
1. rRNA, spikein, optional
2. genome
"""

import os
import sys
import re
import json
import fnmatch
import tempfile
import shutil
import shlex
import subprocess
import logging
import pandas as pd
import pysam
from operator import is_not
from functools import partial
from utils_parser import Json_file
from arguments import args_init
from helper import *


class AlignIndex(object):
    """1. Pick the index for aligner
    2. validate index
    """
    def __init__(self, aligner, index=None, genome=None, 
        rRNA=False, repeat_masked_genome=False, **kwargs):
        """input index, or genome
        index path to the index
        genome name of the genome
        priority: index > genome > rRNA > repeat_mask
        """
        args = args_init(kwargs, align=True)
        args['aligner'] = aligner
        args['index'] = index
        args['genome'] = genome
        args['rRNA'] = rRNA
        self.args = args

        ## determine priority
        ## fetch the index_name
        if isinstance(index, str):
            tag = self.index_validator(index, aligner)
        elif index is None:
            if isinstance(genome, str):
                tag = self.index_finder(aligner, rRNA, repeat_masked_genome)
            elif genome is None:
                raise Exception('either index= or genome=, is required')
            else:
                raise Exception('unknown argument, genome=, expect NoneType or str')
        else:
            raise Exception('unknown argument, index=, expect NoneType or str')
        self.tag = tag


    def get_index(self):
        return self.tag


    def get_index_name(self):
        """Return the name of the index
        dm3, rRNA, ...
        priority:
        index > rRNA > genome
        
        1. bowtie/bowtie2
        /path-to-data/genome/dm3/bowtie_index/genome
        /path-to-data/genome/dm3/bowtie_index/MT_trRNA
        /path-to-data/genome/dm3/dm3_transposon/bowtie_index/dm3_transposon

        2. STAR
        /path-to-data/genome/dm3/STAR_index/genome/
        /path-to-data/genome/dm3/dm3_transposon/STAR_index/
        """
        args = self.args.copy()
        index_base = os.path.basename(args['index'])

        if index_base == 'STAR_index':
            prefix = os.path.basename(os.path.dirname(args['index']))
        elif index_base == 'genome':
            prefix = os.path.basename(os.path.dirname(os.path.dirname(args['index'])))
        else:
            prefix = index_base
        return prefix


    def index_validator(self, index, aligner):
        """Validate the index for aligner
        search the aligner in $PATH
        1. bowtie: bowtie-inspect -s <index>
        2. bowtie2: bowtie2-inspect -s <index>
        3. STAR: <index>/Genome, file exists
        4. ...
        """
        args = self.args.copy()

        # check command
        aligner_exe = which(aligner)
        if not aligner_exe:
            raise Exception('aligner not detected in $PATH: %s' % aligner)

        flag = None
        if index is None:
            flag = None
        if aligner.lower().startswith('bowtie'):
            # bowtie, bowtie2
            c = [aligner_exe + '-inspect', '-s', index]
            p = subprocess.run(c, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if len(p.stdout) > 0:
                flag = index
        elif aligner.lower() == 'star':
            p = os.path.join(index, 'Genome')
            if os.path.exists(p):
                flag = index
        else:
            raise ValueError('unknown aligner: %s' % aligner)

        return flag


    def index_finder(self, aligner, rRNA=False, repeat_masked_genome=False):
        """Find the index for aligner: STAR, bowtie, bowtie2
        if rRNA is True, return the rRNA index only
        if return None, the index not found, or file not exists
        check if repeat masked genome required (for non-TE mapping)
        
        structure of genome_path:
        default: {HOME}/data/genome/{genome_version}/{aligner}/

        /genome_path/
            |- genome
            |- rRNA
            |- MT_trRNA
            |- 

        # priority
        rRNA > repeat_masked >

        """
        args = self.args.copy()

        # rRNA
        if rRNA:
            # choose the first one
            rRNA_list = ['MT_trRNA', 'rRNA']
            index_rRNA = [os.path.join(args['genome_path'], args['genome'], aligner + '_index', i) for i in rRNA_list]
            index_rRNA = [i for i in index_rRNA if self.index_validator(i, aligner)]
            if len(index_rRNA) >= 1:
                index = index_rRNA[0]
            else:
                index = None
        # genome
        else:
            if repeat_masked_genome:
                tag = 'genome_rm'
            else:
                tag = 'genome'
            index = os.path.join(args['genome_path'], args['genome'], aligner + '_index', tag)

        # validate
        if not self.index_validator(index, aligner):
            index = None
        return index


class AlignNode(object):
    """Run alignment for SE reads using bowtie/bowtie2/STAR, BWA, HISAT, kallisto, ...
    support SE mode
    (to-do: PE mode)
    """
    def __init__(self, fq1, path_out, aligner, index, smp_name=None, **kwargs):
        """Parse arguments
        fq1 fastq file or a list of files
        required arguments: fqs, path_out, smp_name, genome, genome, spikein,
        index_ext, threads, unique_only, n_map, aligner, align_to_rRNA,
        """
        ## default arguments
        args = args_init(kwargs, trim=False, align=True, call_peak=False)

        ## check index 
        tmp = args.pop('aligner')
        assert AlignIndex(aligner=aligner, index=index, **args).get_index()

        ## global variables
        self.fq1 = fq1
        self.path_out = path_out
        self.aligner = aligner
        self.index = index
        self.smp_name = smp_name # names for merge name of replicates

        ## return dict
        self.kwargs = args


    def align_init(self, fq_in, path_out):
        """Create directories for each fastq file, 
        args: fq, the fastq file
        args: index, the aligner index file
        args: align_path, output directory of the alignment    

        return files:
        prefix, bam, bed, log, unmap
        """
        args = self.kwargs.copy()

        ## prefix
        fq_prefix = file_prefix(fq_in)[0]
        fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
        if(len(self.fq1) == 1 and not self.smp_name is None):
            fq_prefix = self.smp_name

        fq_type = seq_type(fq_in)

        ## sub directory
        index_name = AlignIndex(aligner=self.aligner, index=self.index, **args).get_index_name()
        align_path = os.path.join(path_out, fq_prefix + '.map_to_' + index_name)
        assert is_path(align_path)

        ## output files
        map_prefix = os.path.join(align_path, '%s.map_%s' % (fq_prefix, index_name))
        unmap_prefix = os.path.join(align_path, '%s.not_%s' % (fq_prefix, index_name))
        map_bam = map_prefix + '.bam'
        map_log = map_prefix + '.log'
        unmap_file = unmap_prefix + '.' + fq_type
        return [fq_prefix, map_bam, map_log, unmap_file]


    def wrap_log(self, log):
        """Wrapper alignment log file, save as json"""
        args = self.kwargs.copy()
        j_file = Alignment_log(log, args['unique_only']).saveas() # save as json


    def bowtie_se(self, fq1_in):
        """Run bowtie for single file
        args: fq, the fastq file
        args: index, the path to the aligner index
        args: reference, the name of the index, keys of index_init
        args: unique_map, booolen, 
        args: align_path, None or str

        arguments:
        reference: genome, genome_rRNA, spikein, spikein_rRNA
        """
        args = self.kwargs.copy()
        bowtie_exe = which('bowtie')

        # output
        fq_prefix, map_bam, map_log, unmap_fq = self.align_init(fq1_in,)

        # determine parameters
        n_map = args['n_map']
        if n_map < 1:
            n_map = 1 # default
        if unique_map:
            para_unique = '-m 1'
        else:
            para_unique = '-v 2 -k %s' % n_map # default: 1

        if seq_type(fq) == 'fasta':
            para_fq = '-f'
        else:
            para_fq = '-q'

        # run
        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('bam file exists: %s' % map_bam)
        else:
            c1 = '%s %s %s -p %s --mm --best --sam --no-unal --un %s %s \
                %s' % (bowtie_exe, para_fq, para_unique, args['threads'], 
                    unmap_fq, index, fq1_in)
            c2 = 'samtools view -bhS -F 0x4 -@ %s -' % args['threads']
            c3 = 'samtools sort -@ %s -o %s -' % (args['threads'], map_bam)
            with open(map_log, 'wt') as ff:
                p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE,
                                      stderr=ff)
                p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout,
                                      stdout=subprocess.PIPE)
                p3 = subprocess.Popen(shlex.split(c3), stdin=p2.stdout)
                px = p3.communicate()

        # process log file
        self.wrap_log(map_log)
        return [map_bam, unmap_fq]


    def bowtie2_se(self, fq1_in):
        """Run bowtie2 for SE reads
        args: fq, the fastq file
        args: index, the path to the aligner index
        args: reference, the name of the index, keys of index_init
        args: unique_map, booolen, 
        args: align_path, None or str

        # unique mapping
        samtools view -q 10 to extract uniquely mapped reads

        arguments:
        reference: genome, genome_rRNA, spikein, spikein_rRNA
        """
        args = self.kwargs.copy()
        bowtie2_exe = which('bowtie2')

        # output directory
        fq_prefix, map_bam, map_log, unmap_fq = self.align_init(fq1_in)

        # determine parameters
        if unique_map:
            para_unique = '-q 10'
        else:
            para_unique = '-q 0'
        
        # multi map
        n_map = args['n_map']
        if n_map == 0:
            # n_map = 1 # default 1, report 1 hit for each read
            # default: #look for multiple alignments, report best, with MAPQ
            para_fq = ''
        else:
            para_fq = '-k %s' % n_map

        # fq type
        if seq_type(fq) == 'fasta':
            para_fq = para_fq + ' -f'
        else:
            para_fq = para_fq + ' -q'

        # file exists
        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('bam file exists: %s' % map_bam)
        else:
            c1 = '%s %s -p %s --very-sensitive-local --mm --no-unal --un %s -x %s -U %s' % (bowtie2_exe, 
                para_fq, args['threads'], unmap_fq, index, fq1_in)
            c2 = 'samtools view -bhS -F 0x4 -@ %s %s -' % (args['threads'], para_unique)
            c3 = 'samtools sort -@ %s -o %s -' % (args['threads'], map_bam)
            with open(map_log, 'wt') as ff:
                p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE,
                                      stderr=ff)
                p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout,
                                      stdout=subprocess.PIPE)
                p3 = subprocess.Popen(shlex.split(c3), stdin=p2.stdout)
                px = p3.communicate()

        # process log file
        self.wrap_log(map_log)
        return [map_bam, unmap_fq]    


    def star_se(self, fq1_in, path_out):
        """Run STAR for SE reads
        args: fq, the fastq file
        args: index, the path to the aligner index
        args: reference, the name of the index, keys of index_init
        args: unique_map, booolen, 
        args: align_path, None or str
        
        use '--outFilterMultimapNmax' to control uniquely mapped reads
        """
        args = self.kwargs.copy()
        star_exe = which('STAR')

        # output directory
        fq_prefix, map_bam, map_log, unmap_fq = self.align_init(fq1_in, path_out)

        # determine parameters
        n_map = args['n_map']
        if n_map > 1:
            n_map = n_map # n_map default: 0
        else:
            n_map = 10 # STAR default: 10
        para_unique = '--outFilterMultimapNmax %s' % n_map

        file_reader = 'zcat' if is_gz(fq1_in) else '-'
        # file exists
        map_prefix = os.path.splitext(map_bam)[0]
        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('bam file exists: %s' % map_bam)
        else:
            c1 = 'STAR --runMode alignReads \
              --genomeDir %s \
              --readFilesIn %s \
              --readFilesCommand %s \
              --outFileNamePrefix %s \
              --runThreadN %s \
              --limitOutSAMoneReadBytes 1000000 \
              --genomeLoad NoSharedMemory  \
              --limitBAMsortRAM 10000000000 \
              --outSAMtype BAM SortedByCoordinate \
              --outFilterMismatchNoverLmax 0.07 \
              --seedSearchStartLmax 20 \
              --outReadsUnmapped Fastx %s %s' % (self.index, fq1_in, file_reader, map_prefix, 
                args['threads'], unmap_fq, para_unique)
            p1 = subprocess.run(shlex.split(c1))
            
            # filter unique mapped reads
            if args['unique_only']: # only unique mapped reads, -q 10
                pysam.view('-bhS', '-q', '10', '-@', str(args['threads']),
                    '-o', map_bam, map_prefix + 'Aligned.sortedByCoord.out.bam',
                    catch_stdout=False)
            else:
                os.rename(map_prefix + 'Aligned.sortedByCoord.out.bam', map_bam)
            os.rename(map_prefix + 'Unmapped.out.mate1', unmap_fq)
            # os.rename(map_prefix + 'Log.final.out', map_log)
            shutil.copyfile(map_prefix + 'Log.final.out', map_log)

        # process log file
        self.wrap_log(map_log)
        return [map_bam, unmap_fq]


    def run(self):
        """Map multiple fastq files to one index
        output directory
        """
        args = self.kwargs.copy()

        # define aligner
        aligner_dict = {
            'bowtie': self.bowtie_se,
            'bowtie2': self.bowtie2_se,
            'STAR': self.star_se}

        aligner_exe = aligner_dict.get(self.aligner, None)
        if not aligner_exe:
            raise Exception('unknown aligner: %s' % self.aligner)

        ## for replicates
        bam_files = []
        unmap_files = []
        for fq in self.fq1:
            ## prefix
            fq_prefix = file_prefix(fq)[0]
            fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
            fq_path_out = os.path.join(self.path_out, fq_prefix)
            fq_bam, fq_unfq = aligner_dict[self.aligner](fq, fq_path_out)
            bam_files.append(fq_bam)
            unmap_files.append(fq_unfq)

        ## merge replicates
        if args['merge_rep'] and len(bam_files) > 1:
            if self.smp_name is None:
                merge_prefix = str_common(os.path.basename(bam_files), prefix=True)
            else:
                merge_prefix = self.smp_name
            ## merge
            index_name = AlignIndex(aligner=self.aligner, index=self.index, **args).get_index_name()
            merge_path = os.path.join(self.path_out, merge_prefix, index_name)
            merge_bam = os.path.join(merge_path, merge_prefix + '.bam')
            assert is_path(merge_path)
            if os.path.exists(merge_bam) and args['overwrite'] is False:
                logging.info('file exists: %s' % merge_bam)
            else:
                tmp = bam_merge(bam_files, merge_bam)
            bam_files.append(merge_bam)

        return [bam_files, unmap_files]


class AlignHub(object):
    """Alignment fastq files to multiple index
    multiple fastq
    multiple index (sequential)
    """

    def __init__(self, fq1, path_out, aligner, smp_name=None, 
        genome=None, **kwargs):
        """Mulitple index"""
        args = args_init(kwargs, align=True)
        #args['fq1'] = args.get('fq1', fq1)
        args['fq1'] = fq1
        args['path_out'] = path_out
        args['aligner'] = aligner
        args['smp_name'] = smp_name
        args['genome'] = genome
        self.args = args


    def index_init(self):
        """Create index, order
        default:
        1. genome_rRNA, genome, spikein_rrNA, spikein
        2. extra_index1
        3. extra_index2
        ...
        """

        args = self.args.copy()

        print(args['aligner'])
        ## nested list
        index_genome = None
        index_rRNA = None
        index_sp = None
        index_sp_rRNA = None

        aligner = args.pop('aligner')
        genome = args.pop('genome')
        rRNA = args.pop('align_to_rRNA')
        spikein = args.pop('spikein')
        
        # group1 - genome
        if isinstance(genome, str):
            index_genome = AlignIndex(aligner=aligner, genome=genome, **args).get_index()
            if rRNA: # align to rRNA
                index_rRNA = AlignIndex(aligner=aligner, genome=genome,
                    rRNA=rRNA, **args).get_index()
        # group1 - spikein
        if isinstance(spikein, str):
            index_sp = AlignIndex(aligner=aligner, genome=spikein, **args).get_index()
            if rRNA: # align to rRNA
                index_sp_rRNA = AlignIndex(aligner=aligner, 
                    genome=spikein, rRNA=rRNA, **args).get_index()
        index_group1 = [index_genome, index_rRNA, index_sp, index_sp_rRNA]

        # group2 - extra
        if args['extra_index'] is None:
            index_group2 = []
        else:
            index_group2 = [AlignIndex(aligner=aligner, index=i, **args).get_index() for i in args['extra_index']]

        # check
        n1 = sum([x is not None for x in index_group1])
        n2 = sum([x is not None for x in index_group2])

        # return
        if index_genome is None and n2 == 0:
            raise Exception('align index error, check: -g, -k, -x')

        return [index_group1, index_group2]

        
    def align_se_batch(self):
        """Align multiple index,
        priority:
        extra_index > genome
        """
        args = self.args.copy()
        index_list = self.index_init() # genome, extra

        ## arguments
        fq1 = args.pop('fq1')
        path_out = args.pop('path_out')
        aligner = args.pop('aligner')
        smp_name = args.pop('smp_name')

        n1 = sum([x is not None for x in index_list[0]]) # genome
        n2 = sum([x is not None for x in index_list[1]]) # extra

        ## extra
        map_files = []
        if n2 > 0:
            fq_x_in = fq1
            for k in index_list[1]:
                if k is None:
                    continue
                bam_x_files, unmap_x_files = AlignNode(
                    fq1=fq_x_in,
                    path_out=path_out,
                    aligner=aligner,
                    index=k,
                    smp_name=smp_name,
                    **args).run()
                fq_x_in = unmap_x_files
                map_files.append(bam_x_files)
            ## check-point
            if len(map_files) == 0:
                raise Exception('align to extra_index failed, -x:')
        else:
        ## genome
            fq_in = fq1
            for i in index_list[0]:
                if i is None:
                    continue
                bam_files, unmap_files = AlignNode(
                    fq1=fq_in,
                    path_out=path_out,
                    aligner=aligner,
                    index=i,
                    smp_name=smp_name,
                    **args).run()
                fq_in = unmap_files
                ## return
                if index_list[0].index(i) == 0:
                    map_files = bam_files # map genome

        return map_files


    def get_bam_files(self):
        """Return the bam files
        genome_rRNA, genome, spikein_rRNA, spikein, extra_index
        """
        args = self.args.copy()
        index_list = self.index_init()

        n1 = sum([x is not None for x in index_list[0]]) # genome
        n2 = sum([x is not None for x in index_list[1]]) # extra

        ## extra
        map_files = []
        return True


    def run(self):
        """run alignment"""
        map_files = self.align_se_batch()


class Alignment(object):
    """Run alignment"""

    def __init__(self, fq1, path_out, aligner, smp_name=None, 
        genome=None, **kwargs):
        """Mulitple index"""
        args1 = args_init(kwargs, align=True)
        args2 = {
            'fq1': fq1,
            'path_out': path_out,
            'aligner': aligner,
            'smp_name': smp_name,
            'genome': genome}
        self.args = {**args1, **args2}

    def run(self):
        # print(self.args)
        tmp = AlignHub(**self.args).run()


class Alignment_log(object):
    """Wrapper log file of aligner, bowtie, bowtie2, STAR
    report: total reads, unique mapped reads, multiple mapped reads

    Bowtie2:

    10000 reads; of these:
      10000 (100.00%) were unpaired; of these:
        166 (1.66%) aligned 0 times
        2815 (28.15%) aligned exactly 1 time
        7019 (70.19%) aligned >1 times
    98.34% overall alignment rate


    Bowtie:

    # reads processed: 10000
    # reads with at least one reported alignment: 3332 (33.32%)
    # reads that failed to align: 457 (4.57%)
    # reads with alignments suppressed due to -m: 6211 (62.11%)

    or:

    # reads processed: 10000
    # reads with at least one reported alignment: 9543 (95.43%)
    # reads that failed to align: 457 (4.57%)


    STAR:
    *final.Log.out

                                 Started job on |       Sep 12 11:08:57
                             Started mapping on |       Sep 12 11:11:27
                                    Finished on |       Sep 12 11:11:29
       Mapping speed, Million of reads per hour |       18.00

                          Number of input reads |       10000
                      Average input read length |       73
                                    UNIQUE READS:
                   Uniquely mapped reads number |       47
                        Uniquely mapped reads % |       0.47%
                          Average mapped length |       51.66
                       Number of splices: Total |       5
            Number of splices: Annotated (sjdb) |       0
                       Number of splices: GT/AG |       3
                       Number of splices: GC/AG |       0
                       Number of splices: AT/AC |       0
               Number of splices: Non-canonical |       2
                      Mismatch rate per base, % |       2.14%
                         Deletion rate per base |       0.04%
                        Deletion average length |       1.00
                        Insertion rate per base |       0.00%
                       Insertion average length |       0.00
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       83
             % of reads mapped to multiple loci |       0.83%
        Number of reads mapped to too many loci |       19
             % of reads mapped to too many loci |       0.19%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.02%
                 % of reads unmapped: too short |       98.31%
                     % of reads unmapped: other |       0.18%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
    """

    def __init__(self, log, unique_only=False):
        self.log = log
        self.unique_only = unique_only
        # stat
        if isinstance(log, Alignment_log):
            self.stat = log.stat
        elif isinstance(log, dict):
            self.stat = log
        elif os.path.isfile(log):
            self.stat = self._log_parser()
        else:
            raise ValueError('not supported file: %s' % log)


    def guess_aligner(self):
        """Guess the aligner of the log file:
        bowtie, bowtie2, STAR, ...
        """
        # read through log file
        log_lines = []
        with open(self.log, 'rt') as ff:
            for r in ff:
                if r.startswith('Warning'):
                    continue
                log_lines.append(r.strip())

        # parsing log file
        line = log_lines[0] # the first line
        if line.startswith('#'):
            log_parser = self._bowtie_parser
        elif 'reads; of these' in line:
            log_parser = self._bowtie2_parser
        elif '|' in line:
            log_parser = self._star_parser
        else:
            raise ValueError('unknown file format: %s' % self.log)
            pass
        return log_parser


    def _is_non_empty(self):
        """Check if log file is empty"""
        if os.path.getsize(self.log) > 0:
            return True
        else:
            return False


    def _bowtie_parser(self):
        """Wrapper bowtie log
        unique, multiple, unmap, map, total
        """
        dd = {}
        with open(self.log, 'rt') as ff:
            for line in ff:
                if not ':' in line or line.startswith('Warning'):
                    continue
                num = line.strip().split(':')[1]
                value = num.strip().split(' ')[0]
                value = int(value)
                if 'reads processed' in line:
                    dd['total'] = value
                elif 'at least one reported alignment' in line:
                    dd['map'] = value
                elif 'failed to align' in line:
                    dd['unmap'] = value
                elif 'alignments suppressed due to -m' in line:
                    dd['multiple'] = value
                else:
                    pass
        # unique_only
        dd['unique'] = dd['map']
        dd['multiple'] = dd.get('multiple', 0) # default 0
        if self.unique_only:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']
        return dd


    def _bowtie2_parser(self):
        """Wrapper bowtie2 log"""
        dd = {}
        with open(self.log, 'rt') as ff:
            for line in ff:
                value = line.strip().split(' ')[0]
                if '%' in value:
                    continue
                value = int(value)
                if 'reads; of these' in line:
                    dd['total'] = value
                elif 'aligned 0 times' in line:
                    dd['unmap'] = value
                elif 'aligned exactly 1 time' in line:
                    dd['unique'] = value
                elif 'aligned >1 times' in line:
                    dd['multiple'] = value
                else:
                    pass
        if self.unique_only:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']
        return dd


    def _star_parser(self):
        """Wrapper STAR *final.log"""
        dd = {}
        with open(self.log, 'rt') as ff:
            for line in ff:
                value = line.strip().split('|')
                if not len(value) == 2:
                    continue
                value = value[1].strip()
                if 'Number of input reads' in line:
                    dd['total'] = int(value)
                elif 'Uniquely mapped reads number' in line:
                    dd['unique'] = int(value)
                elif 'Number of reads mapped to multiple loci' in line:
                    dd['multiple'] = int(value)
                else:
                    pass
        if self.unique_only is True:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']
        return dd


    def _log_parser(self):
        """Read log file as dict
        delimiter:
        bowtie:  ":"
        bowtie2:  ":"
        STAR:  "|"

        extra trimming
        1. trim "(10.00%)" 
        2. trim "blank" at both ends
        """
        log = self.log
        log_parser = self.guess_aligner()
        dd = log_parser()
        return dd


    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.json',
                                            delete=False)
        return tmpfn.name


    def saveas(self, _out=None):
        """Make a copy of statistics of mapping results"""
        log = self.log
        if _out is None:
            # _out = self._tmp()
            _out = os.path.splitext(log)[0] + '.json'

        dd = self.stat

        with open(_out, 'wt') as fo:
            json.dump(dd, fo, indent=4, sort_keys=False)

        return _out


class Alignment_stat(object):
    """Parse mapping reads in directory
    1. for each rep bam, parse json files, 
    2. merge replicates
    """
    def __init__(self, path):
        self.path = path

        if isinstance(path, Alignment_stat):
            self.stat = path.stat
        elif isinstance(path, pd.DataFrame):
            self.stat = path
        elif os.path.isdir(path):
            json_files = self.json_files()
            bam_files = self.bam_files()
            if json_files is None: # no json files
                self.stat = self.merge_stat() 
            elif len(json_files) == 1:
                self.stat = self.single_stat()
            elif json_files:
                self.stat = self.rep_stat()
            elif bam_files:
                self.stat = self.merge_stat()
            else:
                raise ValueError('BAM or json files not found: %s' % path)
        else:
            raise ValueError('unknown format')


    def _is_non_empty(self, fn):
        """Check if log file is empty"""
        if os.path.getsize(fn) > 0:
            return True
        else:
            return False


    def findfiles(self, which, where='.'):
        """Returns list of filenames from `where` path matched by 'which'
        shell pattern. Matching is case-insensitive.
        # findfiles('*.ogg')
        """    
        # TODO: recursive param with walk() filtering
        rule = re.compile(fnmatch.translate(which), re.IGNORECASE)
        hits = [os.path.join(where, name) for name in os.listdir(where) if rule.match(name)]
        return hits


    def json_files(self):
        """Return all json files in path"""
        j = self.findfiles('*.json', self.path)
        j = [i for i in j if self._is_non_empty(i)] # exclude empty json
        if len(j) == 0:
            j = None
        return j


    def bam_files(self):
        """Return all BAM files in path
        skip symlink
        # directory of merged fastq, count bam file
        # to-do: summary replicates (except unmap)
        """
        bam_files = self.findfiles('*.bam', self.path)
        bam_files = [b for b in bam_files if self._is_non_empty(b) and not os.path.islink(b)]
        if len(bam_files) == 0:
            bam_files = None
        return bam_files


    def _json_log_reader(self, fn, to_dataframe=True):
        """Parse json file, save as DataFrame"""
        if fn is None:
            return None
        fn_name = os.path.basename(fn)
        group, aligner = fn.split('.')[-3:-1] # 
        group = re.sub('^map_', '', group)
        dd = Json_file(fn).stat # dict of count
        rpt = dd
        if to_dataframe:
            # total, unique, multiple, unmap
            df = pd.DataFrame(data=[list(dd.values())], columns=list(dd.keys()))
            df.pop('map') # drop "map" column
            df.index = [fn_name]
            rpt = df
        return rpt


    def single_stat(self):
        """Extract alignment log files in index_ext directory"""
        path = self.path.rstrip('/')
        json_files = self.json_files()
        json_files = sorted(json_files, key=len)
        prefix = os.path.basename(self.path)
        dd = self._json_log_reader(json_files[0], False)
        df = pd.DataFrame(dd, index=[prefix])
        return df


    def rep_stat(self):
        """Extract alignment log files in directory
        """
        path = self.path.rstrip('/')
        json_files = self.json_files()
        json_files = sorted(json_files, key=len)
        prefix = os.path.basename(self.path)
        # genome_rRNA, genome, sp_rRNA, sp, unmap, map, total
        if len(json_files) == 2 or len(json_files) == 4:
            pass
        else:
            raise ValueError('number of json files should be; 2|4, [%d]' % len(json_files))

        # genome
        # map rRNA, genome, spikein_rRNA, spikein
        json_genome_rRNA, json_genome = json_files[:2]
        json_sp_rRNA = json_sp = None
        if len(json_files) == 4:
            json_sp_rRNA, json_sp = json_files[2:]
        df1 = self._json_log_reader(json_genome_rRNA, False)
        df2 = self._json_log_reader(json_genome, False)
        df3 = self._json_log_reader(json_sp_rRNA, False)
        df4 = self._json_log_reader(json_sp, False)

        # genome
        n_total = df1['total']
        g_rRNA = df1['unique'] + df1['multiple']
        g_unique = df2['unique']
        g_multi = df2['multiple']
        n_unmap = df2['unmap']

        # spikein rRNA
        if isinstance(df3, dict):
            sp_rRNA = df3['unique'] + df3['multiple']
        else:
            sp_rRNA = 0

        # spikein
        if isinstance(df4, dict):
            sp_unique = df4['unique']
            sp_multi = df4['multiple']
            n_unmap = df4['unmap']
        else:
            sp_unique = sp_multi = 0

        # output
        data = {
            'genome_rRNA': g_rRNA,
            'genome_unique': g_unique,
            'genome_multiple': g_multi,
            'spikein_rRNA': sp_rRNA,
            'spikein_unique': sp_unique,
            'spikein_multi': sp_multi,
            'unmap': n_unmap,
            'total': n_total}
        df = pd.DataFrame(data, index=[prefix])

        # return
        return df


    def merge_stat(self):
        """Stat reads for merged sample
        combine all reads in each replicates
        no json files detedted in merged directory
        search *.csv files in up-level directory
        """
        merge_path_name = os.path.basename(self.path.rstrip('/'))
        parent_path = os.path.dirname(self.path.rstrip('/')) # 
        rep_path = [i for i in os.listdir(parent_path) if not i == merge_path_name]
        rep_csv_files = [os.path.join(parent_path, i + '.mapping_stat.csv') for i in rep_path]
        rep_csv_files = [f for f in rep_csv_files if os.path.isfile(f)]

        if len(rep_csv_files) > 0:
            frames = [pd.read_csv(i, '\t', index_col=0) for i in rep_csv_files]
            df = pd.concat(frames, axis=0)
            # merge
            df_merge = pd.DataFrame(data=[df.sum(axis=0)], columns=list(df.columns.values),
                index=[merge_path_name])
            return df_merge
        else:
            logging.error('%10s | not contain merged mapping files: %s' % ('failed', self.path))
            return None


    def saveas(self, _out=None):
        """Make a copy of statistics of mapping results"""
        path = self.path
        if _out is None:
            prefix = os.path.basename(path)
            _out = os.path.join(os.path.dirname(path),
                                prefix + '.mapping_stat.csv')
        df = self.stat.reset_index().rename(columns={'index': 'name'})
        if isinstance(df, pd.DataFrame):
            df.to_csv(_out, sep='\t', header=True, index=False)
        
        return _out

