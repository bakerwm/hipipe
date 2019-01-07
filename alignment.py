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


class Alignment(object):
    """Run alignment for SE reads using bowtie/bowtie2/STAR, BWA, HISAT, kallisto, ...
    support SE mode
    (to-do: PE mode)
    """
    def __init__(self, fq1, path_out, smp_name, genome, **kwargs):
        """Parse arguments
        required arguments: fqs, path_out, smp_name, genome, genome, spikein,
        index_ext, threads, unique_only, n_map, aligner, align_to_rRNA,
        genome_path, overwrite
        """
        args1 = args_init(kwargs, trim=False, align=True, call_peak=False)
        args2 = {
            'fq1': fq1, 
            'path_out': path_out,
            'smp_name': smp_name,
            'genome': genome}
        args = {**args1, **args2} # update arguments
        # check
        if args['spikein'] == args['genome']:
            args['spikein'] = None
        self.kwargs = args


    def align_init(self, fq, index, reference=None, align_path=None):
        """Create directories for the alignment, 
        args: fq, the fastq file
        args: index, the aligner index file
        args: reference, the name of the index
        args: align_path, output directory of the alignment    

        Alignment, genome versions
        1.genome_rRNA
        2.genome
        3.spikein_rRNA
        4.spikein

        return files:
        prefix, bam, bed, log, unmap
        """
        args = self.kwargs.copy()

        if not reference:
            reference = args['genome'] # default is reference genome
        if not align_path:
            align_path = args['path_out']

        fq_prefix = file_prefix(fq)[0]
        fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
        fq_type = seq_type(fq)

        map_prefix = os.path.join(align_path, '%s.map_%s' % (fq_prefix, reference))
        unmap_prefix = os.path.join(align_path, '%s.not_%s' % (fq_prefix, reference))
        map_bam = map_prefix + '.bam'
        unmap_fq = unmap_prefix + '.%s' % fq_type
        map_log = map_prefix + '.%s.log' % args['aligner']

        return [fq_prefix, map_bam, map_log, unmap_fq, reference]

    
    def index_builder(self, rRNA=False, genome=None):
        """Return the genome index
        args: rRNA, True/False, whether return the rRNA index
        args: genome, None or str, the name of genome, 
        return the index
        """
        args = self.kwargs.copy()
        aligner = args['aligner']

        if not genome:
            genome = args['genome']
        # check aligner
        if aligner == 'bowtie':
            index = Genome(genome, repeat_masked_genome=args['repeat_masked_genome']).bowtie_index(rRNA=rRNA)
        elif aligner == 'bowtie2':
            index = Genome(genome, repeat_masked_genome=args['repeat_masked_genome']).bowtie2_index(rRNA=rRNA)
        elif aligner == 'STAR':
            index = Genome(genome, repeat_masked_genome=args['repeat_masked_genome']).star_index(rRNA=rRNA)
        else:
            logging.error('unknown aligner: %s' % aligner)
            index = None # unknonwn aligner
        return index


    def index_init(self):
        """Prepare the indexes for alignment
        default: genome_rRNA, genome, spikein_rRNA, spikein
        # args required
        
        """
        args = self.kwargs.copy()

        # aligner index
        idx = {
            'genome_rRNA': self.index_builder(rRNA=True),
            'genome': self.index_builder(rRNA=False),
            'sp_rRNA': self.index_builder(rRNA=True, genome=args['spikein']),
            'sp': self.index_builder(rRNA=False, genome=args['spikein'])}

        # validate
        if not args['align_to_rRNA']:
            idx['genome_rRNA'] = None

        if idx['genome'] is None:
            raise Exception('aligner index not found for : %s' % args['genome'])

        # save in dict
        return idx # dictionary


    def wrap_log(self, log):
        """Wrapper alignment log file, save as json"""
        args = self.kwargs.copy()
        j_file = Alignment_log(log, args['unique_only']).saveas() # save as json


    def bowtie_se(self, fq, index, reference=None, unique_map=False, align_path=None):
        """Run bowtie for SE reads
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
        if align_path is None:
            align_path = args['path_out']

        # output directory
        fq_prefix, map_bam, map_log, unmap_fq, reference = self.align_init(fq, index, 
            reference, align_path)

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
                    unmap_fq, index, fq)
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


    def bowtie2_se(self, fq, index, reference=None, unique_map=False, align_path=None):
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
        if not align_path:
            align_path = args['path_out']

        # output directory
        prefix, map_bam, map_log, unmap_fq, reference = self.align_init(fq, 
            index, reference, align_path)
        
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
                para_fq, args['threads'], unmap_fq, index, fq)
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


    def star_se(self, fq, index, reference, unique_map=False, align_path=None):
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
        if not align_path:
            align_path = args['path_out']

        # output directory
        prefix, map_bam, map_log, unmap_fq, reference = self.align_init(fq, 
            index, reference, align_path)
        
        # determine parameters
        n_map = args['n_map']
        if n_map > 1:
            n_map = n_map # n_map default: 0
        else:
            n_map = 10 # STAR default: 10
        para_unique = '--outFilterMultimapNmax %s' % n_map

        fr = 'zcat' if is_gz(fq) else '-'
        # file exists
        map_prefix = os.path.join(align_path, prefix)
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
              --outReadsUnmapped Fastx %s %s' % (index, fq, fr, map_prefix, 
                args['threads'], unmap_fq, para_unique)
            p1 = subprocess.run(shlex.split(c1))
            
            # filter unique mapped reads
            if unique_map: # only unique mapped reads, -q 10
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


    def align_batch_se(self, fq, align_path=None):
        """Align reads to multiple indexes in specific order,
        return align.stat, map_bam, unmap_reads
        determine the index order
        return bam_files
        """
        args = self.kwargs.copy()

        if not align_path:
            align_path = args['path_out']

        # define aligner
        aligner_dict = {
            'bowtie': self.bowtie_se,
            'bowtie2': self.bowtie2_se,
            'STAR': self.star_se}

        aligner_exe = aligner_dict.get(args['aligner'], None) # determine aligner

        if not aligner_exe:
            raise ValueError('unknown aligner: %s' % args['aligner'])

        # get all index in order
        index_dict = self.index_init() # genome_rRNA, genome, sp_rRNA, sp
        bam_files = []
        fq_input = fq

        # 1. genome_rRNA (rRNA: both unique, multiple)
        idx1 = index_dict['genome_rRNA']
        if idx1 is None:
            raise ValueError('genome_rRNA index not found: %s' % args['genome'])
        reference = args['genome'] + '_rRNA'
        bam_idx1, unmap_idx1 = aligner_exe(fq=fq_input, index=idx1, 
            reference=reference, unique_map=False, 
            align_path=align_path)
        fq_input = unmap_idx1

        # 2. genome
        idx2 = index_dict['genome']
        if idx2 is None:
            raise ValueError('genome index not found: %s' % args['genome'])

        reference = args['genome']
        bam_idx2, unmap_idx2 = aligner_exe(fq=fq_input, index=idx2, 
            reference=reference, unique_map=args['unique_only'], 
            align_path=align_path)
        fq_input = unmap_idx2

        if args['spikein']: # add spikein
            # 3. sp_rRNA  (rRNA: both unique, multiple)
            idx3 = index_dict['sp_rRNA']
            reference = args['spikein'] + '_rRNA'
            bam_idx3, unmap_idx3 = aligner_exe(fq=fq_input, index=idx3, 
                reference=reference, unique_map=False,
                align_path=align_path)
            fq_input = unmap_idx3

            # 4. sp (optional)
            idx4 = index_dict['sp']
            reference = args['spikein']
            bam_idx4, unmap_idx4 = aligner_exe(fq=fq_input, index=idx4, 
                reference=reference, unique_map=args['unique_only'],
                align_path=align_path)
            fq_input = unmap_idx4

            bam_files = [bam_idx1, bam_idx2, bam_idx3, bam_idx4]
        else: # no spikein
            bam_idx3 = bam_idx4 = None

        bam_files = [bam_idx1, bam_idx2, bam_idx3, bam_idx4]
        bam_files = list(filter(partial(is_not, None), bam_files)) # remove None

        return bam_files


    def align_extra(self, fq, align_path=None):
        """Align reads to extra index
        such as GFP, white, firefly, transposon
        return bam_files
        """
        args = self.kwargs.copy()

        fq_prefix = file_prefix(fq)[0]
        fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)

        if align_path is None:
            align_path = args['path_out']

        # define aligner
        aligner_dict = {
            'bowtie': self.bowtie_se,
            'bowtie2': self.bowtie2_se,
            'STAR': self.star_se}

        aligner_exe = aligner_dict.get(args['aligner'], None) # determine aligner
        if not aligner_exe:
            raise ValueError('unknown aligner: %s' % args['aligner'])

        # get all index in order
        # index_ext = args['index_ext']
        bam_ext_list = []

        for ext in args['index_ext']:
            if index_validator(ext, args['aligner']):
                print(ext)
                ext_order = args['index_ext'].index(ext) # number
                fq_path = os.path.join(align_path, 'extra_mapping_' + str(ext_order), fq_prefix)
                assert is_path(fq_path)
                reference = os.path.basename(ext)
                bam_ext, unmap_ext = aligner_exe(
                    fq=fq, 
                    index=ext, 
                    reference=reference, 
                    unique_map=True, 
                    align_path=fq_path)
            else:
                bam_ext = None
            bam_ext_list.append(bam_ext)
        return bam_ext_list


    def run_extra(self):
        """Run the alignment for specific fastq file onto extra index
        1. run alignment for each replicate
        2. merge replicates
        3. run log parser, in json format
        4. organize the log files, saved in one report, including the following groups:
        """
        args = self.kwargs.copy()

        bam_out = []

        # run alignment for replicates
        for fq in args['fq1']:
            # logging.info('alignment: %s' % fq)
            fq_prefix = file_prefix(fq)[0]
            fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
            # fq_path = os.path.join(args['path_out'], 'extra_mapping', fq_prefix)
            # assert is_path(fq_path)
            logging.info('align to index_ext: %s' % fq_prefix)
            bam_files = self.align_extra(fq)
            bam_out.append(bam_files) #
            # stat alignment
            Alignment_stat(fq_path).saveas()

        # merge bam files
        if args['merge_rep']: 
            merged_path = os.path.join(args['path_out'], 'extra_mapping', args['smp_name'])
            merged_files = []
            if len(bam_out) > 1: # for multiple bam files
                assert is_path(merged_path)
                for i in range(len(bam_out[0])):
                    rep_bam_files = [b[i] for b in bam_out]
                    merged_suffix = str_common(rep_bam_files, suffix=True)
                    merged_suffix = re.sub('^_[12]|_R[12]', '', merged_suffix)
                    merged_bam_name = args['smp_name'] + merged_suffix
                    merged_bam_file =  os.path.join(merged_path, merged_bam_name)
                    if os.path.exists(merged_bam_file) and args['overwrite'] is False:
                        logging.info('bam file exists: %s' % merged_bam_file)
                    else:
                        tmp = bam_merge(rep_bam_files, merged_bam_file)
                    merged_files.append(merged_bam_file)
                Alignment_stat(merged_path).saveas()
                bam_out.append(merged_files)

        return bam_out


    def run(self):
        """Run the alignment for specific fastq file
        1. run alignment for each replicate
        2. merge replicates
        3. run log parser, in json format
        4. organize the log files, saved in one report, including the following groups:
        genome_rRNA, genome_unique, genome_multi, sp_rRNA, sp_unique, sp_multi, unmap
        """
        args = self.kwargs

        bam_out = []

        # run alignment for replicates
        for fq in args['fq1']:
            logging.info('alignment: %s' % fq)
            args_fq = args.copy()
            fq_prefix = file_prefix(fq)[0]
            fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
            fq_path = os.path.join(args_fq['path_out'], fq_prefix)
            args_fq['path_out'] = fq_path
            assert is_path(fq_path)

            ## save arguments
            args_file = os.path.join(fq_path, fq_prefix + '.arguments.txt')
            args_pickle = os.path.join(fq_path, fq_prefix + '.arguments.pickle')
            if args_checker(args_fq, args_pickle) and args['overwrite'] is False:
                logging.info('files exists, arguments not changed, alignment skipped - %s' % fq_prefix)
                # return True
            else:
                args_logger(args_fq, args_file, True) # update arguments.txt

            bam_files = self.align_batch_se(fq, fq_path)
            bam_out.append(bam_files) # 
            # stat alignment
            Alignment_stat(fq_path).saveas()

        # merge bam files
        if args['merge_rep']: 
            merged_path = os.path.join(args['path_out'], args['smp_name'])
            merged_files = []
            if len(bam_out) > 1: # for multiple bam files
                assert is_path(merged_path)
                for i in range(len(bam_out[0])):
                    rep_bam_files = [b[i] for b in bam_out]
                    merged_suffix = str_common(rep_bam_files, suffix=True)
                    merged_suffix = re.sub('^_[12]|_R[12]', '', merged_suffix)
                    merged_bam_name = args['smp_name'] + merged_suffix
                    merged_bam_file =  os.path.join(merged_path, merged_bam_name)
                    if os.path.exists(merged_bam_file) and args['overwrite'] is False:
                        logging.info('file exists: %s' % merged_bam_file)
                    else:
                        tmp = bam_merge(rep_bam_files, merged_bam_file)
                    merged_files.append(merged_bam_file)
                Alignment_stat(merged_path).saveas()
                bam_out.append(merged_files)

        # make short names for genome bam files
        genome_bam_files = []
        for b in bam_out: # nested array
            bam_from = b[1]
            bam_to = os.path.join(os.path.dirname(bam_from),
                filename_shorter(bam_from))
            if not os.path.exists(bam_to):
                os.symlink(os.path.basename(bam_from), bam_to)
            if not os.path.exists(bam_to + '.bai'):
                if not os.path.exists(bam_from + '.bai'):
                    if os.path.getsize(bam_from) < 1000: # !!!! empty bam files
                        continue
                    else:
                        pysam.index(bam_from) # empty bam
                os.symlink(os.path.basename(bam_from) + '.bai', bam_to + '.bai')
            genome_bam_files.append(bam_to)

        # run extra index mapping
        ext_bam_files = None
        if not args['index_ext'] is None:
            ext_bam_files = self.run_extra()
        
        print(args)
        return [genome_bam_files, ext_bam_files]


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

