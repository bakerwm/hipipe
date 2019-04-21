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
import pickle
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
    """Run alignment
    :args fq1
    :args path_out
    :args aligner
    :args smp_name
    :args genome
    :args spikein
    :args align_to_te
    :args te_index
    :args extra_index
    :args repeat_masked_genome
    :args genome_path
    """
    def __init__(self, fq1, path_out, aligner, smp_name=None, 
        genome=None, align_by_order=True, **kwargs):
        """Mulitple index"""
        args = args_init(kwargs, align=True)

        ## update arguments
        args['fq1'] = fq1
        # args['path_out'] = path_out
        args['aligner'] = aligner
        args['smp_name'] = smp_name
        args['genome'] = genome
        args['align_by_order'] = align_by_order

        ## create index
        args['index_list'] = AlignIndexBuilder(**args).get_index()
        self.fq1 = fq1
        self.path_out = path_out # original path_out
        self.args = args


    def run(self):
        args = self.args.copy()

        ## check arguments
        assert is_path(self.path_out)
        args_file = os.path.join(self.path_out, 'arguments.txt')
        args_pickle = os.path.join(self.path_out, 'arguments.pickle')
        if args_checker(args, args_pickle) and args['overwrite'] is False:
            logging.info('arguments not changed, alignment skipped')
            args['dry_run'] = True
        else:
            args_logger(args, args_file, overwrite=True) # update arguments.txt

        map_files = []
        for fq_in in self.fq1:
            fq_prefix = file_prefix(fq_in)[0]
            fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
            if args['simple_name']:
                fq_prefix = re.sub('.not_\w+', '', fq_prefix)
                fq_prefix = re.sub('.map_\w+', '', fq_prefix)

            fq_path_out = os.path.join(self.path_out, fq_prefix)
            logging.info('fq_file: %s' % fq_prefix) # !!!! logging
            ## update arguments
            args['fq1'] = fq_in
            args['path_out'] = fq_path_out
            bam_files = AlignHub(**args).run()
            map_files.append(bam_files)
            ## summarize alignment log
            align_stat = AlignSummarize(path_out=fq_path_out).run()
            # print(fq_path_out)
            # print(align_stat)

        return map_files


class AlignIndex(object):
    """1. Pick the index for aligner
    2. validate index
    """
    def __init__(self, aligner, genome=None, index=None, 
        align_to_rRNA=False, align_to_te=False, 
        repeat_masked_genome=False, **kwargs):
        """input index, or genome
        index path to the index
        genome name of the genome
        priority: index > genome > rRNA > repeat_mask
        """
        args = args_init(kwargs, align=True)
        # update arguments
        args['aligner'] = aligner
        args['genome'] = genome
        args['index'] = index # ?
        args['align_to_rRNA'] = align_to_rRNA 
        args['align_to_te'] = align_to_te
        args['repeat_masked_genome'] = repeat_masked_genome

        self.args = args

        ## fetch the index_name
        if isinstance(index, str):
            index = index.rstrip('/') # remove '/' at the end
            tag = self.index_validator(index, aligner)
        elif index is None:
            if isinstance(genome, str):
                tag = self.index_finder()
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
            # args['index'] = args['index'].rstrip('/') # remove '/' in the end
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


    def index_finder(self):
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

        # te
        if args['align_to_te']:
            # choose TE index
            index = os.path.join(args['genome_path'], 
                args['genome'], 
                args['aligner'] + '_index', 
                args['genome'] + '_transposon')
        # rRNA
        elif args['align_to_rRNA']:
            # choose the first one
            rRNA_list = ['MT_trRNA', 'rRNA']
            index_rRNA = [os.path.join(args['genome_path'], 
                args['genome'], 
                args['aligner'] + '_index', i) for i in rRNA_list]
            index_rRNA = [i for i in index_rRNA if self.index_validator(i, args['aligner'])]
            if len(index_rRNA) >= 1:
                index = index_rRNA[0]
            else:
                index = None
        # genome
        else:
            if args['repeat_masked_genome']:
                tag = 'genome_rm'
            else:
                tag = 'genome'
            index = os.path.join(args['genome_path'], 
                args['genome'], 
                args['aligner'] + '_index', 
                tag)

        # validate
        if not self.index_validator(index, args['aligner']):
            index = None
        return index


class AlignIndexBuilder(object):
    """Create alignment indexes
    genome + spikein
    transposon
    extra
    : arguments :
    1. genomme (str), spikein (str|None), align_to_rRNA (True|False)
    2. align_to_te (True|False), te_index (str|None) (mapping)
    3. extra_index (str|None)
    """
    def __init__(self, aligner, genome=None, spikein=None, 
        align_to_rRNA=True, align_to_te=True, te_index=None,
        extra_index=None, **kwargs):
        """Mulitple index"""
        args = args_init(kwargs, align=True)

        ## update arguments
        args['aligner'] = aligner
        args['genome'] = genome
        args['spikein'] = spikein
        args['align_to_rRNA'] = align_to_rRNA
        args['align_to_te'] = align_to_te
        args['te_index'] = te_index
        args['extra_index'] = extra_index
        self.args = args


    def index_builder(self):
        """Create index, order
        default:
        1. genome_rRNA, genome, spikein_rrNA, spikein
        2. extra_index1
        3. extra_index2
        ...
        """
        args = self.args.copy()

        align_to_te = args.pop('align_to_te')

        # group1 - extra
        if args['extra_index']: # True
            index_group = [AlignIndex(index=i, **args).get_index() for i in args['extra_index']]
        # group2 - te
        elif align_to_te:
            if args['te_index']:
                index_group = [AlignIndex(index=args['te_index'], **args).get_index()]
            else:
                index_group = [AlignIndex(align_to_te=True, **args).get_index()]
        # group3 - genome + spikein
        elif isinstance(args['genome'], str): # genome, required
            ## nested list
            index_genome = None
            index_rRNA = None
            index_sp = None
            index_sp_rRNA = None
            ## genome
            align_to_rRNA = args.pop('align_to_rRNA')
            index_genome = AlignIndex(align_to_rRNA=False, **args).get_index()
            if align_to_rRNA: # align to rRNA
                index_rRNA = AlignIndex(align_to_rRNA=True, **args).get_index()
            # spikein
            if isinstance(args['spikein'], str):
                args['genome'] = args['spikein'] # update
                index_sp = AlignIndex(align_to_rRNA=False, **args).get_index()
                if align_to_rRNA: # align to rRNA
                    index_sp_rRNA = AlignIndex(align_to_rRNA=True, **args).get_index()
            index_group = [index_rRNA, index_genome, index_sp_rRNA, index_sp]
        else:
            raise Exception('align index failed, check: -g, -k, -x, --align-to-te, --te-index')

        # check
        index_n = sum([x is not None for x in index_group])

        # return
        if index_n == 0:
            raise Exception('align index error, check: -g, -k, -x, --align-to-te, --te-index')

        return index_group


    def get_index(self):
        """Choose index
        priority:
        extra > te > genome + spikein
        """
        return self.index_builder()


class AlignNode(object):
    """Run alignment for SE reads using bowtie/bowtie2/STAR, BWA, HISAT, kallisto, ...
    **Only support only one read file to one index**
    support SE mode
    (to-do: PE mode)
    """
    def __init__(self, fq1, path_out, aligner, index, smp_name=None, **kwargs):
        """Parse arguments
        fq1 fastq file or a list of files
        required arguments: fqs, path_out, smp_name, genome, genome, spikein,
        index_ext, threads, unique_only, n_map, aligner, align_to_rRNA,
        """
        assert isinstance(fq1, str) # only one fastq file

        ## default arguments
        args = args_init(kwargs, align=True)
        
        ## update arguments
        args['fq1'] = fq1
        args['path_out'] = path_out
        args['aligner'] = aligner
        args['index'] = index
        args['smp_name'] = smp_name

        ## check index 
        assert AlignIndex(**args).get_index()

        ## check aligner
        aligner_dict = {
            'bowtie': self.bowtie_se,
            'bowtie2': self.bowtie2_se,
            'STAR': self.star_se}
       
        aligner_exe = aligner_dict.get(aligner, None)
        if not aligner_exe:
            raise Exception('unknown aligner: %s' % aligner)
        self.aligner_exe = aligner_exe

        ## return dict
        self.args = args


    def align_init(self):
        """Create directories for each fastq file, 
        args: fq, the fastq file
        args: index, the aligner index file
        args: align_path, output directory of the alignment    

        return files:
        prefix, bam, bed, log, unmap
        """
        args = self.args.copy()

        ## prefix
        fq_prefix = file_prefix(args['fq1'])[0]
        fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
        if(len(args['fq1']) == 1 and not args['smp_name'] is None):
            fq_prefix = args['smp_name']
        fq_type = seq_type(args['fq1'])

        ## simplify sample name
        if args['simple_name']:
            if '.not_' in fq_prefix:
                fq_prefix = re.sub('.not_\w+', '', fq_prefix)
                fq_prefix = re.sub('.map_\w+', '', fq_prefix)

        ## sub directory
        index_name = AlignIndex(**args).get_index_name()
        align_path = os.path.join(args['path_out'], fq_prefix + '.map_' + index_name)
        assert is_path(align_path)

        ## output files
        # print(fq_prefix)
        map_prefix = os.path.join(align_path, fq_prefix + '.map_' + index_name) #
        unmap_prefix = os.path.join(align_path, fq_prefix + '.not_' + index_name) #
        map_bam = map_prefix + '.bam'
        map_log = map_prefix + '.log'
        unmap_file = unmap_prefix + '.' + fq_type
        return [fq_prefix, map_bam, map_log, unmap_file]


    def wrap_log(self, log):
        """Wrapper alignment log file, save as json"""
        args = self.args.copy()
        j_file = Alignment_log(log, args['unique_only']).saveas() # save as json


    def bowtie_se(self):
        """Run bowtie for single file
        args: fq, the fastq file
        args: index, the path to the aligner index
        args: reference, the name of the index, keys of index_init
        args: unique_map, booolen, 
        args: align_path, None or str

        arguments:
        reference: genome, genome_rRNA, spikein, spikein_rRNA
        """
        args = self.args.copy()
        bowtie_exe = which('bowtie')

        # output
        fq_prefix, map_bam, map_log, unmap_fq = self.align_init()

        # determine parameters
        n_map = args['n_map']
        if n_map < 1:
            n_map = 1 # default
        if args['unique_only']:
            para_unique = '-m 1'
        else:
            para_unique = '-v 2 -k %s' % n_map # default: 1

        if seq_type(args['fq1']) == 'fasta':
            para_fq = '-f'
        else:
            para_fq = '-q'

        # run
        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('bam file exists: %s' % map_bam)
        else:
            c1 = '%s %s %s -p %s --mm --best --sam --no-unal --un %s %s \
                %s' % (bowtie_exe, para_fq, para_unique, args['threads'], 
                    unmap_fq, args['index'], args['fq1'])
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


    def bowtie2_se(self):
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
        args = self.args.copy()
        bowtie2_exe = which('bowtie2')

        # output directory
        fq_prefix, map_bam, map_log, unmap_fq = self.align_init()

        # determine parameters
        if args['unique_only']:
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
        if seq_type(args['fq1']) == 'fasta':
            para_fq = para_fq + ' -f'
        else:
            para_fq = para_fq + ' -q'

        # file exists
        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('bam file exists: %s' % map_bam)
        else:
            c1 = '%s %s -p %s --very-sensitive-local --mm --no-unal --un %s -x %s -U %s' % (bowtie2_exe, 
                para_fq, args['threads'], unmap_fq, args['index'], args['fq1'])
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


    def star_se(self):
        """Run STAR for SE reads
        args: fq, the fastq file
        args: index, the path to the aligner index
        args: reference, the name of the index, keys of index_init
        args: unique_map, booolen, 
        args: align_path, None or str
        
        use '--outFilterMultimapNmax' to control uniquely mapped reads

        # update 2019-03-21
        since 99% of the reads do not map to Blumeria, STAR takes a lot of time trying to squeeze the reads into the wrong genome. The best solution is to make a combined genome of Barley and Blumeria, which will alloy mapping simultaneously to the two genomes, with the best alignments winning.
        Another option (if you insist on mapping to Blumeria alone) is to reduce --seedPerWindowNmax to 20 or even smaller values. More discussion on it here: https://groups.google.com/d/msg/rna-star/hJL_DUtliCY/HtpiePlMBtYJ .
        see: https://github.com/alexdobin/STAR/issues/329#issuecomment-334879474

        """
        args = self.args.copy()
        star_exe = which('STAR')

        # output directory
        fq_prefix, map_bam, map_log, unmap_fq = self.align_init()

        # determine parameters
        n_map = args['n_map']
        if n_map > 1:
            n_map = n_map # n_map default: 0
        else:
            n_map = 10 # STAR default: 10
        # small genome
        seed_max = 50 # STAR default: 50
        if args['small_genome']:
            seed_max = 20 # even smaller

        para_unique = '--outFilterMultimapNmax %s --seedPerWindowNmax %s' % (n_map, seed_max)

        file_reader = 'zcat' if is_gz(args['fq1']) else '-'

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
              --outReadsUnmapped Fastx %s %s' % (args['index'], args['fq1'], file_reader, map_prefix, 
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
            shutil.copyfile(map_prefix + 'Log.final.out', map_log)

        # process log file
        self.wrap_log(map_log)
        return [map_bam, unmap_fq]


    def get_prefix(self):
        """Return the prefix of fastq files
        return: path
        """
        fq_prefix = self.align_init()[0]
        return fq_prefix


    def get_bam(self):
        """Return the name of BAM file
        return: path
        """
        map_bam = self.align_init()[1]
        return map_bam


    def get_log(self):
        """Return the mapping log
        return: path
        """
        map_log = self.align_init()[2]
        return map_log


    def get_unmap(self):
        """Return the unmap files
        return: path
        """
        unmap_fq = self.align_init()[3]
        return unmap_fq


    def run(self):
        """Map multiple fastq files to one index
        output directory
        """
        args = self.args.copy()
        index_name = AlignIndex(**args).get_index_name()
        logging.info('index: %s' % args['index']) # !!!! logging
        logging.info('index_name: %s' % index_name) # !!!! logging

        aligner_exe = self.aligner_exe
        map_bam, unmap_file = aligner_exe()

        # sort bam

        # index bam
        BAM(map_bam).index()

        return [map_bam, unmap_file]


class AlignHub(object):
    """Alignment 1 fastq file to multiple indexes
    1 fastq
    N index (sequential)
    """
    def __init__(self, fq1, path_out, aligner, index_list,
        smp_name=None, align_by_order=True, dry_run=False, **kwargs):
        """Mulitple indexes"""
        assert isinstance(index_list, list)
        args = args_init(kwargs, align=True)

        ## udpate arguments
        args['fq1'] = fq1
        args['path_out'] = path_out
        args['aligner'] = aligner
        args['smp_name'] = smp_name

        self.index_list = index_list
        self.align_by_order = align_by_order
        self.dry_run = dry_run
        self.args = args

        
    def align_se_batch(self):
        """Alignment 1 fastq file to multiple index,
        priority:
        extra_index > genome
        """
        args = self.args.copy()

        ## output
        map_files = []
        for x in self.index_list:
            if x is None:
                continue
            args['index'] = x

            # index order
            x_order = self.index_list.index(x) + 1

            # check exists
            x_bam = AlignNode(**args).get_bam()
            x_unmap = AlignNode(**args).get_unmap()

            if self.dry_run: # do not run alignment
                pass
            else:
                if not os.path.exists(x_bam) or not os.path.exists(x_unmap) or args['overwrite']:
                    x_bam, x_unmap = AlignNode(**args).run()

            ## update arguments
            if self.align_by_order:
                args['fq1'] = x_unmap

            ## save bam files
            map_files.append(x_bam)

        return map_files


    def run(self):
        """run alignment"""
        map_files = self.align_se_batch()

        return map_files


class AlignReader(object):
    """Return the diretory, files of alignment
    BAM
    count

    Example, structure of directory output
    .
    ├── extra
    │   ├── bigWig
    │   ├── count
    │   ├── mapping
    │   └── report
    ├─ gene
    │   ├── bigWig
    │   ├── count
    │   ├── mapping
    │   └── report
    └── transposon
        ├── bigWig
        ├── count
        ├── mapping
        └── report

    
    Example. structure of gene/mapping
    .
    ├── demo_control_rep1
    │   ├── demo_control_rep1.map_dm6
    │   ├── demo_control_rep1.map_MT_trRNA
    │   ├── demo_control_rep1.not_MT_trRNA.map_dm3
    │   ├── demo_control_rep1.not_MT_trRNA.not_dm3.map_MT_trRNA
    │   └── demo_control_rep1.not_MT_trRNA.not_dm3.not_MT_trRNA.map_hg19
    └── demo_control_rep2
        ├── demo_control_rep2.map_dm6
        ├── demo_control_rep2.map_MT_trRNA
        ├── demo_control_rep2.not_MT_trRNA.map_dm3
        ├── demo_control_rep2.not_MT_trRNA.not_dm3.map_MT_trRNA
        └── demo_control_rep2.not_MT_trRNA.not_dm3.not_MT_trRNA.map_hg19

    Example. structure 

    *.bam
    *.json
    *.log

    """

    def __init__(self, x):
        """The path of mapping directory
        mapping/rep1/rep1.map_dm3/rep1.map_dm3.bam
        *.json
        *.log
        """
        assert os.path.isdir(x)

        ## level1
        subdir1 = self.listdir(x) # level 1
        smp_name = self.listdir(x, full_path=False)

        ## level2
        subdir2 = [self.listdir(i, sort_by_len=True) for i in subdir1] # level 2
        map_name = [self.listdir(i, full_path=False, sort_by_len=True) for i in subdir1]

        self.subdir1 = subdir1
        self.subdir2 = subdir2
        self.smp_name = smp_name
        self.map_name = map_name

        ## check the directory
        assert self.check_dir()

        self.root = x


    def check_dir(self):
        """Check the directory 
        The number of BAM and log files in the directory
        """
        bam_files = self.get_bam()

        # number of samples
        smp_num = len(bam_files)

        # number of indexes for each sample
        index_num = [len(i) for i in bam_files]
        index_unique = list(set(index_num))

        # check index number for each sample
        if len(index_unique) == 1:
            return True
        else:
            index_num = [str(i) for i in index_num]
            logging.error('failed, BAM files for each sample: ' + ' '.join(index_num))
            return False


    def listdir(self, x, full_path=True, sort_by_len=False):
        """Return the name of files, dirs in x
        return the full name
        """
        if full_path:
            p = [os.path.join(x, i) for i in os.listdir(x)]
        else:
            p = os.listdir(x)

        if sort_by_len:
            p = sorted(p, key=len)
        else:
            p = sorted(p)

        return p


    def index_id(self, fn):
        """Extract the name of index from filename
        *.map_{index}
        """
        x = os.path.splitext(fn)[1]
        tag = None
        if x.startswith(r'.map_'):
            tag = re.sub(r'.map_', '', x)

        return tag


    def get_sample_name(self):
        """Sample names in level-1 of root
        return the names
        """
        return self.smp_name


    def get_index_name(self):
        """Name of the alignment dir
        map_{index}
        """
        # index_name
        index_name = [[self.index_id(k) for k in i] for i in self.subdir2]

        # check 
        x = index_name[0] #
        if all([i == x for i in index_name]):
            tag = x
        else:
            tag = None
            logging.error('index number not consistent between samples')

        # return 
        return tag


    def get_bam(self):
        """Return BAM files
        self.map/*.bam
        """
        bam_files = []
        for i in self.subdir2:
            i_bam = []
            for k in i:
                k_name = os.path.basename(k)
                bam = os.path.join(k, k_name + '.bam')
                if os.path.exists(bam):
                    i_bam.append(bam)
            bam_files.append(i_bam)

        return bam_files


    def get_log(self):
        """Return the log files
        *.log
        """
        log_files = []
        for i in self.subdir2:
            i_log = []
            for k in i:
                k_name = os.path.basename(k)
                log = os.path.join(k, k_name + '.log')
                if os.path.exists(log):
                    i_log.append(log)
            log_files.append(i_log)

        return log_files


class AlignSummarize(object):
    """Summarize the alignment 
    Number of reads mapped, unmapped, ...
    multiple index

    format-1:
    genome+spikein:
    rRNA1,genome,rRNA2,spikein

    format-2:
    extra,te,...
    list_in_paralle
    """

    def __init__(self, path_out, **kwargs):
        """The directory of alignemnt"""
        self.path_out = path_out
        self.sub_dirs = self.list_dir()


    def list_dir(self):
        """List the sub-folders, indexes,"""
        sub_dirs = [os.path.join(self.path_out, i) for i in os.listdir(self.path_out)]
        # sub_dirs = [os.path.join(self.path_out, i) for i in sub_dirs if os.path.isdir(i)] # directory
        sub_dirs = [i for i in sub_dirs if os.path.isdir(i)]
        sub_dirs = sorted(sub_dirs, key=len) # sort by length
        return sub_dirs


    def align_type(self):
        """Check the type of alignment:
        1. genome + spikein
        2. te, extra, ...
        
        input: ../arguments.pickle
        """
        args_pickle = os.path.join(os.path.dirname(self.path_out), 'arguments.pickle')
        if not os.path.exists(args_pickle):
            raise Exception('not a alignment directory: %s' % args_pickle)
        # args = pickle.load(args_pickle) # dict
        # read arguments
        with open(args_pickle, 'rb') as fi:
            args = pickle.load(fi)
        if args['align_to_te'] or args['extra_index']:
            flag = 2
        else:
            flag = 1
        return flag


    def align_reads(self, x):
        """total, unique, multiple, map, unmap
        x, the path to alignment directory, (sub_dir)
        require *.json file
        """
        x_name = os.path.basename(x)
        json_file = os.path.join(x, x_name + '.json')
        if not os.path.exists(json_file):
            raise Exception('json file not found: %s' % json_file)
        ## load
        dd = Json_file(json_file).json_reader()

        ## return
        return dd
 

    def index_name(self, x):
        """Return the name of index, at the tail of dirname
        .map_dm3
        .map_rRNA
        .map_extra
        ...
        get the name of index
        """
        x_ext = os.path.splitext(x)[1] # 
        if x_ext.startswith('.map_'):
            flag = re.sub(r'.map_', '', x_ext)
        else:
            flag = None
            logging.warning('not a alignment directory: %s' % x)

        return flag


    def align_sum1(self):
        """Summarize the alignment log/json file
        genome + spikein
        rRNA, genome, rRNA, spikein
        uniform the format
        """
        dirs = self.list_dir()
        # length: 2, 4
        if len(dirs) == 1:
            ## 1st, genome
            d_name1 = self.index_name(dirs[0])
            # label
            # total, rRNA, genome, rRNA, genome, unmap
            k_labels = ['00_total',]
            k_label1 = '%02d_%s_rRNA' % (1, d_name1)
            k_label2 = '%02d_%s' % (2, d_name1)
            k_label3 = '%02d_%s_rRNA' % (3, d_name1)
            k_label4 = '%02d_%s' % (4, d_name1)
            #
            for k in [k_label1, k_label2, k_label3, k_label4]:
                k_labels.append(k + '_unique')
                k_labels.append(k + '_multiple')
            k_labels.append('05_unmap')
            # count
            d_reads1 = self.align_reads(dirs[0])
            k_reads = [d_reads1['total'], 
                0, 0, 
                d_reads1['unique'], d_reads1['multiple'], 
                0, 0, 
                0, 0,
                d_reads1['unmap']]
        elif len(dirs) == 2:
            ## 1st, genome
            ## 2nd, spikein
            d_name1 = self.index_name(dirs[0])
            d_name2 = self.index_name(dirs[1])
            # label
            # total, rRNA, genome, rRNA, genome, unmap
            k_labels = ['00_total',]
            k_label1 = '%02d_%s_rRNA' % (1, d_name1)
            k_label2 = '%02d_%s' % (2, d_name1)
            k_label3 = '%02d_%s_rRNA' % (3, d_name2)
            k_label4 = '%02d_%s' % (4, d_name2)
            #
            for k in [k_label1, k_label2, k_label3, k_label4]:
                k_labels.append(k + '_unique')
                k_labels.append(k + '_multiple')
            k_labels.append('05_unmap')
            # count
            d_reads1 = self.align_reads(dirs[0])
            d_reads2 = self.align_reads(dirs[1])
            k_reads = [d_reads1['total'], 
                0, 0, 
                d_reads1['unique'], d_reads1['multiple'], 
                0, 0, 
                d_reads2['unique'], d_reads2['multiple'],
                d_reads2['unmap']]
        elif len(dirs) == 4:
            d_name1 = self.index_name(dirs[0])
            d_name2 = self.index_name(dirs[1])
            d_name3 = self.index_name(dirs[2])
            d_name4 = self.index_name(dirs[3])
            # label
            k_labels = ['00_total',]
            k_label1 = '%02d_%s_%s' % (1, d_name2, d_name1)
            k_label2 = '%02d_%s' % (1, d_name2)
            k_label3 = '%02d_%s_%s' % (1, d_name4, d_name3)
            k_label4 = '%02d_%s' % (1, d_name4)
            # 
            for k in [k_label1, k_label2, k_label3, k_label4]:
                k_labels.append(k + '_unique')
                k_labels.append(k + '_multiple')
            k_labels.append('05_unmap')
            # count
            d_reads1 = self.align_reads(dirs[0])
            d_reads2 = self.align_reads(dirs[1])
            d_reads3 = self.align_reads(dirs[2])
            d_reads4 = self.align_reads(dirs[3])
            k_reads = [d_reads1['total'],
                d_reads1['unique'], d_reads1['multiple'],
                d_reads2['unique'], d_reads2['multiple'],
                d_reads3['unique'], d_reads3['multiple'],
                d_reads4['unique'], d_reads4['multiple'],
                d_reads4['unmap']]
        else:
            k_labels = ['null']
            k_reads = [0]

        # build dict/pandas.DataFrame
        k_dict = dict(zip(k_labels, k_reads))
        return k_dict


    def align_sum2(self):
        """extra mapping reads"""
        dirs = self.list_dir()
        
        # labels
        k_labels = ['00_total']
        k_reads = []

        for d in dirs:
            d_index = dirs.index(d)
            d_name = self.index_name(d)
            d_reads = self.align_reads(d)
            d_labels = ['%02d_%s_%s' % (d_index, d_name, i) for i in ['unique', 'multiple']]
            d_counts = [d_reads['unique'], d_reads['multiple']]
            if d_index == 0: # 1st
                k_reads.append(d_reads['total'])
            k_labels.extend(d_labels)
            k_reads.extend(d_counts)
        # last one
        k_labels.append('unmap')
        k_reads.append(d_reads['unmap'])

        # build dict
        k_dict = dict(zip(k_labels, k_reads))
        return k_dict


    def run(self):
        flag = self.align_type()
        if flag == 1: 
            d = self.align_sum1()
        else:
            d = self.align_sum2()
        # save dict to file
        map_stat_file = self.path_out.rstrip('/') + '.mapping_stat.csv'
        prefix = os.path.basename(self.path_out)
        df = pd.DataFrame(d, index=[prefix])
        df.to_csv(map_stat_file, index=True)

        return d
 

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

