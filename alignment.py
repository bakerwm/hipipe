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
import io
import glob
import json
import fnmatch
import tempfile
import shlex
import subprocess
import logging
import pandas as pd
import pysam
import pybedtools

from utils_parser import Json_file
from helper import *


class Alignment(object):
    """Run alignment for SE reads using bowtie/bowtie2/STAR"""

    def __init__(self, fqs, path_out, smp_name, genome, **kwargs):
        """parse arguments
        required arguments: fqs, path_out, smp_name, genome, genome, spikein,
        index_ext, threads, unique_only, n_map, aligner, align_to_rRNA,
        genome_path, overwrite
        """
        self.fqs = fqs
        self.path_out = path_out
        self.smp_name = smp_name
        self.genome = genome
        self.kwargs = kwargs
        self.args = self._args_init()


    def _args_init(self):
        """Inititate the arguments, assign the default values to arg
        """
        args = self.kwargs
        args['fqs'] = self.fqs
        args['path_out'] = self.path_out
        args['smp_name'] = self.smp_name
        args['genome'] = self.genome
        args['spikein'] = args.get('spikein', None)
        args['index_ext'] = args.get('index_ext', None)
        args['threads'] = args.get('threads', 1)
        args['unique_only'] = args.get('unique_only', False)
        args['n_map'] = args.get('n_map', 0)
        args['aligner'] = args.get('aligner', 'bowtie')
        args['align_to_rRNA'] = args.get('align_to_rRNA', True)
        args['genome_path'] = args.get('genome_path', None)
        args['overwrite'] = args.get('overwrite', False)
        # check
        if args['spikein'] == self.genome:
            args['spikein'] = None #
        return args


    def _path_init(self, fq, index, reference=None, align_path=None):
        """Create folders for the alignment, 
        Alignment, genome versions
        1.genome_rRNA
        2.genome
        3.spikein_rRNA
        4.spikein

        return files:
        prefix, bam, bed, log, unmap
        """
        args = self.args

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
        # map_bed = map_prefix + '.bed'
        map_log = map_prefix + '.%s.log' % args['aligner']
        unmap_fq = unmap_prefix + '.%s' % fq_type

        return [fq_prefix, map_bam, map_log, unmap_fq, reference]

    
    def _index_builder(self, rRNA=False, genome=None):
        """Return the genome index
        """
        args = self.args
        aligner = args['aligner']

        if not genome:
            genome = self.genome
        # check aligner
        if aligner == 'bowtie':
            index = Genome(genome).bowtie_index(rRNA=rRNA)
        elif aligner == 'bowtie2':
            index = Genome(genome).bowtie2_index(rRNA=rRNA)
        elif aligner == 'STAR':
            index = Genome(genome).star_index(rRNA=rRNA)
        else:
            logging.error('unknown aligner: %s' % aligner)
            index = None # unknonwn aligner
        return index


    def _index_list(self):
        """List the align index (es) for the job
        rRNA, reference, genome
        """
        args = self.args

        # aligner index
        idx = {
            'genome_rRNA': self._index_builder(rRNA=True),
            'genome': self._index_builder(rRNA=False),
            'sp_rRNA': self._index_builder(rRNA=True, genome=args['spikein']),
            'sp': self._index_builder(rRNA=False, genome=args['spikein'])}

        # determine
        if not args['align_to_rRNA']:
            idx['genome_rRNA'] = None

        # save in dict
        return idx # dictionary


    def wrap_log(self, log):
        """Wrapper alignment log file, save as json"""
        args = self.args
        j_file = Alignment_log(log, args['unique_only']).saveas() # save as json


    def bowtie_se(self, fq, index, reference=None, unique_map=False, align_path=None):
        """Run bowtie

        arguments:
        reference: genome, genome_rRNA, spikein, spikein_rRNA
        """
        args = self.args

        bowtie_exe = which('bowtie')
        if not align_path:
            align_path = args['path_out']

        # output directory
        prefix, map_bam, map_log, unmap, reference = self._path_init(fq, index, 
            reference, align_path)

        # determine parameters
        if unique_map:
            para_unique = '-m 1'
        else:
            para_unique = '-v 2 -k %s' % args['n_map'] # default: 1

        if seq_type(fq) == 'fasta':
            para_fq = '-f'
        else:
            para_fq = '-q'

        # file exists
        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('bam file exists: %s' % map_bam)
        else:
            c1 = '%s %s %s -p %s --mm --best --sam --no-unal --un %s %s \
                %s' % (bowtie_exe, para_fq, para_unique, args['threads'], 
                    unmap, index, fq)
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
        return [map_bam, unmap]


    def bowtie2_se(self, fq, index, reference=None, unique_map=False, align_path=None):
        """Run bowtie2

        arguments:
        reference: genome, genome_rRNA, spikein, spikein_rRNA
        """
        args = self.args

        bowtie2_exe = which('bowtie2')
        if not align_path:
            align_path = args['path_out']

        # output directory
        prefix, map_bam, map_log, unmap, reference = self._path_init(fq, index, 
            reference, align_path)
        
        # determine parameters
        if unique_map:
            para_unique = '-q 10'
        else:
            para_unique = '-q 0'
        
        # multi map
        n_map = args['n_map']
        if n_map == 0:
            n_map == 1 # default 1, report 1 hit for each read
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
            c1 = '%s %s -p %s --mm --no-unal --un %s -x %s %s' % (bowtie2_exe, 
                para_fq, args['threads'], unmap, index, fq)
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
        return [map_bam, unmap]    


    def star_se(self, fq, index, reference, unique_map=False, align_path=None):
        """Run STAR, default kwargs
        
        args['unique_only'] is TRUE, unique_map=True:
        """
        args = self.args
        star_exe = which('STAR')
        if not align_path:
            align_path = args['path_out']

        # output directory
        prefix, map_bam, map_log, unmap, reference = self._path_init(fq, index, 
            reference, align_path)
        
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
              --genomeLoad LoadAndKeep \
              --limitBAMsortRAM 10000000000 \
              --outSAMtype BAM SortedByCoordinate \
              --outFilterMismatchNoverLmax 0.07 \
              --seedSearchStartLmax 20 \
              --outReadsUnmapped Fastx %s %s' % (index, fq, fr, map_prefix, 
                args['threads'], unmap, para_unique)
            p1 = subprocess.run(shlex.split(c1))
            
            # filter unique mapped reads
            if unique_map: # only unique mapped reads, -q 10
                print('AAAA')
                pysam.view('-bhS', '-q', '10', '-@', str(args['threads']),
                    '-o', map_bam, map_prefix + 'Aligned.sortedByCoord.out.bam',
                    catch_stdout=False)
            else:
                print('BBBB')
                os.rename(map_prefix + 'Aligned.sortedByCoord.out.bam', map_bam)
            os.rename(map_prefix + 'Unmapped.out.mate1', unmap)
            os.rename(map_prefix + 'Log.final.out', map_log)

        # process log file
        self.wrap_log(map_log)
        return [map_bam, unmap]


    def align_se_batch(self, fq, align_path=None):
        """Align reads to multiple indexes in specific order,
        return align.stat, map_bam, unmap_reads
        determine the index order
        return bam_files
        """
        args = self.args

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
        index_dict = self._index_list() # genome_rRNA, genome, sp_rRNA, sp
        bam_files = []
        fq_input = fq

        # 1. genome_rRNA
        idx1 = index_dict['genome_rRNA']
        reference = self.genome + '_rRNA'
        bam_idx1, unmap_idx1 = aligner_exe(fq=fq_input, index=idx1, 
            reference=reference, unique_map=args['unique_only'], 
            align_path=align_path)
        fq_input = unmap_idx1

        # 2. genome
        idx2 = index_dict['genome']
        reference = self.genome
        bam_idx2, unmap_idx2 = aligner_exe(fq_input, idx1, reference, 
            args['unique_only'], align_path)
        fq_input = unmap_idx2

        if args['spikein']: # add spikein
            # 3. sp_rRNA (optional)
            idx3 = index_dict['sp_rRNA']
            reference = args['spikein'] + '_rRNA'
            bam_idx3, unmap_idx3 = aligner_exe(fq_input, idx1, reference, 
                args['unique_only'], align_path)
            fq_input = unmap_idx3

            # 4. sp (optional)
            idx3 = index_dict['sp']
            reference = args['spikein']
            bam_idx4, unmap_idx4 = aligner_exe(fq_input, idx1, reference, 
                args['unique_only'], align_path)
            fq_input = unmap_idx4

            bam_files = [bam_idx1, bam_idx2, bam_idx3, bam_idx4]
        else: # no spikein
            bam_idx3 = bam_idx4 = None

        bam_files = [bam_idx1, bam_idx2, bam_idx3, bam_idx4]

        return bam_files


    def run(self):
        """Run the alignment for specific fastq file
        1. run alignment for each replicate
        2. merge replicates
        3. run log parser, in json format
        4. organize the log files, saved in one report, including the following groups:
        genome_rRNA, genome_unique, genome_multi, sp_rRNA, sp_unique, sp_multi, unmap
        """
        args = self.args

        bam_out = []
        # run alignment for replicates
        for fq in args['fqs']:
            logging.info('alignment: %s' % fq)
            fq_prefix = file_prefix(fq)[0]
            fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
            fq_path = os.path.join(args['path_out'], fq_prefix)
            assert is_path(fq_path)
            bam_files = self.align_se_batch(fq, fq_path)
            bam_out.append(bam_files) # 
            # stat alignment
            Alignment_stat(fq_path).saveas()


        # merge bam files
        merged_path = os.path.join(args['path_out'], args['smp_name'])
        merged_files = []
        assert is_path(merged_path)
        if len(bam_out) > 1: # for multiple bam files
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

        return genome_bam_files



# class Alignment(object):
#     """Run alignment using bowtie/bowtie2/STAR, SE reads"""

#     def __init__(self, fqs, path_out, smp_name, genome, spikein=None,
#                  index_ext=None, multi_cores=1, unique_only=False, 
#                  n_map=0, aligner='bowtie', align_to_rRNA=False, 
#                  path_data=None, overwrite=False):
#         self.fqs = fqs
#         self.path_out = path_out
#         self.smp_name = smp_name
#         self.genome = genome
#         self.index_ext = index_ext
#         self.spikein = spikein
#         self.multi_cores = multi_cores
#         self.unique_only = unique_only
#         self.n_map = n_map
#         self.aligner = aligner
#         self.align_to_rRNA = align_to_rRNA
#         self.path_data = path_data
#         self.overwrite = overwrite
#         ## wrapper parameters
#         self.args = self.get_args()
        


#     def get_args(self):
#         """Parse parameters"""
#         args = {'fqs' : self.fqs,
#                 'path_out' : self.path_out,
#                 'smp_name' : self.smp_name,
#                 'genome' : self.genome,
#                 'spikein' : self.spikein,
#                 'index_ext' : self.index_ext,
#                 'multi_cores' : self.multi_cores,
#                 'unique_only' : self.unique_only,
#                 'n_map' : self.n_map,
#                 'aligner' : self.aligner,
#                 'align_to_rRNA' : self.align_to_rRNA,
#                 'path_data' : self.path_data,
#                 'overwrite' : self.overwrite
#                 }
#         assert isinstance(args['fqs'], list)
#         assert is_path(args['path_out'])
#         assert isinstance(args['smp_name'], str)
#         assert isinstance(args['genome'], str)
#         assert isinstance(args['multi_cores'], int)
#         assert isinstance(args['unique_only'], bool)
#         assert isinstance(args['aligner'], str)
#         assert isinstance(args['align_to_rRNA'], bool)
#         assert isinstance(args['overwrite'], bool)
#         return args



#     def get_align_index(self):
#         """Get genome index"""
#         args = self.args

#         # spikein index
#         if args['spikein'] == args['genome'] or args['spikein'] is None:
#             idxes = []
#         else:
#             sp_idx = idx_picker(
#                 args['spikein'], 
#                 path_data=args['path_data'], 
#                 aligner=args['aligner'])
#             idxes = [sp_idx, ]

#         # genome index
#         if args['align_to_rRNA'] is True:
#             sg_idx = idx_grouper(
#                 args['genome'],
#                 path_data=args['path_data'], 
#                 aligner=args['aligner'])
#             idxes.extend(sg_idx) # add genome list
#         else:
#             sg_idx = idx_picker(
#                 args['genome'], 
#                 path_data=args['path_data'], 
#                 aligner=args['aligner'])
#             idxes.append(sg_idx) # add genome item

#         # remove possible duplicates, zero-records
#         idxes = list(filter(None.__ne__, idxes)) # remove None
#         idxes = list(dict.fromkeys(idxes)) # keep orders
#         if len(idxes) == 0:
#             raise ValueError('genome index not exists: %s' % args['path_data'])
#         return idxes



#     def init_dir(self, fq, idx, fq_path=None):
#         """Prepare directory, file name for each fastq and index"""
#         assert os.path.exists(fq)
#         assert isinstance(idx, str)
#         args = self.args
#         if fq_path is None:
#             fq_path = args['path_out']

#         fq_prefix = file_prefix(fq)[0] #
#         fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
#         fq_type = seq_type(fq)
#         idx_name = os.path.basename(idx)
#         # fq_path = os.path.join(args['path_out'], fq_prefix)
#         assert is_path(fq_path) # create sub-directory
        
#         map_suffix = os.path.join(fq_path, '%s.map_%s' % (fq_prefix, idx_name))
#         unmap_suffix = os.path.join(fq_path, '%s.not_%s' % (fq_prefix, idx_name))
#         map_bam = map_suffix + '.bam'
#         map_bed = map_suffix + '.bed'
#         map_log = map_suffix + '.%s.log' % args['aligner']
#         unmap_fq = unmap_suffix + '.%s' % fq_type
#         return [fq_prefix, map_bam, map_bed, map_log, unmap_fq]       



#     def bam_index(self, bam):
#         """Run index"""
#         assert isinstance(bam, str)
#         assert os.path.exists(bam)
#         bai = bam + '.bai'
#         if self.overwrite is False and os.path.exists(bai):
#             pass
#         else:
#             pysam.index(bam)
#         return bai



#     def bam_to_bed(self, bam):
#         """Convert bam to bed"""
#         assert isinstance(bam, str)
#         assert os.path.exists(bam)
#         bed = os.path.splitext(bam)[0] + '.bed'
#         if self.overwrite is False and os.path.exists(bed):
#             pass
#         else:
#             pybedtools.BedTool(bam).bam_to_bed().saveas(bed)
#         return bed



#     def wrap_log(self, log):
#         """Wrapper alignment log file, save as json"""
#         assert isinstance(log, str)
#         assert os.path.getsize(log) > 0
#         args = self.args
#         j_file = Alignment_log(log, args['unique_only']).saveas() # save as json



#     def bowtie_se(self, fq, idx, fq_path=None):
#         """Run bowtie, default kwargs"""
#         assert os.path.exists(fq)
#         assert isinstance(idx, str)
#         args = self.args
#         prefix, map_bam, map_bed, map_log, unmap_fq = self.init_dir(fq, idx, fq_path)

#         para_fq = '-f' if seq_type(fq) == 'fasta' else '-q'
#         n_map = args['n_map']
#         if args['unique_only'] is True:
#             para_bowtie = '-m 1'
#         else:
#             n_map = 1 if n_map == 0 else n_map
#             para_bowtie = '-k %s' % n_map

#         if os.path.exists(map_bam) and args['overwrite'] is False:
#             logging.info('file exists, alignment skipped: %s' % map_bam)
#         else:
#             c1 = 'bowtie %s %s -p %s --mm --best --sam --no-unal --un %s %s \
#                   %s' % (para_fq, para_bowtie, args['multi_cores'], unmap_fq, 
#                          idx, fq)
#             c2 = 'samtools view -bhS -F 0x4 -@ %s -' % args['multi_cores']
#             c3 = 'samtools sort -@ %s -o %s -' % (args['multi_cores'], map_bam)
#             with open(map_log, 'wt') as ff:
#                 p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE,
#                                       stderr=ff)
#                 p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout,
#                                       stdout=subprocess.PIPE)
#                 p3 = subprocess.Popen(shlex.split(c3), stdin=p2.stdout)
#                 px = p3.communicate()

#         # process log file
#         self.wrap_log(map_log)
#         return [map_bam, unmap_fq]



#     def bowtie2_se(self, fq, idx, fq_path=None):
#         """Run bowtie2, default kwargs"""
#         assert os.path.exists(fq)
#         assert isinstance(idx, str)
#         args = self.args
#         prefix, map_bam, map_bed, map_log, unmap_fq = self.init_dir(fq, idx, fq_path)

#         para_fq = '-f' if seq_type(fq) == 'fasta' else '-q'
#         n_map = args['n_map']
#         para_bowtie2 = ''
#         if args['unique_only'] is True:
#             para_bowtie2 = '-q 10'
#         else:
#             n_map = 1 if n_map == 0 else n_map
#             para_fq += ' -k %s' % n_map

#         if os.path.exists(map_bam) and args['overwrite'] is False:
#             logging.info('file exists, alignemnt skipped: %s' % map_bam)
#         else:
#             c1 = 'bowtie2 %s -p %s --mm --no-unal --un %s \
#                   -x %s %s' % (para_fq, args['multi_cores'], unmap_fq, idx, fq)
#             c2 = 'samtools view -bhS -F 0x4 -@ %s %s -' % (args['multi_cores'],
#                                                            para_bowtie2)
#             c3 = 'samtools sort -@ %s -o %s -' % (args['multi_cores'], map_bam)
#             with open(map_log, 'wt') as ff:
#                 p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE, 
#                                       stderr=ff)
#                 p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, 
#                                       stdout=subprocess.PIPE)
#                 p3 = subprocess.Popen(shlex.split(c3), stdin=p2.stdout)
#                 px = p3.communicate()

#         # process log file
#         self.wrap_log(map_log)
#         return [map_bam, unmap_fq]



#     def star_se(self, fq, idx, fq_path=None):
#         """Run STAR, default kwargs"""
#         assert os.path.exists(fq)
#         assert isinstance(idx, str)
#         args = self.args
#         prefix, map_bam, map_bed, map_log, unmap_fq = self.init_dir(fq, idx, fq_path)

#         n_map = args['n_map']
#         if args['unique_only'] is True:
#             n_map = 1
#         else:
#             n_map = 20 if n_map == 0 else n_map
#         para_star = '--outFilterMultimapNmax %s' % n_map

#         # if args['unique_only'] is True:
#         #     para_star = '--outFilterMismatchNoverLmax 0.07 --outFilterMultimapNmax 1'

#         freader = 'zcat' if is_gz(fq) else '-'
#         map_prefix = os.path.splitext(map_bam)[0]

#         if os.path.exists(map_bam) and args['overwrite'] is False:
#             logging.info('file exists, alignemnt skipped: %s' % map_bam)
#         else:
#             c1 = 'STAR --runMode alignReads \
#               --genomeDir %s \
#               --readFilesIn %s \
#               --readFilesCommand %s \
#               --outFileNamePrefix %s \
#               --runThreadN %s \
#               --limitOutSAMoneReadBytes 1000000 \
#               --genomeLoad LoadAndRemove \
#               --limitBAMsortRAM 10000000000 \
#               --outSAMtype BAM SortedByCoordinate \
#               --outFilterMismatchNoverLmax 0.07 \
#               --seedSearchStartLmax 20 \
#               --outReadsUnmapped Fastx %s' % (idx, fq, freader, map_prefix,
#                                               args['multi_cores'], para_star)
#             p1 = subprocess.run(shlex.split(c1))
#             # rename exists file
#             # sort bam files
#             pysam.sort('-o', map_bam, 
#                 '-@', str(args['multi_cores']),
#                 map_prefix + 'Aligned.sortedByCoord.out.bam')
#             # os.rename(map_prefix + 'Aligned.sortedByCoord.out.bam', map_bam)
#             os.rename(map_prefix + 'Unmapped.out.mate1', unmap_fq)
#             os.rename(map_prefix + 'Log.final.out', map_log)
            
#         # process log file
#         self.wrap_log(map_log)
#         return [map_bam, unmap_fq]



#     def align_se_batch(self, fq, idxes, fq_path=None):
#         """Run alignment in batch mode, for multiple genome indexes
#         one fastq file to mulitple indexes
#         """
#         assert os.path.exists(fq)
#         assert isinstance(idxes, list)
#         args = self.args
#         if fq_path is None:
#             fq_path = args['path_out']

#         # choose aligner
#         if args['aligner'].lower() == 'star':
#             align_se = self.star_se
#         elif args['aligner'].lower() == 'bowtie2':
#             align_se = self.bowtie2_se
#         elif args['aligner'].lower() == 'bowtie':
#             align_se = self.bowtie_se
#         else:
#             raise ValueError('unknown aligner: %s' % args['aligner'])

#         # iterate indexes
#         bam_files = []
#         fq_input = fq
#         for idx in idxes:
#             bam_map, fq_unmap = align_se(fq_input, idx, fq_path)
#             fq_input = fq_unmap
#             bam_files.append(bam_map)

#         # output
#         return bam_files



#     def run(self):
#         """Run alignments"""
#         args = self.args
#         fqs = args['fqs']
        
#         # check indexes
#         if args['index_ext'] is None:
#             idxes = self.get_align_index()
#         else:
#             idxes = args['index_ext']
#             idxes_flag = [is_idx(f, args['aligner']) for f in idxes]
#             if all(idxes_flag) is False:
#                 logging.error('%10s | -x option error' % 'failed')
#                 logging.error('%10s : %s' % ('Aligner', args['aligner']))
#                 logging.error('%10s : %s' % ('index', idxes))
#                 sys.exit(1)

#         # run alignment
#         out_bam_files = []
#         for fq in fqs:
#             logging.info('mapping: %s ' % fq)
#             fq_prefix = file_prefix(fq)[0]
#             fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
#             fq_path = os.path.join(args['path_out'], fq_prefix)
#             bam_files = self.align_se_batch(fq, idxes, fq_path)
#             out_bam_files.append(bam_files) # bam files
#             # rep path
#             bam_path_out = os.path.dirname(bam_files[0])
#             Alignment_stat(bam_path_out).saveas()

#         # merge bam files
#         merged_path_out = os.path.join(args['path_out'], args['smp_name'])
#         merged_bam_files = [] # map to multiple indexes
#         if len(out_bam_files) > 1: #
#             assert is_path(merged_path_out)
#             for i in range(len(out_bam_files[0])): # 
#                 se_bam_files = [b[i] for b in out_bam_files]
#                 merged_suffix = str_common(se_bam_files, suffix=True)
#                 merged_suffix = re.sub('^_[12]|_R[12]', '', merged_suffix)
#                 merged_bam_name = args['smp_name'] + merged_suffix
#                 merged_bam_file = os.path.join(merged_path_out, merged_bam_name)
#                 merged_bed_file = re.sub('.bam$', '.bed', merged_bam_file)
#                 if os.path.exists(merged_bam_file) and args['overwrite'] is False:
#                     logging.info('file exists: %s' % merged_bam_name)
#                 else:
#                     tmp = bam_merge(se_bam_files, merged_bam_file)
#                     pybedtools.BedTool(merged_bam_file).bam_to_bed().saveas(merged_bed_file)
#                 merged_bam_files.append(merged_bam_file)
#             # merge_map_wrapper(path_out_merge)
#             Alignment_stat(merged_path_out).saveas()

#             out_bam_files.append(merged_bam_files)

#         # create short names for each genome bam
#         genome_bam_files = [f[-1] for f in out_bam_files] # last one
#         genome_bed_files = []
#         for bam_from in genome_bam_files:
#             bam_to = os.path.join(os.path.dirname(bam_from),
#                                   filename_shorter(bam_from))
#             if not os.path.exists(bam_to):
#                 os.symlink(os.path.basename(bam_from), bam_to)
#             if not os.path.exists(bam_to + '.bai'):
#                 if not os.path.exists(bam_from + '.bai'):
#                     pysam.index(bam_from)
#                 os.symlink(os.path.basename(bam_from) + '.bai',
#                            bam_to + '.bai')
#             # rename *.bed
#             bed_from = os.path.splitext(bam_from)[0] + '.bed'
#             bed_to = os.path.splitext(bam_to)[0] + '.bed'
#             if os.path.exists(bed_from) and not os.path.exists(bed_to):
#                 os.symlink(os.path.basename(bed_from), bed_to)
#             genome_bed_files.append(bed_to)
#         return genome_bam_files



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
            if json_files:
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
        if isinstance(df3, pd.DataFrame):
            sp_rRNA = df3['unique'] + df3['multiple']
        else:
            sp_rRNA = 0

        # spikein
        if isinstance(df4, pd.DataFrame):
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
            logging.error('%10s | not contain mapping files: %s' % ('failed', path))
            return None


    def saveas(self, _out=None):
        """Make a copy of statistics of mapping results"""
        path = self.path
        if _out is None:
            prefix = os.path.basename(path)
            _out = os.path.join(os.path.dirname(path),
                                prefix + '.mapping_stat.csv')
        df = self.stat
        df.to_csv(_out, sep='\t', header=True, index=True)
        
        return _out



## EOF
