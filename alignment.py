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
    """Run alignment using bowtie/bowtie2/STAR, SE reads"""

    def __init__(self, fqs, path_out, smp_name, genome, spikein=None,
                 multi_cores=1, unique_only=False, aligner='bowtie', 
                 align_to_rRNA=False, path_data=None, overwrite=False):
        self.fqs = fqs
        self.path_out = path_out
        self.smp_name = smp_name
        self.genome = genome
        self.spikein = spikein
        self.multi_cores = multi_cores
        self.unique_only = unique_only
        self.aligner = aligner
        self.align_to_rRNA = align_to_rRNA
        self.path_data = path_data
        self.overwrite = overwrite
        ## wrapper parameters
        self.args = self.get_args()
        


    def get_args(self):
        """Parse parameters"""
        args = {'fqs' : self.fqs,
                'path_out' : self.path_out,
                'smp_name' : self.smp_name,
                'genome' : self.genome,
                'spikein' : self.spikein,
                'multi_cores' : self.multi_cores,
                'unique_only' : self.unique_only,
                'aligner' : self.aligner,
                'align_to_rRNA' : self.align_to_rRNA,
                'path_data' : self.path_data,
                'overwrite' : self.overwrite
                }
        assert isinstance(args['fqs'], list)
        assert is_path(args['path_out'])
        assert isinstance(args['smp_name'], str)
        assert isinstance(args['genome'], str)
        assert isinstance(args['multi_cores'], int)
        assert isinstance(args['unique_only'], bool)
        assert isinstance(args['aligner'], str)
        assert isinstance(args['align_to_rRNA'], bool)
        assert isinstance(args['overwrite'], bool)
        return args



    def get_align_index(self):
        """Get genome index"""
        args = self.args

        # spikein index
        if args['spikein'] == args['genome'] or args['spikein'] is None:
            idxes = []
        else:
            sp_idx = idx_picker(
                args['spikein'], 
                path_data=args['path_data'], 
                aligner=args['aligner'])
            idxes = [sp_idx, ]

        # genome index
        if args['align_to_rRNA'] is True:
            sg_idx = idx_grouper(
                args['genome'],
                path_data=args['path_data'], 
                aligner=args['aligner'])
            idxes.extend(sg_idx) # add genome list
        else:
            sg_idx = idx_picker(
                args['genome'], 
                path_data=args['path_data'], 
                aligner=args['aligner'])
            idxes.append(sg_idx) # add genome item

        # remove possible duplicates, zero-records
        idxes = list(filter(None.__ne__, idxes)) # remove None
        idxes = list(dict.fromkeys(idxes)) # keep orders
        if len(idxes) == 0:
            raise ValueError('genome index not exists: %s' % args['path_data'])
        return idxes



    def init_dir(self, fq, idx, fq_path=None):
        """Prepare directory, file name for each fastq and index"""
        assert os.path.exists(fq)
        assert isinstance(idx, str)
        args = self.args
        if fq_path is None:
            fq_path = args['path_out']

        fq_prefix = file_prefix(fq)[0] #
        fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
        fq_type = seq_type(fq)
        idx_name = os.path.basename(idx)
        # fq_path = os.path.join(args['path_out'], fq_prefix)
        assert is_path(fq_path) # create sub-directory
        
        map_suffix = os.path.join(fq_path, '%s.map_%s' % (fq_prefix, idx_name))
        unmap_suffix = os.path.join(fq_path, '%s.not_%s' % (fq_prefix, idx_name))
        map_bam = map_suffix + '.bam'
        map_bed = map_suffix + '.bed'
        map_log = map_suffix + '.%s.log' % args['aligner']
        unmap_fq = unmap_suffix + '.%s' % fq_type
        return [fq_prefix, map_bam, map_bed, map_log, unmap_fq]       



    def bam_index(self, bam):
        """Run index"""
        assert isinstance(bam, str)
        assert os.path.exists(bam)
        bai = bam + '.bai'
        if self.overwrite is False and os.path.exists(bai):
            pass
        else:
            pysam.index(bam)
        return bai



    def bam_to_bed(self, bam):
        """Convert bam to bed"""
        assert isinstance(bam, str)
        assert os.path.exists(bam)
        bed = os.path.splitext(bam)[0] + '.bed'
        if self.overwrite is False and os.path.exists(bed):
            pass
        else:
            pybedtools.BedTool(bam).bam_to_bed().saveas(bed)
        return bed



    def wrap_log(self, log):
        """Wrapper alignment log file, save as json"""
        assert isinstance(log, str)
        assert os.path.getsize(log) > 0
        args = self.args
        j_file = Alignment_log(log, args['unique_only']).saveas() # save as json



    def bowtie_se(self, fq, idx, fq_path=None):
        """Run bowtie, default kwargs"""
        assert os.path.exists(fq)
        assert isinstance(idx, str)
        args = self.args
        prefix, map_bam, map_bed, map_log, unmap_fq = self.init_dir(fq, idx, fq_path)

        para_fq = '-f' if seq_type(fq) == 'fasta' else '-q'
        para_bowtie = '-v 2 -m 1' if args['unique_only'] is True else '-v 2 -k 1'

        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('file exists, alignment skipped: %s' % map_bam)
        else:
            c1 = 'bowtie %s %s -p %s --mm --best --sam --no-unal --un %s %s \
                  %s' % (para_fq, para_bowtie, args['multi_cores'], unmap_fq, 
                         idx, fq)
            c2 = 'samtools view -bhS -F 0x4 -@ %s -' % args['multi_cores']
            c3 = 'samtools sort -@ %s -o %s -' % (args['multi_cores'], map_bam)
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



    def bowtie2_se(self, fq, idx, fq_path=None):
        """Run bowtie2, default kwargs"""
        assert os.path.exists(fq)
        assert isinstance(idx, str)
        args = self.args
        prefix, map_bam, map_bed, map_log, unmap_fq = self.init_dir(fq, idx, fq_path)

        para_fq = '-f' if seq_type(fq) == 'fasta' else '-q'
        para_bowtie2 = '-q 10' if args['unique_only'] is True else ''

        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('file exists, alignemnt skipped: %s' % map_bam)
        else:
            c1 = 'bowtie2 %s -p %s --mm --no-unal --un %s \
                  -x %s %s' % (para_fq, args['multi_cores'], unmap_fq, idx, fq)
            c2 = 'samtools view -bhS -F 0x4 -@ %s %s -' % (args['multi_cores'],
                                                           para_bowtie2)
            c3 = 'samtools sort -@ %s -o %s -' % (args['multi_cores'], map_bam)
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



    def star_se(self, fq, idx, fq_path=None):
        """Run STAR, default kwargs"""
        assert os.path.exists(fq)
        assert isinstance(idx, str)
        args = self.args
        prefix, map_bam, map_bed, map_log, unmap_fq = self.init_dir(fq, idx, fq_path)

        para_star = '--outFilterMismatchNoverLmax 0.05 --seedSearchStartLmax 20'
        if args['unique_only'] is True:
            para_star = '--outFilterMismatchNoverLmax 0.07 --outFilterMultimapNmax 1'

        freader = 'zcat' if is_gz(fq) else '-'
        map_prefix = os.path.splitext(map_bam)[0]

        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('file exists, alignemnt skipped: %s' % map_bam)
        else:
            c1 = 'STAR --runMode alignReads \
              --genomeDir %s \
              --readFilesIn %s \
              --readFilesCommand %s \
              --outFileNamePrefix %s \
              --runThreadN %s \
              --limitOutSAMoneReadBytes 1000000 \
              --genomeLoad LoadAndRemove \
              --limitBAMsortRAM 10000000000 \
              --outSAMtype BAM SortedByCoordinate \
              --outReadsUnmapped Fastx %s' % (idx, fq, freader, map_prefix,
                                              args['multi_cores'], para_star)
            p1 = subprocess.run(shlex.split(c1))
            # rename exists file
            os.rename(map_prefix + 'Aligned.sortedByCoord.out.bam', map_bam)
            os.rename(map_prefix + 'Unmapped.out.mate1', unmap_fq)
            os.rename(map_prefix + 'Log.final.out', map_log)
            
        # process log file
        self.wrap_log(map_log)
        return [map_bam, unmap_fq]



    def align_se_batch(self, fq, idxes, fq_path=None):
        """Run alignment in batch mode, for multiple genome indexes
        one fastq file to mulitple indexes
        """
        assert os.path.exists(fq)
        assert isinstance(idxes, list)
        args = self.args
        if fq_path is None:
            fq_path = args['path_out']

        # choose aligner
        if args['aligner'].lower() == 'star':
            align_se = self.star_se
        elif args['aligner'].lower() == 'bowtie2':
            align_se = self.bowtie2_se
        elif args['aligner'].lower() == 'bowtie':
            align_se = self.bowtie_se
        else:
            raise ValueError('unknown aligner: %s' % args['aligner'])

        # iterate indexes
        bam_files = []
        fq_input = fq
        for idx in idxes:
            bam_map, fq_unmap = align_se(fq_input, idx, fq_path)
            fq_input = fq_unmap
            bam_files.append(bam_map)

        # output
        return bam_files



    def run(self):
        """Run alignments"""
        args = self.args
        fqs = args['fqs']
        idxes = self.get_align_index()
        out_bam_files = []
        for fq in fqs:
            logging.info('mapping: %s ' % fq)
            fq_prefix = file_prefix(fq)[0]
            fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
            fq_path = os.path.join(args['path_out'], fq_prefix)
            bam_files = self.align_se_batch(fq, idxes, fq_path)
            out_bam_files.append(bam_files) # bam files
            # rep path
            bam_path_out = os.path.dirname(bam_files[0])
            Alignment_stat(bam_path_out).saveas()

        # merge bam files
        merged_path_out = os.path.join(args['path_out'], args['smp_name'])
        merged_bam_files = [] # map to multiple indexes
        if len(out_bam_files) > 1: #
            assert is_path(merged_path_out)
            for i in range(len(out_bam_files[0])): # 
                se_bam_files = [b[i] for b in out_bam_files]
                merged_suffix = str_common(se_bam_files, suffix=True)
                merged_suffix = re.sub('^_[12]|_R[12]', '', merged_suffix)
                merged_bam_name = args['smp_name'] + merged_suffix
                merged_bam_file = os.path.join(merged_path_out, merged_bam_name)
                merged_bed_file = re.sub('.bam$', '.bed', merged_bam_file)
                if os.path.exists(merged_bam_file) and args['overwrite'] is False:
                    logging.info('file exists: %s' % merged_bam_name)
                else:
                    tmp = bam_merge(se_bam_files, merged_bam_file)
                    pybedtools.BedTool(merged_bam_file).bam_to_bed().saveas(merged_bed_file)
                merged_bam_files.append(merged_bam_file)
            # merge_map_wrapper(path_out_merge)
            Alignment_stat(merged_path_out).saveas()

            out_bam_files.append(merged_bam_files)

        # create short names for each genome bam
        genome_bam_files = [f[-1] for f in out_bam_files] # last one
        genome_bed_files = []
        for bam_from in genome_bam_files:
            bam_to = os.path.join(os.path.dirname(bam_from),
                                  filename_shorter(bam_from))
            if not os.path.exists(bam_to):
                os.symlink(os.path.basename(bam_from), bam_to)
            if not os.path.exists(bam_to + '.bai'):
                if not os.path.exists(bam_from + '.bai'):
                    pysam.index(bam_from)
                os.symlink(os.path.basename(bam_from) + '.bai',
                           bam_to + '.bai')
            # rename *.bed
            bed_from = os.path.splitext(bam_from)[0] + '.bed'
            bed_to = os.path.splitext(bam_to)[0] + '.bed'
            if os.path.exists(bed_from) and not os.path.exists(bed_to):
                os.symlink(os.path.basename(bed_from), bed_to)
            genome_bed_files.append(bed_to)
        return genome_bam_files



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
        elif isinstance(log, io.TextIOWrapper):
            self.stat = self._log_parser()
        elif os.path.isfile(log):
            self.stat = self._log_parser()
        else:
            raise ValueError('not supported file')



    def _is_file(self):
        """Check the log file is exists, not empty
        """
        if os.path.isfile(self.log):
            return True
        else:
            return False



    def _is_non_empty(self):
        """Check if log file is empty"""
        if os.path.getsize(self.log) > 0:
            return True
        else:
            return False



    def _bowtie_log(self):
        """Wrapper bowtie log"""
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
        dd['multiple'] = 0 #no multiple
        # if self.unique_only is True or 'multiple' in dd:
        #     pass
        # else:
        #     dd['multiple'] = 0
        dd['unmap'] = dd['total'] - dd['map']
        return dd



    def _bowtie2_log(self):
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
        if self.unique_only is True:
            dd['map'] = dd['unique']
            # dd['multiple'] = 0
        else:
            # unique and multiple
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['map']
        return dd



    def _star_log(self):
        """Wrapper STAR *Final.log"""
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
            # dd['multiple'] = 0
        else:
            # unique and multiple
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['map']
        return dd


    def _log_parser(self):
        """Read log file as dict
        delimiter:
        bowtie:  ":"
        bowtie2:  ":"
        STAR:  "|"

        extra trimming
        1. trim "(10.00%)" 
        2. trim "blank" at both tails
        """
        log = self.log
        log_lines = []
        if isinstance(log, dict):
            return Alignment_log(log)
        elif isinstance(log, io.TextIOWrapper):
            for r in log:
                if r.startswith('Warning'): # skip warnings
                    continue
                log_lines.append(r.strip())
        elif os.path.exists(log):
            with open(self.log, 'rt') as ff:
                for r in ff:
                    if r.startswith('Warning'):
                        continue
                    log_lines.append(r.strip())
        else:
            raise ValueError('unknown file format')
        
        # parsing log file
        line = log_lines[0] # the first line
        if line.startswith('#'):
            dd = self._bowtie_log()
        elif 'reads; of these' in line:
            dd = self._bowtie2_log()
        elif '|' in line:
            dd = self._star_log()
        else:
            raise ValueError('unknown file format')

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
            json.dump(dd, fo, indent=4, sort_keys=True)

        return _out



class Alignment_stat(object):
    """Parse mapping reads in directory
    1. for each rep bam, parse json files, 
    del. 2. merged bam files, count bam lines
    2. merge replicates
    """
    def __init__(self, path):
        self.path = path

        if isinstance(path, Alignment_stat):
            self.stat = path.stat
        elif isinstance(path, pd.DataFrame):
            self.stat = path
        elif os.path.isdir(path):
            json_files = self._list_files('.json')
            bam_files = self._list_bam_files()
            if not json_files is None:
                self.stat = self.align_log_stat()
            elif not bam_files is None:
                # self.stat = self.bam_stat()
                self.stat = self.merge_rep_stat()
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



    def _list_files(self, suffix='*'):
        """List json files in mapping directory"""
        path = self.path
        json_files = [os.path.join(path, f) for f in os.listdir(path) 
            if f.endswith(suffix)]
        json_files = [f for f in json_files if self._is_non_empty(f) is True]
        if len(json_files) > 0:
            return json_files
        else:
            return None



    def _json_log_reader(self, fn):
        """Parse json file, save as DataFrame"""
        # fn_name = os.path.splitext(os.path.basename(fn))[0]
        fn_name = os.path.basename(fn)
        group, aligner = fn.split('.')[-3:-1]
        group = re.sub('^map_', '', group)
        dd = Json_file(fn).stat # dict of count
        # total, unique, multiple, unmap
        df = pd.DataFrame(data=[list(dd.values())], columns=list(dd.keys()))
        df.pop('map') # drop "map" column
        df.index = [fn_name]
        return df



    # directory of replicate fastq file contains *.json file
    def align_log_stat(self):
        """Extract alignment log from one directory"""
        path = self.path.rstrip('/')
        json_files = self._list_files('.json')

        if len(json_files) == 0:
            return None
        elif len(json_files) == 1:
            # map_genome
            df = self._json_log_reader(log_files[0])
            return df
        else:
            json_files = sorted(json_files, key=len)
            # map to genome
            dfg = self._json_log_reader(json_files[-1])
            dfg = dfg.loc[:, ['unique', 'multiple', 'unmap']]

            # map to others
            frames = []
            for json_file in json_files[:-1]:
                df = self._json_log_reader(json_file)
                log_name = df.index.tolist()[0] # first one
                group = log_name.split('.')[-3] # map_group
                group = re.sub('map_', '', group)
                df = df.assign(mapped=df.unique + df.multiple)
                df.columns.values[-1] = group # rename column
                frames.append(df.loc[:, group])
            # all rows
            dfx = pd.concat(frames, axis=1)
            dfx.index = list(dfg.index.values)
           
            # all stat
            df = pd.concat([dfx, dfg], axis=1)
            df = df.assign(total=pd.Series(df.sum(axis=1)).values)

            # index
            df.index = [os.path.basename(path)]
            return df



    # directory of merged fastq, count bam file
    # to-do: summary replicates (except unmap)
    def _list_bam_files(self):
        """Skip sym-link file"""
        path = self.path
        bam_files = [b for b in sorted(self._list_files('.bam'), key=len) 
            if self._is_non_empty(b) and not os.path.islink(b)]
        if len(bam_files) > 0:
            return bam_files
        else:
            return None



    # def bam_stat(self):
    #     path = self.path
    #     prefix = os.path.basename(path)
    #     bam_files = self._list_bam_files()
        
    #     frames = []
    #     for bam in bam_files:
    #         bam_count = pysam.AlignmentFile(bam, 'rb').count()
    #         group = bam.split('.')[-2] # map_group
    #         group = re.sub('map_', '', group)
    #         if len(bam_files) > 1 and bam_files.index(bam) == 0 and group == 'genome':
    #             group = 'spikein'
    #         df = pd.DataFrame({'data': bam_count}, index=[prefix])
    #         df.columns = [group]
    #         frames.append(df)
    #     # merge
    #     df_all = pd.concat(frames, axis=1)
    #     return df_all



    def merge_rep_stat(self):
        """Stat reads for merged sample
        combine stat of replicates
        search the parent_directory of merged_dir
        """
        path = self.path.rstrip('/')
        parent_path = os.path.dirname(path)
        rep_dirs = [os.path.join(parent_path, f) for f in os.listdir(parent_path)]
        rep_dirs = [f for f in rep_dirs if os.path.isdir(f)]

        frames = []
        for rep_dir in rep_dirs:
            json_files = [os.path.join(rep_dir, f) for f in os.listdir(rep_dir)
                if f.endswith('.json')]
            if len(json_files) > 0:
                frames.append(Alignment_stat(rep_dir).stat)
            else:
                continue # next
        # output
        if len(frames) > 0:
            # rep
            df = pd.concat(frames, axis=0)
            rep_names = list(df.index.values)
            rep_suffix = str_common(rep_names, suffix=True)
            df.index = [f.rstrip(rep_suffix) for f in rep_names]
            # merge
            df_merge = pd.DataFrame(data=[df.sum(axis=0)], columns=list(df.columns.values))
            df_merge.index = [os.path.basename(path)]
            # return pd.concat([df, df_merge], axis=0)
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
        df = df.reset_index()

        default_kwargs = dict(sep='\t', header=True, index=False)
        if isinstance(_out, io.TextIOWrapper):
            print(df.to_string(index=False, header=False, justify='left'))
        else:
            df.to_csv(_out, **default_kwargs)
        
        return _out



## EOF
