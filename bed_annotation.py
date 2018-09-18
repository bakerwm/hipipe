#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
get annotation of bed
"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


import os
import sys
import re
import argparse
import pathlib
import tempfile
import pybedtools
import pandas as pd
# import numpy as np
# import goldclip
from goldclip.bin.bed_fixer import Bed_parser



class Bed_anno(object):
    """Annotate BED records using BEDTools/pybedtools
    advance:
    - ANNOVAR, http://annovar.openbioinformatics.org/en/latest/
    - annotatePeaks.pl, from HOMER
    """



    def __init__(self, bed, genome, group='homer', path_data=None):
        self.bed = bed
        self.genome = genome
        self.group = group
        self.path_data = path_data


        if isinstance(bed, Bed_anno):
            self.anno = bed.anno
        elif isinstance(bed, pd.DataFrame):
            self.anno = bed
        elif self._is_bed(bed):
            self.anno = self.bed2anno().anno
        elif self._is_bam(bed):
            self.bed = self._bam2bed(self.bed)
            self.anno = self.bed2anno().anno
            os.remove(self.bed)
        else:
            raise ValueError('unknown file format')



    def _is_bam(self, fn):
        """Check input file is BAM file"""
        if fn.lower().endswith('.bam'):
            return True
        else:
            return False



    def _is_bed(self, fn):
        """Check input file is BED file"""
        if fn.lower().endswith('.bed'):
            return True
        else:
            return False



    def _tmp(self):
        """Create temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.bed',
                                            delete=False)
        return tmpfn.name



    def saveas(self, _out=None):
        """Save DataFrame to file"""
        if _out is None:
            _out = self._tmp()

        df = self.anno # DataFrame
        df.to_csv(_out, sep=' ', header=False, index=False)
        return _out



    def _bam2bed(self, fn, overwrite=False):
        """Convert bam to bed, uising pybedtools"""
        if self._is_bed(fn):
            return fn
        elif self._is_bam(fn):
            fn_name = os.path.splitext(os.path.basename(fn))[0] + '.bed'
            fn_bed = os.path.join(tempfile.mkdtemp(), fn_name)
            pybedtools.BedTool(fn).bam_to_bed().saveas(fn_bed)
            return fn_bed



    def anno_picker(self):
        """Pick annotation files
        genome: [dm3|hg19|mm10]
        group: basic, homer
        """
        group = self.group
        genome = self.genome
        path_data = self.path_data

        if path_data is None:
            path_data = os.path.join(pathlib.Path.home(), 'data', 'genome')
        anno_dir = os.path.join(path_data, genome, 'annotation_and_repeats')
        if group == 'basic':
            group_basic = [
                'genicPiRNA',
                'nonGenicPiRNA',
                'TE',
                'tRNA',
                'rRNA',
                'miRNA',
                'targetscan',
                'sncRNA',
                '3u',
                '5u',
                'exon',
                'intron', 
                'igr']
            anno = [os.path.join(anno_dir, 'basic', genome + '.'  + i + '.bed') for i in group_basic]
        if group == 'homer':
            group_homer = [
                'tts', 
                'rRNA',
                'pseudo', 
                'promoters', 
                'ncRNA', 
                'utr3',
                'utr5', 
                'coding', 
                'introns', 
                'intergenic']
            anno = [os.path.join(anno_dir, 'homer', 'ann.' + i + '.bed') for i in group_homer]
        # retrun only exists
        # anno = [f for f in anno if os.path.exists(f)]
        anno_fixed = []
        for f in anno:
            if os.path.exists(f):
                f_tmp = f + '.tmp'
                if not os.path.exists(f_tmp):
                    os.replace(f, f_tmp)
                    Bed_parser(f_tmp).bed_fixer().saveas(f)
                anno_fixed.append(f)
        if len(anno_fixed) == 0:
            raise ValueError('illegeal annotation files: %s' % anno_dir)
            # sys.exit('annotation bed files not found: %s' % genome)
        return anno_fixed



    def bed2anno(self):
        """Annotate BED file"""
        bed = self.bed

        bed_prefix = os.path.splitext(os.path.basename(bed))[0]
        annos = self.anno_picker()

        df = pd.DataFrame(columns=[bed_prefix])
        da = pybedtools.BedTool(bed)
        da_count = da.count() # all count
        da_bed6 = True if da.to_dataframe().shape[1] >= 6 else False
        for f in annos:
            tag = os.path.basename(f).split(r'.')[-2] # category
            db = pybedtools.BedTool(f)
            if da_bed6 is True:
                a_not_b = da.intersect(db, v=True, s=True)
            else:
                a_not_b = da.intersect(db, v=True)
            df.loc[tag] = da.count() - a_not_b.count()
            da = a_not_b
        # others
        df.loc['other'] = da_count - df.sum(axis=0)
        df['sample'] = bed_prefix
        # index
        df = df.reset_index()
        df.columns = ['group', 'count', 'sample']
        df = df.assign(genome=self.genome)
        df = df[['sample', 'genome', 'group', 'count']]
        return Bed_anno(df, self.genome, self.group, self.path_data)



## EOF
