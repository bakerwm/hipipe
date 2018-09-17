#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Call RT-stop from bed file
save as BED format

"""

import os
import sys
import re
import shlex
import subprocess
import logging
import pysam
import pybedtools
from helper import *


class RTStop(object):
    """Call RT-stop from BED file
    criteria:
    1-nt upstream of read 
    intersect / union between replicates

    input: a bed file
    output: RTStop in BED format

    >>> from call_RTStop import RTStop

    # call RTStop from one bed file
    >>> df = RTStop('demo.bed').rtstop

    # call RTStop from multiple bed files
    >>> fns = ['f1.bed', 'f2.bed']
    >>> df = RTStop(threshold=1).merge(fns, how='inner').rtstop


    # save to file
    >>> RTStop('demo.bed', threshold=2).saveas('out.RTStop.bed')

    """

    def __init__(self, bed=None, threshold=1):
        """Parsing BED file, extract RTStop, save as BED format"""
        self.bed = bed
        self.threshold = threshold
        if bed is None:
            self.rtstop = pd.DataFrame()
        elif isinstance(bed, RTStop):
            self.rtstop = bed.rtstop
        elif isinstance(bed, pd.DataFrame):
            self.rtstop = bed
        elif isinstance(bed, io.TextIOWrapper):
            self.rtstop = self._bed2rtstop().rtstop
        elif os.path.exists(bed):
            self.rtstop = self._bed2rtstop().rtstop
        else:
            raise ValueError('not suppoted file')



    def _tmp(self):
        """Create temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.tmp',
                                            delete=False)
        tmpfn = tmpfn.name
        return tmpfn



    def saveas(self, _out=None):
        """Make a copy of the rtstop records"""
        if _out is None:
            _out = self._tmp()

        df = self.rtstop
        default_kwargs = dict(sep='\t', header=False, index=False)
        df.to_csv(_out, **default_kwargs)
        return _out




    def to_bed6(self):
        """Convert rtstop to BED6 format"""
        df = self.rtstop
        return self._rtstop2bed(df)


    def _is_bed6(self, df=None):
        """Check DataFrame is BED6 format
        Input dataframe contains at least 6-column
        support extra columns
        """
        if df is None:
            df = self.bed # input

        flag = 0
        if len(df.columns) < 6:
            flag += 1
        if pd.api.types.is_string_dtype(df.iloc[:, 0]) is False:
            flag += 1
        if pd.api.types.is_numeric_dtype(df.iloc[:, 1]) is False:
            flag += 1
        if pd.api.types.is_numeric_dtype(df.iloc[:, 2]) is False:
            flag += 1
        if pd.api.types.is_string_dtype(df.iloc[:, 3]) is False:
            flag += 1
        if pd.api.types.is_numeric_dtype(df.iloc[:, 4]) is False:
            flag += 1
        if pd.api.types.is_string_dtype(df.iloc[:, 5]) is False:
            flag += 1
        
        # summary
        if flag > 0:
            return False
        else:
            return True



    def _rtstop2bed(self, df=None):
        """Convert RTStop to BED format
        input: DataFrame, grouped by 'chr, RTStop strand'
        output: DataFrame, BED6, 'count' in extra column
        """
        if df is None:
            df = self.rtstop # DataFrame, groupby ['chrom', 'RTStop', 'strand']
        df = df.reset_index() # DataFrame
        df.insert(1, 'start', df['RTStop'] - 1)
        df.insert(3, 'name', 'RTStop')
        df.insert(4, 'score', 100)
        df.rename(index=str, columns={'RTStop': 'end'}, inplace=True)
        return RTStop(df)



    def _bed2rtstop(self, bed=None):
        """Convert bed to rtstop, save as BED format"""
        if bed is None:
            bed = self.bed # file, io, df
        threshold = self.threshold #

        if os.path.exists(bed):
            df = pybedtools.BedTool(bed).to_dataframe()
        elif isinstance(bed, io.TextIOWrapper):
            df = '' # pybedtols
        else:
            raise ValueError('unknown values')

        assert self._is_bed6(df) # DataFrame

        ## RT reads
        df['RTStop'] = df.apply(lambda row: row.start - 1 if(row.strand == '+')
                                else row.end + 1, axis=1)
        df['RTStart'] = df['RTStop'] - 1
        df['name'] = df['name'] if 'name' in df.columns else 'RTRead'
        df['score'] = df['score'] if 'score' in df.columns else '100'
        # select 6 columns
        df2 = df[['chrom', 'RTStart', 'RTStop', 'name', 'score', 'strand']]

        ## RT stops
        df_cnt = df2.groupby(['chrom', 'RTStop', 'strand']).size() # Series
        df_cnt = df_cnt[df_cnt >= threshold]
        df_cnt = df_cnt.to_frame() # DataFrame

        return RTStop(df_cnt)



    def _merge_dataframe(self, dfs, how='inner'):
        """Merge a list of DataFrame
        Input: DataFrame, contains only one column
        Output: DataFrame, multiple columns
        """
        assert isinstance(dfs, list)
        if len(dfs) == 1:
            return dfs[0]
        else:
            df = dfs.pop(0) # select the first one
            assert isinstance(df, pd.DataFrame)
            for n in dfs:
                if not isinstance(n, pd.DataFrame):
                    continue # skip non-DataFrame items
                df = pd.merge(df, n, how=how, left_index=True, right_index=True)
            # fill na
            df = df.fillna(0) # replace NaN by 0
            return RTStop(df)



    # def merge(self, fns=None, how='inner'):
    #     """Merge RTStops from multiple BED files"""
    #     df = self.rtstop # DataFrame

    #     if fns is None:
    #         return RTStop(df)
    #     else:
    #         assert isinstance(fns, list)
    #         listn = [self._bed2rtstop(i).rtstop for i in fns]
    #         dfn = self._merge_dataframe(listn, how=how)

    #         # make stats for each replicates
    #         dfn = dfn.rtstop # DataFrame
    #         # change column names
    #         dfn.columns = ['rep_%s' % str(int(i) + 1) for i in range(len(fns))]
    #         rep_sum = dfn.sum(axis=1).astype('int')
    #         rep_mean = dfn.mean(axis=1).map('{:.4f}'.format)
    #         rep_std = dfn.std(axis=1).map('{:.4f}'.format)
    #         dfn = dfn.assign(sum=rep_sum, mean=rep_mean, std=rep_std)
    #         # recover to dataframe
    #         return self._rtstop2bed(dfn)



    def merge(self, fns, how='inner'):
        """Merge RTStops from multiple BED files"""
        assert isinstance(fns, list)
        listn = [self._bed2rtstop(i).rtstop for i in fns]
        dfn = self._merge_dataframe(listn, how=how)

        # make stats for each replicates
        dfn = dfn.rtstop # DataFrame
        dfn.columns = ['rep_%s' % str(int(i) + 1) for i in range(len(fns))]
        rep_sum = dfn.sum(axis=1).astype('int')
        rep_mean = dfn.mean(axis=1).map('{:.4f}'.format)
        rep_std = dfn.std(axis=1).map('{:.4f}'.format)
        dfn = dfn.assign(sum=rep_sum, mean=rep_mean, std=rep_std)

        # recover to dataframe
        return self._rtstop2bed(dfn)

