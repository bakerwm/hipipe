#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
import io
import re
import json
import glob
import tempfile
import pysam
import pandas as pd


class Cutadapt_log(object):
    """Wrapper cutadapt log file"""

    def __init__(self, log):
        self.log = log
        # stat
        if isinstance(log, Cutadapt_log):
            self.stat = stat.stat
        elif isinstance(log, dict):
            self.stat = stat
        elif isinstance(log, io.TextIOWrapper):
            self.stat = self._log_parser()
        elif os.path.isfile(log):
            self.stat = self._log_parser()
        else:
            raise ValueError('not supported file')


    def _log_parser(self):
        """Wrapper log file"""
        dd = {}
        with open(self.log, 'rt') as ff:
            for line in ff:
                if line.startswith('This is cutadapt'):
                    sep = line.strip().split(' ')
                    dd['version'] = sep[3]
                    dd['python'] = sep[6]
                elif 'Command line parameters' in line:
                    dd['cmd'] = line.strip().split(':')[1]
                elif 'Total reads processed' in line:
                    value = line.strip().split(':')[1]
                    value = re.sub(',', '', value.strip())
                    dd['total'] = int(value)
                elif 'Reads written (passing filters)' in line:
                    value = line.strip().split(':')[1]
                    value = value.strip().split(' ')[0]
                    value = re.sub(',', '', value)
                    dd['clean'] = int(value)
                else:
                    continue
        pct = float(dd['clean']) / float(dd['total']) * 100
        dd['pct'] = '%.1f%%' % pct
        return dd



    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.json',
                                            delete=False)
        return tmpfn.name



    def saveas(self, _out=None):
        """Make a copy of statistics of mapping results"""
        if _out is None:
            _out = os.path.splitext(self.log)[0] + '.json'
            # _out = self._tmp()

        dd = self.stat

        with open(_out, 'wt') as fo:
            json.dump(dd, fo, indent=4)

        return _out



class Json_file(object):
    """Parsing Json and dict file"""

    def __init__(self, fn):
        self.fn = fn
        if isinstance(fn, Json_file):
            self.stat = fn.stat
        elif isinstance(fn, dict):
            self.stat = fn
        elif os.path.exists(fn):
            self.stat = self.json_reader()
        else:
            raise ValueError('unknown file format:')


    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.json',
                                            delete=False)
        return tmpfn.name


    def json_reader(self):
        """Load json file as dict"""
        fn = self.fn
        if os.path.isfile(fn) and os.path.getsize(fn) > 0:
            with open(fn, 'rt') as ff:
                return json.load(ff)



    def json_writer(self, to=None):
        """Write dict to file in json format"""
        fn = self.fn

        if to is None:
            to = self._tmp()

        if isinstance(fn, Json_file):
            fn = fn.fn
        elif isinstance(fn, dict):
            fn = fn
        with open(to, 'wt') as ff:
            json.dump(fn, ff, indent=4, sort_keys=True)






##--------------------##
## figure 1
def trim_wrapper(path, smp_name='demo'):
    """
    trimming and remove duplicates
    input: /path_out/input_reads/
    """
    json_files = sorted(glob.glob(path + '/*.cutadapt.json'))
    da = []
    for j in json_files:
        id = re.sub(r'.cutadapt.json', '', os.path.basename(j))
        nodup = os.path.join(os.path.dirname(j), id + '.reads.txt') # total reads, nodup
        d = json_reader(j)
        with open(nodup) as f:
            d['nodup'] = next(f).rstrip()
        tooshort = int(d['raw']) - int(d['clean'])
        dup = int(d['clean']) - int(d['nodup'])
        dn = pd.DataFrame({'group': ['raw', 'too_short', 'PCR_dup', 'no_dup'],
            id: [d['raw'], tooshort, dup, d['nodup']]})
        dn.set_index('group', inplace = True)
        da.append(dn)
    df = pd.concat(da, axis = 1)
    df = df.apply(pd.to_numeric)
    # add merge data
    df.insert(0, smp_name, df.sum(axis = 1))
    return df

# path_trim = os.path.join(path_out, 'input_reads')
# df = trim_wrapper(path_trim, smp_name)
# print(df)


## mapping pct
def map_wrapper(path, smp_name='demo'):
    """
    mapping to various genome
    input: /path_out/genome_mapping/
    """
    m_files = glob.glob(os.path.join(path, '*.mapping_stat.csv'))
    m_files = sorted(m_files, key=len)
    ma = []
    for m in m_files:
        # skip merge stat
        m_prefix = re.sub(r'.mapping_stat.csv', '', os.path.basename(m))
        if m_prefix == smp_name:
            continue
        dm = pd.read_csv(m, ',').filter(items=['group', 'read'])
        dm.set_index('group', inplace=True)
        dm2 = dm.rename(columns={'read': m_prefix})
        ma.append(dm2)
    df = pd.concat(ma, axis=1)
    # add merge data
    df.insert(0, smp_name, df.sum(axis=1))
    return df







