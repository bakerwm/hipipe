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

import datetime
import logging

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


## design
class DesignReader(object):
    """Parsing tab file for arguments
    return dict 
    optional:
    return: json file, dict
    """

    def __init__(self, file, j=None):
        self.file = file
        self.j = j
        self.d = self.parseTxt()

    def parseFqfiles(self, x, path):
        """Parse fastq files within path, using the prefix
        PE reads
        SE reads

        *.fastq
        *.fq
        *.fastq.gz
        *.fq.gz

        _1.
        _2.

        """

        all_files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

        # filter fastq
        p1 = re.compile(r'(_[12])?.f(ast)?q(.gz)?$')
        hit_files = [f for f in all_files if p1.search(f) and x in f]
        
        # determine SE or PE
        fq_r1 = []
        fq_r2 = []
        p2 = re.compile(r'_[rR]?2.f(ast)?q(.gz)?$') # read1
        for f in sorted(hit_files):
            if p2.search(f):
                fq_r2.append(f)
            else:
                fq_r1.append(f)

        if len(fq_r1) > 3:
            logging.warning('too many records matched : {} \n{}'.format(x, fq_r1))

        if len(fq_r2) > 0 and not len(fq_r1) == len(fq_r2):
            logging.error('read1 and read2 not matched: \n{}\n{}'.format(fq_r1, fq_r2))

        if len(fq_r1) == 0 and len(fq_r2) == 0:
            logging.error('no fastq files matched: {} {}'.format(x))

        # print('\t'.join(fq_r1))
        return [fq_r1, fq_r2]


    def parseTxt(self):
        """parse txt and create dict

        control_name    treatment_name  output genome

        """
        dd = {}
        with open(self.file, 'rt') as fi:
            nline = 1
            for line in fi:
                if line.strip().startswith('#') or not line.strip():
                    nline += 1
                    continue
                try:
                    ctl, tre, fq_dir, out, genome, spikein = line.strip().split('\t')
                except:
                    logging.error('unknown format in line-{}: {}'.format(nline, line))
                    continue

                # parse fastq files
                k = '{}_vs_{}'.format(ctl, tre)
                ctl_fq = self.parseFqfiles(ctl, fq_dir)
                tre_fq = self.parseFqfiles(tre, fq_dir)

                if ctl == tre or ctl_fq == tre_fq:
                    logging.warning('identical names detected in line-{}: {}'.format(nline, line))
                    continue

                if k in dd:
                    logging.warning('duplicated design in line-{}: {}'.format(nline, line))
                    continue

                # tre_fq = None

                dd[k] = {
                    'control': ctl_fq,
                    'control_name': ctl,
                    'treatment': tre_fq,
                    'treatment_name': tre,
                    'path_out': out,
                    'genome': genome,
                    'spikein': spikein,
                }
                nline +=1

        return dd

    
    def to_dict(self):
        # return self.parseTxt()
        return self.d


    def to_json(self, file=None):
        # d = self.to_dict()

        if file is None:
            tmp_prefix = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
            file = '{}.design.json'.format(tmp_prefix)

        with open(file, 'wt') as fo:
            json.dump(self.d, fo, indent=4, sort_keys=True)

        return file

