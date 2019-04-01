#!/usr/bin/env python3

"""Download data from OSS

Date: 2019-03-29
Author: Wang Ming
"""
 
__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2019-03-29"
__version__ = "0.1"

import os
import subprocess
import logging
import argparse


def get_args():
    """Download OSS data
    recognize the data by access_id
    """
    parser = argparse.ArgumentParser(prog='oss_downloader.py', 
                                     description='download data using ossutil')
    parser.add_argument('-i', required=True,
        help='the text file saving id and keys')
    parser.add_argument('-o', required=True,
        help='the output directory saving data')
    parser.add_argument('-k', nargs='+', 
        help='A list of keys from -i file')
    parser.add_argument('-n', type=int, default=2,
        help='the number of jobs for each run, default: [2]')
    parser.add_argument('-e', type=int, default=1,
        help='the endpoint, 1=huadong1, 2=huadong2, 3=huabei1, 4=huabei2, default:[4]')
    parser.add_argument('-t', action='store_true',
        help='do not run the download program, just test')
    parser.add_argument('--ossutil', default=None,
        help='path to the ossutil command, default: None, detect in $ENV')
    args = parser.parse_args()
    return args


class OSS(object):
    """oss download
    OSS Endpoint:
    1. 华东 1: oss-cn-hangzhou.aliyuncs.com
    2. 华东 2: oss-cn-shanghai.aliyuncs.com
    3. 华北 1: oss-cn-qingdao.aliyuncs.com
    4. 华北 2: oss-cn-beijing.aliyuncs.com
    """
    def __init__(self, x, path_out, endpoint=4, jobs=2):
        """x, path to records like this
        AccessKeyId: LTAIB4MLnr0DO0dO
        AccessKeySecret: oOHkqfLddz9ik7PM4KaPpucbXzJcHY
        预设OSS路径: oss://novo-data-nj/customer-iEtiIt6f/
        区域: 华东1(杭州)
        """
        ## endpoint 
        ep = {
            1: 'oss-cn-hangzhou.aliyuncs.com',
            2: 'oss-cn-shanghai.aliyuncs.com',
            3: 'oss-cn-qingdao.aliyuncs.com',
            4: 'oss-cn-beijing.aliyuncs.com',
            }

        self.x = x
        self.path_out = path_out
        self.jobs = jobs
        self.endpoint = ep[endpoint]
        self.ossutil = '/data/biosoft/ossutil_aliyuncs/ossutil'


    def config(self, i, k, oss_dir):
        """Config OSS keys"""
        logging.info('config: %s' % i)
        cmd1 = ' '.join([self.ossutil, 'config -e', self.endpoint, '-i', i, '-k', k])
        cmd2 = ' '.join([self.ossutil, 'ls', oss_dir])
        n1 = self.run_shell_cmd(cmd1)
        n2 = self.run_shell_cmd(cmd2)
        return n2


    def copy_files(self, dir_from, dir_to):
        """Download file to directory"""
        cmd1 = ' '.join([self.ossutil, 'cp -r -f --jobs', str(self.jobs), 
            '--parallel', str(self.jobs), dir_from, dir_to + '/'])
        n2 = self.run_shell_cmd(cmd1)


    def run_shell_cmd(self, cmd): 
        """Run shell command"""
        try:
            p = subprocess.run(cmd, shell=True,
                stdout=subprocess.PIPE,
                # stderr=subprocess.STDOUT,
                universal_newlines=True,
                preexec_fn=os.setsid)
            logging.info('run_shell_cmd: CMD={}'.format(cmd))
            return True
        except:
            # raise Exception('Killed CMD={}\nSTDOUT={}'.format(
            #      cmd, ret))
            logging.Warning('Falied CMD: %s' % cmd)
            return False


    def parse_txt(self):
        """Read the i, k, oss_dir, region, ..."""

        rec = []
        with open(self.x, 'rt') as fi:
            lines = filter(None, (line.rstrip() for line in fi))

            for line in lines:
                if line.startswith('AccessKeyId'):
                    acc_id = line.rstrip().split(' ')[1]
                    acc_key = next(lines).rstrip().split(' ')[1]
                    oss_dir = next(lines).rstrip().split(' ')[1]
                    rec.append([acc_id, acc_key, oss_dir])
        return rec


    def download(self, k=None, test=False):
        """Download files"""
        # data list
        rec = self.parse_txt()

        # iterate
        for a in rec:
            acc_id, acc_key, oss_dir = a
            if isinstance(k, list) and not acc_id in k:
                # logging.info('AccessKeyId sippped: %s' % acc_id)
                continue # skip
            # config
            self.config(acc_id, acc_key, oss_dir)
            # run copy
            if not test:
                self.copy_files(oss_dir, self.path_out)


def main():
    logging.basicConfig(level=logging.DEBUG,
        format = '[%(asctime)s] %(message)s',
        datefmt = '%Y-%m-%d %H:%M:%S')

    args = get_args()
    OSS(args.i, args.o, endpoint=args.e).download(args.k, test=args.t)


if __name__ == '__main__':
    main()


