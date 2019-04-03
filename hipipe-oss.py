#!/usr/bin/env python3
#-*- coding: utf8 -*-

"""
Download OSS objects using ossutil or OSS python SDK package oss2

see the following links for more details:

1. https://help.aliyun.com/document_detail/32027.html
2. https://aliyun-oss-python-sdk.readthedocs.io/en/stable/api.html

"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2019-03-29"
__version__ = "0.1"


import os
import pathlib
import oss2
import argparse
import logging
import subprocess


def get_args():
    """Download OSS data
    recognize the data by access_id
    """
    parser = argparse.ArgumentParser(
        prog='hipipe-oss.py', 
        description='download data using oss2 python script or ossutil')
    parser.add_argument('-i', required=True,
        help='text file saving id and keys')
    parser.add_argument('-o', required=True,
        help='the output directory saving data')
    parser.add_argument('-k', nargs='+', 
        help='A list of keys from -i file')
    parser.add_argument('--ossutil', default=None,
        help='absolute path to ossutil command')
    parser.add_argument('--dry-run', action='store_true', dest='dry_run',
        help='do not run the download program, just test')
    args = parser.parse_args()
    return args


class OSS(object):
    """Download/Upload files to OSS on aliyun"""

    def __init__(self, ki, ks, endpoint, bucket_name, sub_dir=None):
        """Download files/directory from OSS
        
        Data_text:
        AccessKeyId: LTAISMrCJ5zDsDVx
        AccessKeySecret: PiSswEuXgRSBQebBcQh3YwtW9aY9Nc
        预设OSS路径: oss://pangoo-bj/customer-c6k9xwmI/
        区域: 华北2(北京) oss-cn-beijing.aliyuncs.com

        arguments:
        ki: AccessKeyId
        ks: AccessKeySecret
        endpoint:  oss-cn-beijing.aliyuncs.com
        bucket_name: pangoo-bj
        sub_dir: customer-c6k9xwmI/

        """
        self.ki = ki
        self.ks = ks
        self.endpoint = endpoint
        self.bucket_name = bucket_name
        self.sub_dir = sub_dir
        ## get bucket
        self.bucket = self.get_bucket()


    def get_bucket(self):
        """Create bucket Class"""
        auth = oss2.Auth(self.ki, self.ks)
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        return bucket


    def dir_list(self, remote_dir, local_dir, dry_run=False):
        """List all files/directories from root (sub_dir=None) or from sub_dir
        Download data from Sequencing company
        root directory is not readable.
        """
        # level of directory
        remote_dir = remote_dir.rstrip('/')

        oss_file_size = 0
        for obj in oss2.ObjectIterator(self.bucket, prefix = '%s/' % remote_dir):
            oss_file_size = oss_file_size + 1
            print(obj.key)
            # if obj.size > 1410000000:
            #     continue
            if dry_run: # skip download
                continue
            self.download_to_local(self.bucket, obj.key, obj.key, local_dir)

        print('%s : file size : %s' % (remote_dir, oss_file_size))


    def download_to_local(self, bucket, object_name, local_file, local_path):
        """object_name, 
        file
        directory/
        """
        # absolute path
        local_file_path = os.path.join(local_path, local_file)

        # filename + directory
        local_file_name = os.path.basename(local_file_path)
        local_file_dirname = os.path.dirname(local_file_path)

        # create directories
        if os.path.exists(local_file_dirname) is False:
            # pathlib.Path(local_file_dirname).mkdir(parents=True, exist_ok=True)
            oss2.utils.makedir_p(local_file_dirname)

        ## for directory
        if object_name.endswith('/'): # directory
            oss2.utils.makedir_p(local_file_path) # create directory
        else:
            ## download
            print('file: %s' % object_name)
            bucket.get_object_to_file(object_name, local_file_path)


    def download(self, local_path='./', dry_run=False):
        # download direcotry
        if self.sub_dir is None:
            check_dir = '/' # root directory
        else:
            check_dir = self.sub_dir
        self.dir_list(check_dir, local_path, dry_run)


class OSSUtil(object):
    """oss download
    OSS Endpoint:
    1. 华东 1: oss-cn-hangzhou.aliyuncs.com
    2. 华东 2: oss-cn-shanghai.aliyuncs.com
    3. 华北 1: oss-cn-qingdao.aliyuncs.com
    4. 华北 2: oss-cn-beijing.aliyuncs.com
    """
    def __init__(self, acc_id, acc_secret, endpoint, remote_path, jobs=2,
        ossutil='/data/biosoft/ossutil_aliyuncs/ossutil'):
        """Download files using ossutil
        acc_id: AccessKeyId: LTAIB4MLnr0DO0dO
        acc_secret: AccessKeySecret: oOHkqfLddz9ik7PM4KaPpucbXzJcHY
        remote_path: 预设OSS路径: oss://novo-data-nj/customer-iEtiIt6f/
        retion: 区域: 华东1(杭州)
        endpoint: oss-cn-hangzhou.aliyuncs.com
        """
        self.acc_id = acc_id
        self.acc_secret = acc_secret
        self.endpoint = endpoint
        self.remote_path = remote_path
        self.jobs = jobs
        self.ossutil = ossutil


    def config(self):
        """Config OSS keys"""
        logging.info('config AccessKeyId: %s' % self.acc_id)
        cmd1 = ' '.join([self.ossutil, 'config -e', self.endpoint, '-i', self.acc_id, '-k', self.acc_secret])
        # n1 = self.run_shell_cmd(cmd1)
        logging.info('run_shell_cmd: {}'.format(cmd1))
        subprocess.run(cmd1, shell=True)


    def list_files(self):
        """List files in remote directory"""
        logging.info('list files:')
        cmd1 = ' '.join([self.ossutil, 'ls', self.remote_path])
        logging.info('run_shell_cmd: {}'.format(cmd1))
        subprocess.run(cmd1, shell=True)


    def copy_files(self, local_path):
        """Download file to directory"""
        cmd1 = ' '.join([self.ossutil, 'cp -r -f --jobs', str(self.jobs), 
            '--parallel', str(self.jobs), self.remote_path, self.local_path + '/'])
        logging.info('run_shell_cmd: {}'.format(cmd1))
        subprocess.run(cmd1, shell=True)


    # def run_shell_cmd(self, cmd): 
    #     """Run shell command"""
    #     try:
    #         p = subprocess.run(cmd, shell=True,
    #             stdout=subprocess.PIPE,
    #             # stderr=subprocess.STDOUT,
    #             universal_newlines=True,
    #             preexec_fn=os.setsid)
    #         logging.info('run_shell_cmd: CMD={}'.format(cmd))
    #         return True
    #     except:
    #         # raise Exception('Killed CMD={}\nSTDOUT={}'.format(
    #         #      cmd, ret))
    #         logging.Warning('Falied CMD: %s' % cmd)
    #         return False


    def download(self, local_path, dry_run=False):
        """Download files"""
        # config
        self.config()
        
        # run
        if dry_run:
            # list files
            self.list_files()
        else:
            # copy files
            self.copy_files(local_path)


class OSSList(object):

    def __init__(self, x):

        self.x = x # fie
        # self.choose_keys = choose_keys
        """
        Data_text:
        AccessKeyId: LTAISMrCJ5zDsDVx
        AccessKeySecret: PiSswEuXgRSBQebBcQh3YwtW9aY9Nc
        预设OSS路径: oss://pangoo-bj/customer-c6k9xwmI/
        区域: 华北2(北京) oss-cn-beijing.aliyuncs.com
        """

    def parse_txt(self):
        """Read the i, k, oss_dir, region, ..."""

        dd = {}        
        with open(self.x, 'rt') as fi:
            lines = filter(None, (line.rstrip() for line in fi))
            for line in lines:
                if line.startswith('AccessKeyId'):
                    # AccessKeyId: LTAISMrCJ5zDsDVx
                    acc_id = line.rstrip().split(' ')[1]
                    if not acc_id in dd:
                        dd[acc_id] = {} # new dict
                    dd[acc_id]['acc_id'] = acc_id
                if line.startswith('AccessKeySecret'):
                    # AccessKeySecret: PiSswEuXgRSBQebBcQh3YwtW9aY9Nc
                    acc_secret = line.rstrip().split(' ')[1]
                    dd[acc_id]['acc_secret'] = acc_secret
                if 'oss://' in line:
                    # 预设OSS路径: oss://pangoo-bj/customer-c6k9xwmI/
                    folder_url = line.rstrip().split(' ')[1]
                    bucket, folder_name = folder_url.rstrip('/').split('/')[-2:] # 
                    dd[acc_id]['folder_url'] = folder_url
                    dd[acc_id]['bucket'] = bucket
                    dd[acc_id]['folder_name'] = folder_name
                if '区域' in line:
                    # 区域: 华北2(北京)
                    region = line.rstrip().split(' ')[1]
                    region = region.rstrip().split('(')[0]
                    dd[acc_id]['endpoint'] = self.get_endpoint(region)
        return dd


    def get_endpoint(self, region):
        """Determine the endpoint url by name or region
        """
        # see aliyun website:
        # https://www.alibabacloud.com/help/doc-detail/31837.htm
        ep = {
        '华东1': 'oss-cn-hangzhou.aliyuncs.com',
        '华东2': 'oss-cn-shanghai.aliyuncs.com',
        '华北1': 'oss-cn-qingdao.aliyuncs.com',
        '华北2': 'oss-cn-beijing.aliyuncs.com',
        '华北3': 'oss-cn-zhangjiakou.aliyuncs.com',
        '华北5': 'oss-cn-huhehaote.aliyuncs.com',
        '华南1': 'oss-cn-shenzhen.aliyuncs.com',
        'China East 1': 'oss-cn-hangzhou.aliyuncs.com',
        'China East 2': 'oss-cn-shanghai.aliyuncs.com',
        'China North 1': 'oss-cn-qingdao.aliyuncs.com',
        'China North 2': 'oss-cn-beijing.aliyuncs.com',
        'China North 3': 'oss-cn-zhangjiakou.aliyuncs.com',
        'China North 5': 'oss-cn-huhehaote.aliyuncs.com',
        'China South 1': 'oss-cn-shenzhen.aliyuncs.com',
        }

        if region in ep:
            endpoint = ep[region]
        else:
            endpoint = None
        return endpoint

    def get_keys(self, choose_keys=None):
        return self.parse_txt()


def main():
    """Run command"""
    logging.basicConfig(level=logging.DEBUG,
        format = '[%(asctime)s] %(message)s',
        datefmt = '%Y-%m-%d %H:%M:%S')

    args = get_args()

    # list of keys
    dk = OSSList(args.i).get_keys()

    # filter keys
    dk_choose = {}
    if args.k is None:
        dk_choose = dk.copy()
    else:
        for k in args.k:
            if k in dk:
                dk_choose[k] = dk[k]
            else:
                logging.Warning('keys not found in -i: ' + k)

    # download files
    for i in dk_choose:
        di = dk_choose[i] # dict
        acc_id = di['acc_id']
        acc_secret = di['acc_secret']
        endpoint = di['endpoint']
        bucket_name = di['bucket']
        remote_path = di['folder_url']
        subdir_name = di['folder_name']
        if args.ossutil is None:
            OSS(acc_id, acc_secret, endpoint, bucket_name, subdir_name).download(args.o, args.dry_run)
        else:
            if os.path.exists(args.ossutil):
                OSSUtil(acc_id, acc_secret, endpoint, remote_path).download(args.o, args.dry_run)
            else:
                logging.error('command not exists: %s' % args.ossutil)


if __name__ == '__main__':
    main()

