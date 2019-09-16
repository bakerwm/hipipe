#!/us/bin/env python3

"""Download sequencing data (.sra) from SRA or ENA 

NCBI/SRA: 
url: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR915/SRR9158298/SRR9158298.sra
aspera: anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR949/SRR949627/SRR949627.sra


EBI/ENA: 
url: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR008/ERR008901/ERR008901_1.fastq.gz
aspera: era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR008/ERR008901/ERR008901_1.fastq.gz  

exceptions:
6-digits: SRR915829
SRR915/SRR915829/SRR915829.fastq.gz

7-digits: SRR9158298
SRR915/008/SRR9158298/SRR9158298.fastq.gz


## FTP download
## Aspera download

http://bioinfostar.com/2017/12/23/How-to-download-SRA-data-zh_CN/
http://bioinfostar.com/2017/12/24/How-to-download-SRA-data-en/
https://www.internationalgenome.org/faq/how-download-ena-files-using-aspera/


## aspera connect download

ascp -i bin/aspera/etc/asperaweb_id_dsa.openssh -Tr -Q -l 100M -L- fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR008/ERR008901/ERR008901_1.fastq.gz ./

-i PRIVATE-KEY-FILE             Private-key file name (id_rsa)
-T                              Disable encryption
-l MAX-RATE                     Max transfer rate
-L LOCAL-LOG-DIR                Local logging directory path

"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2019-09-16"
__version__ = "0.1"

import os
import sys
import re
import argparse
import subprocess
import ftplib
import logging

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)

log = logging.getLogger(__name__)


def get_args():
    """Mapping SE read or one of PE reads to reference genome
    using bowtie, STAR, ... (universal)
    """
    parser = argparse.ArgumentParser(prog='sra_download.py', 
                                     description='download SRA file')
    parser.add_argument('-i', '--input', nargs='+', dest='input', required=True,
        help='SRR123456 or a file contain list of SRR id ')
    parser.add_argument('-o', '--out', dest='out', default=None,
        help='Output directory to save the files, default: [pwd]')
    parser.add_argument('-s', '--source', dest='source', default='ncbi',
        choices=['ncbi', 'ebi'],
        help='Choose the source to download files, ncbi or ebi, default: [ncbi]')
    parser.add_argument('-t', dest='transport', default='fasp',
        choices=['fasp', 'ftp', 'both'],
        help='Transport: one of: fasp; ftp; both. default: [fasp]')
    parser.add_argument('-l', '--max-rate', dest='max_rate', default='24m',
        help='Max transfer rate, G/g, M/m, K/k, or bytes, default: [24m]')
    parser.add_argument('--log-level', default='INFO', dest='log_level',
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')

    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)

    return args


def run_shell_cmd(cmd):
    """This command is from 'ENCODE-DCC/atac-seq-pipeline'
    https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_common.py
    """
    p = subprocess.Popen(['/bin/bash','-o','pipefail'], # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid) # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(pid, pgid, rc,
        stderr.strip(), stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')


def id_parser(l):
    """Check a list of records
    can be mixed: SRR, files
    """
    out = []
    for n in l:
        out += id_parser_n(n)
    return out


def id_parser_n(x):
    """Check x is SRR123456 or a file with SRR123456
    only first column of file parsed
    """
    out = []
    if os.path.exists(x):
        with open(x, 'rt') as fi:
            for line in fi:
                s = line.strip().split()[0]
                if s.startswith('#') or s == '':
                    continue
                if id_validator(s):
                    out.append(s)
                else:
                    log.error('unknown format: {}'.format(s))
    else:
        if id_validator(x):
            out = [x]
        else:
            log.error('unknown format: {}'.format(x))

    return out


def id_validator(x, url_check=False):
    """Check the id
    SRR123456
    SRR1234567
    ERR|DRR|SRR
    """
    # id format
    p = re.compile('^[SED]RR[0-9]{6,7}$')
    flag = re.match(p, x)

    if url_check:
        # check url exists
        # ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR008/ERR008901/ERR008901_1.fastq.gz
        url = id_url(x, fasp=False, ENA=True)
        flag2 = ftp_dir_checker(url)
        # update flag
        flag = flag and flag2

    return flag


def id_formatter(x, ENA=False):
    """Convert SRR id to /SRR008/SRR008001/
    SRA: /SRR/SRR915/SRR9158298/SRR9158298.sra
    
    ENA: 
        /SRR123/SRR123456/SRR123456.fastq.gz
        /SRR915/008/SRR9158298/SRR9158298.fastq.gz
    se|pe: only return se
    """
    p1 = x[:3] # SRR
    p2 = x[:6] # SRR915
    px = '{:03d}'.format(int(x[-1])) # 008

    p2 = p2 + px if ENA and len(x) > 9 else p2
    fn = x + '.fastq.gz' if ENA else x + '.sra'
    id_path = r'/'.join(['/', p1, p2, x, fn])

    return id_path


def id_url(x, fasp=1, ENA=0):
    """ENA, 0:NCBI, 1:ENA
    server:
        1  EBI: era-fasp@fasp.sra.ebi.ac.uk:
        0  NCBI: anonftp@ftp.ncbi.nlm.nih.gov:
    root: 
        1  EBI: /vol1/fastq
        0  NCBI: /sra/sra-instant/reads/ByRun/sra
    """
    # server
    if fasp:
        server = 'era-fasp@fasp.sra.ebi.ac.uk:' if ENA else 'anonftp@ftp.ncbi.nlm.nih.gov:'
    else:
        server = 'ftp://ftp.sra.ebi.ac.uk' if ENA else 'ftp://ftp.ncbi.nlm.nih.gov'

    # prefix
    path_prefix = '/vol1/fastq' if ENA else '/sra/sra-instant/reads/ByRun/sra'

    # suffix
    path_suffix = id_formatter(x, ENA=ENA)

    return server + path_prefix + path_suffix


def ftp_dir_checker(url):
    """List files in ftp directory:
    support ncbi|ebi, ftp
    """
    # remove prefix
    url2 = url.replace(r'ftp://', '')

    if url2.startswith('ftp.'):
        # parse the url
        # server, path, file
        server, remote_path = url2.split('/', 1)
        remote_dir = os.path.dirname(remote_path) # 
        remote_filename = os.path.basename(remote_path) # 

        # login ftp
        # check dir exists, return files in path_dir
        try:
            ftp = ftplib.FTP(server)
            ftp.login()
            ftp.cwd(remote_dir) # switch to directory
            flag = 1
            ftp.quit()
        except:
            flag = 0
    else:
        logging.error('url is not FTP: ' + url)
        flag = 0
            
    return flag


def ftp_dir_list(url):
    """List files in ftp url
    support ncbi|ebi, ftp
    """
    # remove prefix
    url2 = url.replace(r'ftp://', '')

    # parse the url
    # server, path, file
    server, remote_path = url2.split('/', 1)
    remote_dir = os.path.dirname(remote_path) # 
    remote_filename = os.path.basename(remote_path)  # 

    # login ftp
    # check dir exists, return files in path_dir
    remote_filenames = []
    try:
        ftp = ftplib.FTP(server)
        ftp.login()
        ftp.cwd(remote_dir) # switch to directory
        ftp.retrlines('NLST', remote_filenames.append) # list of file names
        ftp.quit()
    except:
        raise Exception('failed to list files: ' + url)


    return [remote_filenames, server, remote_dir]


def ftp_download(url, outdir, max_rate='24m'):
    """Download files using FTP
    wget --limit-rate=5m -O file.txt url

    1. list files
    2. download

    SRA: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR915/SRR9158298/SRR9158298.sra
    ENA: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR008/ERR008901/ERR008901_1.fastq.gz

    """
    remote_filenames, server, remote_dir = ftp_dir_list(url)

    if not os.path.basename(url) in remote_filenames:
        log.error('file not exists in remote server: ' + url)
        remote_filenames = [] # blank
    else:
        try:
            # connect to FTP server
            ftp = ftplib.FTP(server)
            ftp.login()
            ftp.cwd(remote_dir) # switch to directory
        except:
            raise Exception('failed to connect FTP server: ' + server)

        count = 0
        for f in remote_filenames:
            count += 1
            log.info('[{}/{}] - {}'.format(count, len(remote_filenames), f))
            local_file = os.path.join(outdir, f)
            if os.path.exists(local_file):
                log.info('file exists: ' + local_file)
            else:
                handle = open(local_file, 'wb')
                ftp.retrbinary('RETR {}'.format(f), handle.write, 8*1024)

        ftp.quit()


def aspera_bin():
    """Return the ascp, *openssh 
    default: 
    ascp: ~/.aspera/connect/bin/ascp
    *openssh: ~/.aspera/connect/etc/asperaweb_id_dsa.openssh
    """
    home = os.path.expanduser('~')
    ascp = os.path.join(home, '.aspera', 'connect', 'bin', 'ascp')
    openssh = os.path.join(home, '.aspera', 'connect', 'etc', 'asperaweb_id_dsa.openssh')
    flag = [ascp, openssh]

    if not os.path.exists(ascp) or not os.path.exists(openssh):
        flag = 0
        log.error('aspera installation not correct, expect: ~/.aspera/connect')

    return flag


def aspera_download(url, outdir, max_rate='24m'):
    """Download files by ascp   
    example:
    ascp -i <asperaweb_id_dsa.openssh with path> -k1 -Tr l100m 
    anonftp@ftp.ncbi.nlm.nih.gov:/<files to transfer> <local destination>

    -i <asperaweb_id_dsa.openssh with path> = fully qualified path & file name where
    this public key file is located. This file is part of Aspera Connect distribution 
    and is usually located in the ‘etc’ subdirectory.

    -T                              Disable encryption
    -k RESUME-LEVEL                 Resume criterion: 0,3,2,1
    -l MAX-RATE                     Max transfer rate
                                    RATE: G/g(gig),M/m(meg),K/k(kilo)
    -@ RANGE-LOW:RANGE-HIGH         Transfer only a range of bytes within file

    """
    aspera_out = aspera_bin()
    flag = 0
    if aspera_out:
        ascp, ascp_openssh = aspera_bin()
        # parameter
        # -T -k 1 -l 5m
        key = '-i ' + ascp_openssh
        para = '-T -k 1 -l ' + max_rate
        # run
        cmd = ' '.join([ascp, key, para, url, outdir])
        run_shell_cmd(cmd)
        # output file
        local_file = os.path.join(outdir, os.path.basename(url))
        if os.path.exists(local_file):
            flag = 1

    return flag


def file_downloader(x, outdir, max_rate='24m', fasp=True, ENA=False):
    """Download sra files from EBI|ENA
    using fasp or FTP 
    """
    downloader = aspera_download if fasp else ftp_download
    url = id_url(x, fasp, ENA)

    # parse remote files
    url_ftp = id_url(x, fasp=False, ENA=ENA)

    remote_filenames = ftp_dir_list(url_ftp)[0]

    for f in remote_filenames:
        url_f = os.path.join(os.path.dirname(url), f)
        local_file = os.path.join(outdir, f)
        if os.path.exists(local_file):
            log.info('[{}] local file exists, downloading skipped ...'.format(f))
        else:
            # downloader(url_f, outdir, max_rate)
            print(url_f)


def main():
    args = get_args()

    # outdir
    if not args.out:
        args.out = os.getcwd()

    if not os.path.exists(args.out):
        os.makedirs(args.out)

    # ids
    ids = id_parser(args.input)

    # check aspera
    fasp = True if args.transport == 'fasp' else False
    ena = args.source == 'ebi'

    # counter
    count = 0
    for i in ids:
        count += 1
        print('[{}/{}] {}'.format(count, len(ids), i))
        # file_downloader(i, args.out, args.max_rate, fasp, ena)


if __name__ == '__main__':
    main()


# 2019-09-15
# 
# old functions:

# def sra_download(x, outdir, max_rate='5m', fasp=True, ENA=False):
#     """Download file using FTP
#     try: aspera -> ftp
#     try: SRA > ENA
#     """
#     # get url
#     url = id_url(x, fasp, ENA)

#     # choose the source: SRA or ENA
#     # if ENA, check SE or PE
#     if ENA:
#         url_ftp = id_url(x, fasp=False, ENA=True)
#         remote_filenames = ftp_dir_list(url_ftp)[0]

#         # SE|PE
#         if len(remote_filenames) == 2: # PE reads
#             url_1 = url.replace('.fastq.gz', '_1.fastq.gz')
#             url_2 = url.replace('.fastq.gz', '_2.fastq.gz')
#             url = [url_1, url_2]

#     # fasp, ftp
#     if fasp:
#         if isinstance(url, list):
#             aspera_download(url[0], outdir, max_rate)
#             aspera_download(url[1], outdir, max_rate)
#         else:
#             aspera_download(url, outdir, max_rate)
#     else:
#         if isinstance(url, list):
#             ftp_download(url[0], outdir, max_rate)
#             ftp_download(url[1], outdir, max_rate)
#         else:
#             ftp_download(url, outdir, max_rate)


# def ncbi_download(x, outdir, max_rate='5m', fasp=True):
#     """Download sra files from NCBI|SRA, 
#     using fasp or FTP 
#     """
#     downloader = aspera_download if fasp else ftp_download

#     url = id_url(x, fasp, ENA=False)
#     local_file = os.path.join(outdir, x + '.sra')

#     if os.path.exists(local_file):
#         log.info('[{}] local file exists, downloading skipped ...'.format(x + '.sra'))
#     else:
#         downloader(url, outdir, max_rate)
    

# def ebi_download(x, outdir, max_rate='5m', fasp=True):
#     """Download sra files from EBI|ENA
#     using fasp or FTP 
#     """
#     downloader = aspera_download if fasp else ftp_download

#     url = id_url(x, fasp, ENA=True)
#     remote_filenames = ftp_dir_list(url)[0]

#     for f in remote_filenames:
#         url_f = os.path.join(os.path.dirname(url), f)
#         local_file = os.path.join(outdir, f)
#         if os.path.exists(local_file):
#             log.info('[{}] local file exists, downloading skipped ...'.format(f))
#         else:
#             downloader(url_f, outdir, max_rate)
