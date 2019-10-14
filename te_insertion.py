#!/usr/bin/env python

"""
1. input paired-end reads
2. align 'paired-end' reads to reference genome
3. call TE insertions using TEMP
"""

import os
import sys
import re
import gzip
import pysam
import subprocess
import argparse
import logging


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout,
    level = logging.DEBUG)
log = logging.getLogger(__name__)


# input
# fq
# out-dir

def parse_arguments():
    parser = argparse.ArgumentParser(prog='te_ins.py',
        description='')
    parser.add_argument('--fq1', nargs='+', type=str, required=True,
        help='List of FASTQ files (R1)')
    parser.add_argument('--fq2', nargs='+', type=str, required=True,
        help='List of FASTQ files (R2)')
    parser.add_argument('-o', '--out-dir', dest='out_dir', required=True,
        help='Output directory')
    parser.add_argument('--frag-length', type=int, default=0, 
        dest='frag_length', help='The length of the fragment of \
        library (insert size), 0=(samtools stats). [0]')
    parser.add_argument('--threads', type=int, default=1,
        help='Number of threads to use')
    args = parser.parse_args()
    return args


def prep_data(genome='dm3'):
    """Only support dm3, dm6, now
    Input: bwa_index, te_fa, te_anno
    """
    pass




def samtools_stats(x):
    """Extract insert size of PE BAM file
    samtools stats {in.bam} | grep 'insert size average' 

    output of command: samtools stats
    Summary Numbers
    grep ^SN | cut -f 2- | grep 'insert size average'
    """
    stat = pysam.stats('-@', '8', x)
    d = {}
    for line in stat.split('\n'):
        # check
        if not line.startswith('SN'):
            continue
        # stat
        sn, group, num = line.strip().split('\t')[0:3]
        group = group.strip(':')
        d[group] = num
    return d


def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def extact_fq(fq, out, nt=25, tail=False):
    """Extract N-nt from reads"""
    n, slen, qlen = 0, 0, 0
    # open file
    if fq.endswith('.gz'):
        f1open = gzip.open
    else:
        f1open = open

    # open file
    if out.endswith('.gz'):
        f2open = gzip.open
    else:
        f2open = open

    # suffix
    suffix = '_1'
    if tail:
        suffix = '_2'
    suffix = ''

    flag = 0
    with f1open(fq, 'rt') as fi, f2open(out, 'wt') as fo:
        for name, seq, qual in readfq(fi):
            nlen = len(seq)
            if nlen < 2 * nt:
                continue # skip
            if tail:
                s1 = seq[:nt]
                q1 = qual[:nt]
            else:
                s1 = seq[-nt:]
                q1 = qual[-nt:]
            flag += 1
            fo.write('\n'.join(['@' + str(flag) + suffix, s1, '+', q1]) + '\n')


def run_shell_cmd(cmd):
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


def strip_ext_fastq(fastq):
    return re.sub(r'\.(fastq|fq|Fastq|Fq)\.gz$','',
                    str(fastq))


def strip_merge_fastqs_prefix(fastq):
    return re.sub(r'^merge\_fastqs\_R\d\_','',str(fastq))


def bwa_aln(fq, index, out_dir, threads=8):
    basename = os.path.basename(strip_ext_fastq(fq))
    prefix = os.path.join(out_dir,
        strip_merge_fastqs_prefix(basename))
    sai = '{}.sai'.format(prefix)

    cmd = 'bwa aln -n 3 -l 100 -R 1000 -q 5 -k 2 -t {} {} {} > {}'.format(
        threads,
        index,
        fq,
        sai)
    if os.path.exists(sai):
        logging.info('file exists, skip aln: {}'.format(sai))
    else:
        run_shell_cmd(cmd)

    return sai


def bwa_pe(fq1, fq2, index, out_dir, threads=8):
    basename = os.path.basename(strip_ext_fastq(fq1))
    prefix = os.path.join(out_dir, 
        strip_merge_fastqs_prefix(basename))
    sam = '{}.sam'.format(prefix)
    bam = '{}.sorted.bam'.format(prefix)

    if os.path.exists(bam):
        logging.info('bam exists, skip alignment: ' + bam)
    else:
        temp_files = []
        sai1 = bwa_aln(fq1, index, out_dir, threads=threads)
        sai2 = bwa_aln(fq2, index, out_dir, threads=threads)
    
        # merge
        cmd = 'bwa sampe {} {} {} {} {} | samtools view -@ {} -Sub - | samtools sort -@ {} -o {} -'.format(
            index,
            sai1,
            sai2,
            fq1,
            fq2,
            threads,
            threads,
            bam)
        if os.path.exists(bam):
            logging.info('file exists, skip sampe: {}'.format(bam))
        else:
            run_shell_cmd(cmd)

    return bam


def run_ins(bam, out_dir, frag_length=500, threads=8):
    """Run TEMP for TE Insertion identification"""
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    temp_dir = '/home/wangming/work/yu_2019/projects/public_data/TE_insertion/src/TEMP'
    # te_fa = '/home/wangming/data/genome/dm3/FL10B/FL10B_and_transposon/FL10B_transposon.fa'
    # te_bed = '/home/wangming/data/genome/dm3/annotation_and_repeats/te_annotation/dm3.te.bed'
    te_fa = '/home/wangming/data/genome/dm6/dm6_transposon/dm6_transposon.fa'
    te_bed = '/home/wangming/data/genome/dm6/annotation_and_repeats/te_annotation/dm6.te.bed'
    te_ins_sh = os.path.join(temp_dir, 'scripts', 'TEMP_Insertion.sh')

    cmd = 'bash {} -i {} -s {} -o {} -r {} -t {} -f {} -x 30 -m 3 -c {}'.format(
        te_ins_sh, 
        os.path.abspath(bam),
        os.path.join(temp_dir, 'scripts'),
        out_dir,
        te_fa,
        te_bed,
        frag_length,
        threads)
    run_shell_cmd(cmd)


def run_te_insertion(fq1, fq2, out_dir, frag_length=0, threads=8):
    # index='/home/wangming/data/genome/dm3/bwa_index/dm3'
    index='/home/wangming/data/genome/dm6/bwa_index/dm6'

    ## prepare dir
    fq_dir = os.path.join(out_dir, 'data')
    aln_dir = os.path.join(out_dir, 'alignment')
    te_dir = os.path.join(out_dir, 'te_insertion')
    if not os.path.exists(fq_dir):
        os.makedirs(fq_dir)
    if not os.path.exists(aln_dir):
        os.makedirs(aln_dir)
    if not os.path.exists(te_dir):
        os.makedirs(te_dir)

    basename = os.path.basename(strip_ext_fastq(fq1))
    prefix = strip_merge_fastqs_prefix(basename)

    ## align
    bam = bwa_pe(fq1, fq2, index, aln_dir, threads=threads)

    ## insertion
    # determin the insert size
    # if frag_length > 0, use "samtools stats"
    # else: frag_length
    if frag_length == 0:
        statsDict = samtools_stats(bam)
        frag_length = statsDict.get('insert size average', 300)
    
    run_ins(bam, te_dir, frag_length, threads)
    
    
def main():
    args = parse_arguments()
    
    for i, j in zip(args.fq1, args.fq2): 
        run_te_insertion(i, j, args.out_dir, args.frag_length, args.threads)


if __name__ == '__main__':
    main()


