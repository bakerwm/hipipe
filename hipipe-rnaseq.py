#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
## RNAseq-pipeline

+ 1. qc (fastqc report)

  - a. fastqc + rmarkdown report (html)

+ 2. trimming (optional)

  - a. NSR library (-6,6)  

  - b. standard illumina RNAseq (0)

+ 3. mapping 

  - 0. remove MT_trRNA (optional)

  - a. map to genome (STAR)  

  - b. map to transposon (STAR, consensus sequence, optional)  

  - c. map to extra sequence (STAR, firefly, ..., sequence + annotation)

+ 4. quantification (featureCounts) 

  - a. count.txt

+ 5. report

  - a. fastqc report

  - b. trimming report  

  - c. mapping report (.txt)  

  - d. count table (.stat)


## step 2

+ 6. de analysis

DESeq2, require >= 2 replicates for each sample

  - a. DESeq2 output

  - b. report (high-quality plots *.pdf, *.rda)

  - c. bigWig (compare between samples)

  
directory structure:

mapping
├── gene
│   ├── qc
│   ├── mapping
│   ├── bigWig
│   ├── count
│   └── report
├── transposon
│   ├── qc
│   ├── mapping
│   ├── bigWig
│   ├── count
│   └── report
└── extra_genes
    ├── mapping
    ├── bigWig
    ├── count
    └── report

de_analysis (A vs B)
├── gene
│   ├── bigWig
│   ├── de_analysis (DESeq2 output)
│   ├── GO_analysis (plots, tables, clusterProfiler, Reactome)
│   └── report (DE gene list)
├── transposon
│   ├── bigWig
│   ├── de_analysis
│   └── report
└── extra_genes
    ├── bigWig
    ├── de_analysis
    └── report

"""

import os
import sys
import pathlib
import argparse
import collections
import shlex
import subprocess
import pysam
from arguments import args_init
from helper import *
from alignment import Alignment, Alignment_log, Alignment_stat
from hipipe_reporter import QC_reporter, Alignment_reporter


def get_args():
    parser = argparse.ArgumentParser(prog='RNAseq-pipeline', 
                                     description='DEseq analysis')
    parser.add_argument('-c', nargs='+', required=True, metavar='control_fq',
        help='FASTQ files of Control sample replicates')
    parser.add_argument('-t', nargs='+', required=True, metavar='treatment_fq',
        help='FASTQ files of Treatment sample, replicates')
    parser.add_argument('-o', '--path-out', default=None, dest='path_out',
        help='The directory to save results, default: [cwd]')
    parser.add_argument('-g', '--genome', required=True, default='hg19', 
        choices=['dm3', 'hg19', 'hg38', 'mm10', 'mm9'],
        help='Reference genome : dm3, hg19, hg39, mm10, default: hg19')
    parser.add_argument('-k', '--spikein', default=None, 
        choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('-x', '--index_ext', nargs='+', dest='index_ext',
        help='Provide extra alignment index(es) for alignment, support multiple\
        indexes. eg. Transposon, tRNA, rRNA and so on.')
    parser.add_argument('--gtf', default=None,
        help='genome annotation file in GTF format, from ensembl')
    parser.add_argument('-C', metavar='Control_NAME', default=None,
        help='Name control samples')
    parser.add_argument('-T', metavar='Treatment_NAME', default=None,
        help='Name control samples')
    parser.add_argument('-p', '--threads', type=int, default=8,
        help='Number of threads to use, default: 8')
    parser.add_argument('-s', metavar='strandness', type=int,
        default=1, choices=[0, 1, 2],
        help='strandness, 1=sens, 2=anti, 0=ignore, default:1\
        This is for featureCounts, dUTP strand-specific mode, \
        read2 is sense strand.')
    parser.add_argument('--aligner', default='bowtie', 
        choices=['bowtie', 'bowtie2', 'STAR'],
        help='Choose which aligner to use. default: bowtie')
    parser.add_argument('--bin-size', dest='bin_size', metavar='binsize', 
        type=int, default=50,
        help='binsize of bigWig, default: 50')
    parser.add_argument('--genome-path', dest='genome_path',
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    args = parser.parse_args()
    return args


def init_rnaseq_project(x, analysis_type=1):
    """Create directories for RNAseq project
    x str, path to directory of project
    analysis_type int, 1=mapping, 2=de_analysis
    group 1, 1=gene, 2=te, 3=other
    """
    assert isinstance(x, str)
    path_dict = {
        1: ['mapping', 'bigWig', 'count', 'report'],
        2: ['bigWig', 'de_analysis', 'report']
    }

    group_dict = {
        1: 'gene',
        2: 'transposon',
        3: 'extra_genes'
    }

    # choose A or B
    path = path_dict.get(analysis_type, None)
    if path is None:
        raise Exception('rnaseq-group, expect 1 or 2, not None')

    # create directories
    prj_path = []
    out_dict = collections.defaultdict(dict)
    for group in [1, 2, 3]:
        # group
        group_name = group_dict.get(group, None)
        if group_name is None:
            raise Exception('[1,2,3] for group supported')
        # path
        for p in path:
            path_path = os.path.join(x, group_name, p)
            prj_path.append(path_path)
            is_path(path_path)

        # save as dict    
        for n, p in zip(path, prj_path):
            out_dict[group_name][n] = p

    # return values
    return out_dict


def run_shell_cmd(cmd, log): 
    try:
        p = subprocess.Popen(cmd, shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            preexec_fn=os.setsid)
        pid = p.pid
        pgid = os.getpgid(pid)
        logging.info('run_shell_cmd: PID={}, CMD={}'.format(pid, cmd))
        ret = ''
        while True:
            line = p.stdout.readline()
            if line=='' and p.poll() is not None:
                break
            # logging.debug('PID={}: {}'.format(pid,line.strip('\n')))
            # save log
            with open(log, 'wt') as fo:
                if line:
                    print('PID={}: {}'.format(pid,line.strip('\n')))
                    fo.write(line)
                    ret += line
        p.communicate() # wait here
        if p.returncode > 0:
            raise subprocess.CalledProcessError(
                p.returncode, cmd)
        return ret.strip('\n')
    except:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
            p.terminate()
        except:
            pass
        raise Exception('Killed PID={}, PGID={}\nCMD={}\nSTDOUT={}'.format(
            pid, pgid, cmd, ret))


def run_featureCounts(gtf, bam_files, out, strandness=0, 
    threads=8, overwrite=False):
    """Count reads on each genes using featureCounts
    gtf, annotation file in GTF format
    bam, one or multiple alignemnt in BAM format
    out, count.txt
    featureCounts -M -O --fraction -T <cpu> -a <gtf> -s <0|1|2> <out.txt> bam.files
    """
    ## init
    assert os.path.exists(gtf)
    assert isinstance(bam_files, list)
    assert isinstance(threads, int)
    logging.info('running featureCounts')
    fc_exe = which('featureCounts')
    if fc_exe is None:
        raise Exception('featureCounts, not found in your $PATH')
    
    ## check BAM files (should be indexed)
    bam_check_flag = [BAM(b).index() for b in bam_files]
    if all(bam_check_flag) is False:
        raise Exception('failed to index bam files')

    ## save BAM prefix to file
    bam_prefix = [file_prefix(b)[0] for b in bam_files]
    bam_list = os.path.join(os.path.dirname(out), 'bam_ids.txt')
    with open(bam_list, 'wt') as fo:
        for b in bam_prefix:
            fo.write(b + '\n')

    ## prepare command
    bam_line = ' '.join(bam_files)
    fc_cmd = '%s -M -O --fraction -g gene_id -t exon -T %s -a %s -s %s \
        -o %s %s' % (fc_exe, threads, gtf, strandness, out, bam_line)
    fc_log = os.path.splitext(out)[0] + '.featureCounts.log'

    if os.path.exists(out) and overwrite is False:
        logging.info('file exists: %s' % out)
    else:
        run_shell_cmd(fc_cmd, fc_log)

        if not os.path.exists(out):
            raise Exception('featureCounts output not correct: %s' % out)
    return out


def mapping_gene(fq_files, smp_name, args):
    """Mapping reads to genome
    control or treatment
    args dict, the arguments of pipeline
    """
    project_path = init_rnaseq_project(args['path_out'], analysis_type=1)
    mapping_gene_path = project_path['gene']

    ## qc-report
    qc_path = os.path.join(mapping_gene_path['report'], 'qc')
    QC_reporter(fq_files, qc_path).run()

    ## update args
    args_map = args.copy() # make a copy dict
    tmp1 = args_map.pop('fq1', None) # remove 'fq1' from args
    tmp2 = args_map.pop('path_out', None) # remove 'path_out' from args
    tmp3 = args_map.pop('smp_name', None) # remove 'smp_name' form args

    map_bam_files, map_bam_ext_files = Alignment(
        fq1=fq_files, 
        path_out=mapping_gene_path['mapping'], 
        smp_name=smp_name, 
        **args_map).run()

    ## mapping-report
    map_report_path = os.path.join(mapping_gene_path['report'], 'mapping')
    map_path_list = [mapping_gene_path['mapping'], ]
    Alignment_reporter(map_path_list, map_report_path).run()
    
    ## create bigWig files
    map_bw_path = mapping_gene_path['bigWig']
    for bam in map_bam_files:
        bam2bigwig(
            bam=bam, 
            genome=args['genome'], 
            path_out=map_bw_path, 
            strandness=args['s'], 
            binsize=args['bin_size'],
            overwrite=args['overwrite'])

    ## count
    count_path = mapping_gene_path['count']
    count_file = os.path.join(count_path, 'count.txt')
    ## only kepp replicates
    map_bam_files = [f for f in map_bam_files if '_rep' in f]
    run_featureCounts(
        gtf=args['gtf'], 
        bam_files=map_bam_files, 
        out=count_file, 
        strandness=args['s'], 
        threads=args['threads'], 
        overwrite=args['overwrite'])


def mapping_te(args):
    """Mapping reads to transposon censensus
    only support dm3, currently
    args dict, the arguments of pipeline
    """
    if not args['genome'] == 'dm3':
        logging.warning('Only support TE analysis for dm3 genome')
        return None

    project_path = init_rnaseq_project(args['path_out'], group=2) # group =1, 2
    mapping_te_path = project_path['transposon']

    logging.info('Transposon analysis for TE')
    te_mapping_path = mapping_te_path['mapping']
   
    ## mapping
    # tre_bam_files, tre_bam_ext_files = Alignment(
    #     fq1=tre_fqs, path_out=tre_path, smp_name=args['T'], **args_map).run()
    te_map_bam_files = [b[0] for b in ctl_bam_ext_files] + [b[0] for b in tre_bam_ext_files]

    ## count
    count_path = mapping_te_path['count']
    count_file = os.path.join(count_path, 'count.txt')
    te_gtf = Genome(args['genome'], genome_path=args['genome_path']).te()
    # keep only replicates bam files
    te_map_bam_files = [f for f in te_map_bam_files if '_rep' in f]
    run_featureCounts(gtf=te_gtf, bam_files=te_map_bam_files, 
        out=count_file, strandness=args['s'], threads=args['threads'],
        overwrite=args['overwrite'])


def run_deseq2(control, treatment, path_out, genome, group='gene'):
    """Run DESeq2 witb count.txt matrix input
    control str, featureCounts output
    treatment str, featureCounts output
    path_out str, path to output
    group str, gene|transposon|extra_genes
    using custom R script
    """
    project_path = init_rnaseq_project(path_out, analysis_type=2)
    de_path = project_path[group] # gene|transposon|extra_genes
    
    Rscript_exe = which('Rscript')
    if Rscript_exe is None:
        raise Exception('Rscript, command not found in your $PATH')

    deseq2_script = os.path.join(sys.path[0], 'run_deseq2.R')
    if not os.path.exists(deseq2_script):
        raise Exception('R script not found: %s' % deseq2_script)

    deseq2_path = de_path['de_analysis']
    count_matrix = os.path.join(deseq2_path, 'count.txt')

    ## merge two matrix using pandas
    df1 = pd.read_csv(control, '\t', header=0, comment='#', index_col=0)
    df2 = pd.read_csv(treatment, '\t', header=0, comment='#', index_col=0)
    ## select count
    df2x = df2.iloc[:, 5:]
    df = pd.concat([df1, df2x], axis=1, join='outer').reset_index()
    df.to_csv(count_matrix, sep='\t', header=True, index=False)

    deseq2_cmd = ' '.join([Rscript_exe, deseq2_script, count_matrix, 
        genome, deseq2_path])
    deseq2_log = os.path.join(deseq2_path, 'run_DESeq2.log')
    run_shell_cmd(deseq2_cmd, deseq2_log)


def map_stat(path):
    """Stat mapping reads for each replicate
    and merge files
    input: path to control/treatment
    directory structure:
    directory
      |- control
      |    |-merged
      |    |-rep1
      |    |-rep2
      |
      |-treatment
      |    |-merged
      |    |-rep1
      |    |-rep2
    """
    path_dirs = [os.path.join(path, f) for f in os.listdir(path)]
    ## exclude extra_mapping
    path_dirs = [f for f in path_dirs if os.path.isdir(f) and 
        not os.path.basename(f).startswith('extra_mapping')]
    frames = [Alignment_stat(f).stat for f in path_dirs]

    frames = [f for f in frames if isinstance(f, pd.DataFrame)]
    if len(frames) > 0:
        df = pd.concat(frames, axis=0)
        return df
    else:
        return None


def main():
    """Main for RNAseq analysis pipeline"""
    args = args_init(vars(get_args()), trim=False, align=True)

    ## default values
    args['unique_only'] = True
    args['align_to_rRNA'] = True
    if args['gtf'] is None:
        args['gtf'] = Genome(args['genome'], 
            genome_path=args['genome_path']).gene_gtf('ensembl')

    
    ############################################################################
    ## gene analysis
    ############################################################################
    # project_path = init_rnaseq_project(args['path_out'], analysis_type=1)

    # qc stat
    # mapping
    # mapping stat
    # bigWig
    # count
    # de_analysis

    ## control, args['c']
    ctl_args = args.copy()
    ctl_fqs = ctl_args['c']
    ctl_prefix = str_common([os.path.basename(f) for f in ctl_fqs])
    ctl_prefix = ctl_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
    if args['C'] is None:
        args['C'] = ctl_prefix
    ## path to control
    ctl_args['path_out'] = os.path.join(args['path_out'], args['C'])
    # mapping_gene(ctl_fqs, args['C'], ctl_args)

    ## treatment, args['t']
    tre_args = args.copy()
    tre_fqs = tre_args['t']
    tre_prefix = str_common([os.path.basename(f) for f in tre_fqs])
    tre_prefix = tre_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
    if args['T'] is None:
        args['T'] = tre_prefix
    ## path to treatment
    tre_args['path_out'] = os.path.join(args['path_out'], args['T'])

    ## de_analysis
    ## gene, transposon, extra_genes
    de_args = args.copy()
    de_path = os.path.join(args['path_out'], args['C'] + '_vs_' + args['T'])
    is_path(de_path)

    ## gene
    count_ctl = os.path.join(ctl_args['path_out'], 'gene', 'count', 'count.txt')
    count_tre = os.path.join(tre_args['path_out'], 'gene', 'count', 'count.txt')
    gene_de_path = de_path
    run_deseq2(
        control=count_ctl,
        treatment=count_tre,
        path_out=de_path,
        genome=args['genome'],
        group='gene')

    





if __name__ == '__main__':
    main()


# EOF
