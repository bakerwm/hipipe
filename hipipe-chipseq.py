#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

## ChIPseq-pipeline

+ 1. trimming

```
$ hipipe-trim.py -i in.fq -o path_out
```

+ 2. align IP/Input to genome

Alignment mode:
0, unique_only
1. unique_only
2. randomely assigned multimapper
4. CSEM allocated multimapper

Using Bowtie2, keep both unique, multiple reads

--very-sensitive-local

```
$ hipipe-align.py -i in.clean.fq -o path_out -n smp -g dm3
```

**Fruitfly**

  - a. map to genome, Bowtie2 (dm3)  
  - b. map to transposon consensus
  - c. featureCount, quantification  
  - d. DESeq2 analysis (from matrix)
  - e. make plots: scatter plot, MA plot, volcano plot  
  - f. scatter plot: RPM
  - g. summary report:
       trimming report: raw, too-short, clean
       mapping report: unique, multiple, MT_trRNA, unmap

+ 3. DEanalysis

DESeq2 analysis

# 4. vier

Create bigWig for merged samples.

"""

import os
import sys
import pathlib
import argparse
import shlex
import subprocess
import pysam
from helper import *
from alignment import Alignment, Alignment_log, Alignment_stat
from hipipe_macs2 import Macs2, Venv


def get_args():
    parser = argparse.ArgumentParser(prog='ChIPseq-pipeline', 
                                     description='ChIPseq analysis')
    parser.add_argument('-c', nargs='+', required=True, metavar='control_fq',
        type=argparse.FileType('r'),
        help='FASTQ files of Control sample replicates')
    parser.add_argument('-t', nargs='+', required=True, metavar='treatment_fq',
        type=argparse.FileType('r'),
        help='FASTQ files of Treatment sample, replicates')
    parser.add_argument('-o', default=None, metavar='OUTPUT',  
        help='The directory to save results, default: [cwd]')
    parser.add_argument('-g', required=True, default='hg19', 
        metavar='GENOME', choices=['dm3', 'hg19', 'hg38', 'mm10', 'mm9'],
        help='Reference genome : dm3, hg19, hg39, mm10, default: hg19')
    # parser.add_argument('-k', default=None, 
    #     metavar='Spike-in', choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
    #     help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('-C', metavar='Control', default=None,
        help='Name control samples')
    parser.add_argument('-T', metavar='Treatment', default=None,
        help='Name treatment samples')
    parser.add_argument('-p', dest='threads', metavar='threads', type=int, 
        default=8, help='Number of threads to use, default: 8')
    parser.add_argument('--aligner', default='bowtie2', 
        choices=['bowtie', 'bowtie2', 'star'],
        help='Choose which aligner to use. default: bowtie2')
    parser.add_argument('-x', nargs='+', metavar='align_index',
        help='A fasta file used to do post-genomic mapping and analysis. For \
        example, -x FL10B.fa, after mapping reads to genome, reads are mapped to \
        FL10B sequence.')
    parser.add_argument('--bin-size', dest='bin_size', metavar='binsize', 
        type=int, default=50,
        help='binsize of bigWig, default: 50')
    parser.add_argument('--path_data', 
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    args = parser.parse_args()
    return args


def prepare_project(path):
    """Prepare subdirs for project"""
    # assert os.path.isdir(path)
    assert isinstance(path, str)
    subdirs = ['genome_mapping', 'bigWig', 'macs2_output', 
        'transposon_analysis', 'src']
    prj_dirs = [os.path.join(path, f) for f in subdirs]
    tmp = [is_path(f) for f in prj_dirs]
    # create dict
    prj_dict = {}
    for k, v in zip(subdirs, prj_dirs):
        prj_dict[k] = v
    return prj_dict


def main():
    args = get_args()
    if args.o is None:
        args.o = str(pathlib.Path.cwd())

    # prep-dirs
    # subdirs = ['genome_mapping', 'count', 'bigWig', 'report', 'src']
    prj_path = prepare_project(args.o)

    # path_out
    #    |-genome_mapping
    #    |-bigWig
    #    |-macs2_output
    #    |-transposon_analysis
    #    |-src

    ## Alignment ##

    # control
    ctl_fqs = [f.name for f in args.c]
    if args.C is None:
        ctl_prefix = str_common([os.path.basename(f) for f in ctl_fqs])
        ctl_prefix = ctl_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
        args.C = ctl_prefix
    ctl_path = os.path.join(prj_path['genome_mapping'], args.C)
    ctl_bam_files = Alignment(
        ctl_fqs, ctl_path, 
        smp_name=args.C,
        genome=args.g,
        spikein=None,
        index_ext=args.x,
        threads=args.threads,
        unique_only=True,
        n_map=1,
        aligner=args.aligner,
        align_to_rRNA=True,
        merge_rep=False,
        path_data=args.path_data,
        overwrite=args.overwrite).run()

    # treatment
    tre_fqs = [f.name for f in args.t]
    if args.T is None:
        tre_prefix = str_common([os.path.basename(f) for f in tre_fqs])
        tre_prefix = tre_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
        args.T = tre_prefix
    tre_path = os.path.join(prj_path['genome_mapping'], args.T)
    tre_bam_files = Alignment(
        tre_fqs, tre_path, 
        smp_name=args.T,
        genome=args.g,
        spikein=None,
        index_ext=args.x,
        threads=args.threads,
        unique_only=True,
        n_map=1,
        aligner=args.aligner,
        align_to_rRNA=True,
        merge_rep=False,
        path_data=args.path_data,
        overwrite=args.overwrite).run()


    # ## create bigWig files ##
    # map_bam_files = ctl_bam_files + tre_bam_files
    # bw_path = prj_path['bigWig']
    # for bam in map_bam_files:
    #     bam2bigwig(
    #         bam=bam, 
    #         genome=args.g, 
    #         path_out=bw_path, 
    #         strandness=args.s, 
    #         binsize=args.bin_size, 
    #         overwrite=args.overwrite)

    # ## mapping stat ##
    # map_stat_path = prj_path['report']
    # map_stat_file = os.path.join(map_stat_path, 'mapping.stat')
    # ctl_map = map_stat(ctl_path)
    # tre_map = map_stat(tre_path)
    # df_map = pd.concat([ctl_map, tre_map], axis=0).reset_index()
    # df_map = df_map.sort_values(['index'])
    # df_map.to_csv(map_stat_file, sep='\t', header=True, index=False)

    # report
    # mapping report
    # de analysis report
    # plots
    # gene lists
    # run DE analysis

    ################
    ## call peaks ##
    ################
    # p = Macs2(ip=tre, control=ctl, genome='dm3', output=out, prefix=None)
    # p.callpeak()
    # p.bdgcmp()
    # d = p.get_effect_size()
    # print(d)
    # # p.bdgcmp(opt='ppois')
    ## for each replicates
    for tre_bam in tre_bam_files:
        i = tre_bam_files.index(tre_bam)
        if i >= len(ctl_bam_files):
            i = 0 # the first one
        ctl_bam = ctl_bam_files[i]
        # output directory
        tre_bam_prefix = file_prefix(tre_bam)[0]
        tre_bam_path = os.path.join()
        Macs2(ip=tre_bam, control=ctl_bam, genome=args.g, output=out, prefix=prefix)


    print('control')
    print(ctl_bam_files)
    print('treatment')
    print(tre_bam_files)









    # #########################
    # ## Transposon analysis ##
    # #########################
    # logging.info('## For Transposon analysis ##')
    # te_path = prj_path['transposon_analysis']
    # ## mapping ##
    # te_mapping_path = os.path.join(te_path, 'genome_mapping')
    # assert is_path(te_mapping_path)
    # # control
    # ctl_fqs = [f.name for f in args.a]
    # ctl_prefix = str_common([os.path.basename(f) for f in ctl_fqs])
    # ctl_prefix = ctl_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
    # if args.A is None:
    #     args.A = ctl_prefix
    # te_ctl_path = os.path.join(te_mapping_path, args.A)
    # te_ctl_bam_files = Alignment(
    #     fqs=ctl_fqs, 
    #     path_out=te_ctl_path, 
    #     smp_name=args.A,
    #     genome=args.g,
    #     spikein=args.k,
    #     index_ext=args.x,
    #     multi_cores=args.p,
    #     unique_only=args.unique_only, 
    #     aligner=args.aligner,
    #     align_to_rRNA=args.align_to_rRNA,
    #     path_data=args.path_data,
    #     overwrite=args.overwrite).run()

    # # treatment
    # tre_fqs = [f.name for f in args.b]
    # tre_prefix = str_common([os.path.basename(f) for f in tre_fqs])
    # tre_prefix = tre_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
    # if args.B is None:
    #     args.B = tre_prefix
    # te_tre_path = os.path.join(te_mapping_path, args.B)
    # te_tre_bam_files = Alignment(
    #     fqs=tre_fqs, 
    #     path_out=te_tre_path, 
    #     smp_name=args.B,
    #     genome=args.g,
    #     spikein=args.k,
    #     index_ext=args.x,
    #     multi_cores=args.p,
    #     unique_only=args.unique_only, 
    #     aligner=args.aligner,
    #     align_to_rRNA=args.align_to_rRNA,
    #     path_data=args.path_data,
    #     overwrite=args.overwrite).run()


    # # ## create bigWig files ##
    # te_map_bam_files = te_ctl_bam_files + te_tre_bam_files
    # te_bw_path = os.path.join(te_path, 'bigWig')
    # assert is_path(te_bw_path)
    # # for bam in te_map_bam_files:
    # #     bam2bigwig(
    # #         bam=bam, 
    # #         genome=args.g, 
    # #         path_out=bw_path, 
    # #         strandness=args.s, 
    # #         binsize=args.bin_size, 
    # #         overwrite=args.overwrite)

    # # ## count ##
    # te_count_path = os.path.join(te_path, 'count')
    # assert is_path(te_count_path)
    # te_count_file = os.path.join(te_count_path, 'count.txt')
    # # !!!!
    # te_gtf = '/home/data/genome/dm3/dm3_transposon/dm3.transposon.gtf'
    # # only kepp replicate samples
    # te_map_bam_files = [f for f in te_map_bam_files if '_rep' in f]
    # te_count_file = fc_run(te_gtf, te_map_bam_files, te_count_file,
    #     args.s, overwrite=args.overwrite)

    # # ## DE analysis ##
    # # using R code #
    # # de_run(args.A, args.B, count_file)
    # te_de_path = os.path.join(te_path, 'de_analysis')
    # run_deseq2 = '/home/wangming/work/wmlib/hipipe/run_deseq2.R'
    # c1 = '/usr/bin/Rscript %s %s %s' % (run_deseq2, te_count_file, te_path)
    # subprocess.run(shlex.split(c1))


    # ## mapping stat ##
    # te_stat_path = os.path.join(te_path, 'report')
    # # map_stat_path = prj_path['report']
    # te_stat_file = os.path.join(te_stat_path, 'mapping.stat')
    # te_ctl_map = map_stat(te_ctl_path)
    # te_tre_map = map_stat(te_tre_path)
    # df_map = pd.concat([te_ctl_map, te_tre_map], axis=0).reset_index()
    # df_map = df_map.sort_values(['index'])
    # print(df_map)
    # df_map.to_csv(te_stat_file, sep='\t', header=True, index=False)



if __name__ == '__main__':
    main()


# EOF
