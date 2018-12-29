#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

## RNAseq-pipeline

+ 1. trimming

NSR library: cut 6 nt at both ends

```
$ hipipe-trim.py -i in.fq -o path_out
```

+ 2. mapping

keep both unique, multiple reads

```
$ hipipe-align.py -i in.clean.fq -o path_out -n smp -g dm3
```

**Fruitfly**

  - a. map to genome, STAR (dm3)  
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


def get_args():
    parser = argparse.ArgumentParser(prog='RNAseq-pipeline', 
                                     description='DEseq analysis')
    parser.add_argument('-a', nargs='+', required=True, metavar='control_fq',
        type=argparse.FileType('r'),
        help='FASTQ files of Control sample replicates')
    parser.add_argument('-b', nargs='+', required=True, metavar='treatment_fq',
        type=argparse.FileType('r'),
        help='FASTQ files of Treatment sample, replicates')
    parser.add_argument('-o', default=None, metavar='OUTPUT',  
        help='The directory to save results, default: [cwd]')
    parser.add_argument('-g', required=True, default='hg19', 
        metavar='GENOME', choices=['dm3', 'hg19', 'hg38', 'mm10', 'mm9'],
        help='Reference genome : dm3, hg19, hg39, mm10, default: hg19')
    parser.add_argument('-k', default=None, 
        metavar='Spike-in', choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('-x', nargs='+', metavar='align_index',
        help='Provide extra alignment index(es) for alignment, support multiple\
        indexes. eg. Transposon, tRNA, rRNA and so on.')
    parser.add_argument('--gtf', required=True, 
        help='genome annotation file in GTF format, from ensembl')
    parser.add_argument('-A', metavar='Control_NAME', default=None,
        help='Name control samples')
    parser.add_argument('-B', metavar='Control_NAME', default=None,
        help='Name control samples')
    parser.add_argument('-p', metavar='threads', type=int, default=8,
        help='Number of threads to use, default: 8')
    parser.add_argument('-s', metavar='strandness', type=int,
        default=1, choices=[0, 1, 2],
        help='strandness, 1=sens, 2=anti, 0=ignore, default:1\
        This is for featureCounts, dUTP strand-specific mode, \
        read2 is sense strand.')
    parser.add_argument('--unique-only', action='store_true',
        dest='unique_only',
        help='if specified, keep unique mapped reads only')
    parser.add_argument('--aligner', default='bowtie', 
        choices=['bowtie', 'bowtie2', 'STAR'],
        help='Choose which aligner to use. default: bowtie')
    parser.add_argument('--align-to-rRNA', dest='align_to_rRNA',
        action='store_true',
        help='if specified, align to rRNA before genome')
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
    subdirs = ['genome_mapping', 'count', 'bigWig', 'de_analysis', 'report', 
               'transposon_analysis', 'src']
    prj_dirs = [os.path.join(path, f) for f in subdirs]
    tmp = [is_path(f) for f in prj_dirs]
    # create dict
    prj_dict = {}
    for k, v in zip(subdirs, prj_dirs):
        prj_dict[k] = v
    return prj_dict



def _is_bam_indexed(bam, overwrite=False):
    """Check if *.bai file exists"""
    bai = bam + '.bai'
    if os.path.exists(bai) and overwrite is False:
        return True
    else:
        pysam.index(bam)
        if os.path.exists(bai):
            return True
        else:
            return False



def fc_run(gtf, bam, fout, strandness=0, threads=8, overwrite=False):
    """Using featureCounts to count reads for genes
    para: -M -O --fraction -T <p> -a <gtf> -s 2 <out.count> bam.files
    dUTP mode strandness
    -s 2 is sense-strand
    """
    assert os.path.exists(gtf)
    assert isinstance(bam, list)
    assert isinstance(threads, int)
    fc = which('featureCounts')
    # strandness = 2 if anti is True else 1
    if fc is None:
        raise ValueError('%10s | command not found: featureCounts' % 'failed')

    logging.info('counting: featureCounts')
    flog = fout + '.featureCounts.log'
    bam_line = ' '.join(bam)
    # check bam files
    bam_flag = [_is_bam_indexed(f) for f in bam]
    if all(bam_flag) is False:
        raise ValueError('%10s | generating bai files, wrong' % 'failed')

    # return bam names
    bam_names = [file_prefix(f)[0] for f in bam]
    with open(fout + '.bam_ids.txt', 'wt') as fo:
        for b in bam_names:
            fo.write(b + '\n')

    # return cmd        
    c1 = '%s -M -O --fraction -g gene_id -t exon -T %s -a %s -s %s \
        -o %s %s' % (fc, threads, gtf, strandness, fout, bam_line)

    with open(fout + '.cmd.txt', 'wt') as fo:
        fo.write(c1 + '\n')

    # run cmd
    print(c1)
    if os.path.exists(fout) and overwrite is False:
        logging.info('file exists: %s' % fout)
        return fout
    else:
        with open(flog, 'wt') as fo:
            p = subprocess.run(shlex.split(c1), stdout=fo, stderr=fo)
        if os.path.exists(fout):
            return fout
        else:
            logging.error('%10s | featureCounts output not correct' % 'failed')
            raise ValueError('featureCounts error: %s' % fout)



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
    path_dirs = [f for f in path_dirs if os.path.isdir(f)]
    frames = [Alignment_stat(f).stat for f in path_dirs]
    frames = [f for f in frames if isinstance(f, pd.DataFrame)]
    if len(frames) > 0:
        df = pd.concat(frames, axis=0)
        return df
    else:
        return None


def main():
    args = get_args()
    if args.o is None:
        args.o = str(pathlib.Path.cwd())

    # prep-dirs
    prj_path = prepare_project(args.o)

    # path_out
    #    |-genome_mapping
    #    |-count
    #    |-bigWig
    #    |-report
    #    |-src

    ## mapping ##

    # control
    ctl_fqs = [f.name for f in args.a]
    ctl_prefix = str_common([os.path.basename(f) for f in ctl_fqs])
    ctl_prefix = ctl_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
    if args.A is None:
        args.A = ctl_prefix
    ctl_path = os.path.join(prj_path['genome_mapping'], args.A)
    ctl_bam_files, ctl_bam_ext_files = Alignment(
        fqs=ctl_fqs, 
        path_out=ctl_path, 
        smp_name=args.A,
        genome=args.g,
        spikein=args.k, 
        index_ext=args.x,
        multi_cores=args.p,
        unique_only=args.unique_only, 
        aligner=args.aligner,
        align_to_rRNA=args.align_to_rRNA,
        path_data=args.path_data,
        overwrite=args.overwrite).run()

    # treatment
    tre_fqs = [f.name for f in args.b]
    tre_prefix = str_common([os.path.basename(f) for f in tre_fqs])
    tre_prefix = tre_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
    if args.B is None:
        args.B = tre_prefix
    tre_path = os.path.join(prj_path['genome_mapping'], args.B)
    tre_bam_files, tre_bam_ext_files = Alignment(
        fqs=tre_fqs, 
        path_out=tre_path, 
        smp_name=args.B,
        genome=args.g,
        spikein=args.k, 
        index_ext=args.x,
        multi_cores=args.p,
        unique_only=args.unique_only, 
        aligner=args.aligner,
        align_to_rRNA=args.align_to_rRNA,
        path_data=args.path_data,
        overwrite=args.overwrite).run()

    # ## create bigWig files ##
    map_bam_files = ctl_bam_files + tre_bam_files
    bw_path = prj_path['bigWig']
    for bam in map_bam_files:
        bam2bigwig(
            bam=bam, 
            genome=args.g, 
            path_out=bw_path, 
            strandness=args.s, 
            binsize=args.bin_size, 
            overwrite=args.overwrite)

    ## count ##
    count_path = prj_path['count']
    count_file = os.path.join(count_path, 'count.txt')
    # only kepp replicate samples
    map_bam_files = [f for f in map_bam_files if '_rep' in f]
    count_file = fc_run(args.gtf, map_bam_files, count_file,
        args.s, overwrite=args.overwrite)

    # ## DE analysis ##
    # using R code #
    # de_run(args.A, args.B, count_file)
    de_path = prj_path['de_analysis']
    run_deseq2 = '/home/wangming/work/wmlib/hipipe/run_deseq2.R'
    c1 = '/usr/bin/Rscript %s %s %s %s' % (run_deseq2, count_file, args.g, args.o)
    subprocess.run(shlex.split(c1), stdout=subprocess.PIPE)

    ## mapping stat ##
    map_stat_path = prj_path['report']
    map_stat_file = os.path.join(map_stat_path, 'mapping.stat')
    ctl_map = map_stat(ctl_path)
    tre_map = map_stat(tre_path)
    df_map = pd.concat([ctl_map, tre_map], axis=0).reset_index()
    df_map = df_map.sort_values(['index'])
    print(df_map)
    df_map.to_csv(map_stat_file, sep='\t', header=True, index=False)

    #########################
    ## Transposon analysis ##
    #########################
    if args.g == 'dm3': # support dm3 TE analysis only
        logging.info('## For Transposon analysis ##')
        te_path = prj_path['transposon_analysis']
        te_mapping_path = os.path.join(te_path, 'genome_mapping')
        assert is_path(te_mapping_path)

        te_ctl_path = os.path.join(ctl_path, 'extra_mapping')
        te_tre_path = os.path.join(tre_path, 'extra_mapping')
       
        ## TE mapping
        te_map_bam_files = ctl_bam_ext_files + tre_bam_ext_files
        # flatten the nested lists
        te_map_bam_files = [item for sublist in te_map_bam_files for item in sublist]

        ## make bigWig
        te_bw_path = os.path.join(te_path, 'bigWig')
        assert is_path(te_bw_path)
        # for bam in te_map_bam_files:
        #     bam2bigwig(
        #         bam=bam, 
        #         genome=args.g, 
        #         path_out=bw_path, 
        #         strandness=args.s, 
        #         binsize=args.bin_size, 
        #         overwrite=args.overwrite)

        ## count ##
        te_count_path = os.path.join(te_path, 'count')
        assert is_path(te_count_path)
        te_count_file = os.path.join(te_count_path, 'count.txt')
        # !!!!
        te_gtf = '/home/data/genome/dm3/dm3_transposon/dm3.transposon.gtf'

        # only kepp replicate samples
        te_map_bam_files = [f for f in te_map_bam_files if '_rep' in f]
        te_count_file = fc_run(te_gtf, te_map_bam_files, te_count_file,
            args.s, overwrite=args.overwrite)

        # ## DE analysis ##
        te_de_path = os.path.join(te_path, 'de_analysis')
        run_deseq2 = '/home/wangming/work/wmlib/hipipe/run_deseq2.R'
        c1 = '/usr/bin/Rscript %s %s %s %s' % (run_deseq2, te_count_file, args.g, te_path)
        subprocess.run(shlex.split(c1))

        ## mapping stat ##
        te_stat_path = os.path.join(te_path, 'report')
        # map_stat_path = prj_path['report']
        te_stat_file = os.path.join(te_stat_path, 'mapping.stat')
        te_ctl_map = map_stat(te_ctl_path)
        te_tre_map = map_stat(te_tre_path)
        df_map = pd.concat([te_ctl_map, te_tre_map], axis=0).reset_index()
        df_map = df_map.sort_values(['index'])
        print(df_map)
        df_map.to_csv(te_stat_file, sep='\t', header=True, index=False)


if __name__ == '__main__':
    main()


# EOF
