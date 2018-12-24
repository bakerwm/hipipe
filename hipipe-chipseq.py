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
        choices=['bowtie', 'bowtie2', 'STAR'],
        help='Choose which aligner to use. default: bowtie2')
    parser.add_argument('-x', metavar='align_index',
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


# def chipseq_te(input_bam, ip_bam):
#     """Run ChIPseq analysis on TE consensus sequences,
#     including GFP, white, Firefly sequences"""
#     args = get_args()
#     if args.o is None:
#         args.o = str(pathlib.Path.cwd())

#     # subdirs
#     prj_path = prepare_project(args.o)
#     te_path = prj_path['transposon_analysis']

#     # make bigWig
#     # normalized by macs2 output, effective tags
#     # search IP and Input BAM files

#     # control
#     ctl_fqs = [f.name for f in args.c]
#     if args.C is None:
#         ctl_prefix = str_common([os.path.basename(f) for f in ctl_fqs])
#         ctl_prefix = ctl_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
#         args.C = ctl_prefix
#     ctl_path = os.path.join(prj_path['genome_mapping'], args.C)


#     # treatment



def bigwig2track_single(ip_bw, input_bw, fasize, n, pdf_out):
    """Create track view plots for TE consensus"""

    # R script
    r_code = os.path.splitext(pdf_out)[0] + '.R'

    c1 = """
#!/usr/bin/Rscript
library(goldclipReport)
library(rtracklayer)
library(trackViewer)
library(ggridges)
    """

    c2 = 'ip_bw    <- %s' % ip_bw
    c3 = 'input_bw <- %s' % input_bw
    c4 = 'fasize   <- %s' % fasize
    c5 = 'n        <- "%s"' % n
    c6 = 'pdf_out  <- %s' % pdf_out
    c7 = """
df_list <- chipseq_bw_parser(ip_bw, fasize, input_bw, n)

plist <- lapply(df_list, function(d){
  coverage_plot_single(d, fill.color = "orange",
                      exclude.minus.scores = TRUE)
})

plot_n_pages(plist, nrow = 5, ncol = 2, pdf_out = pdf_out)
    """

    c = '\n'.join([c1, c2, c3, c4, c5, c6, c7])

    with open(r_code, 'wt') as fo:
        fo.write(c)

    # subprocess.run(shlex.split('Rscripts %s') % r_code)






def chipseq_genome():
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
    ctl_bam_files, ext_ctl_bam_files = Alignment(
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
    tre_bam_files, ext_tre_bam_files = Alignment(
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
    ## for each replicates
    for tre_bam in tre_bam_files:
        i = tre_bam_files.index(tre_bam)
        if i >= len(ctl_bam_files):
            i = 0 # the first one
        ctl_bam = ctl_bam_files[i]
        # output directory
        tre_bam_prefix = file_prefix(tre_bam)[0]
        tre_bam_path = os.path.join(prj_path['macs2_output'], tre_bam_prefix)
        p = Macs2(ip=tre_bam, control=ctl_bam, genome=args.g, output=tre_bam_path, 
            prefix=tre_bam_prefix)
        # call peaks
        p.callpeak()
        p.bdgcmp(opt='ppois')
        p.bdgcmp(opt='FE')
        p.bdgcmp(opt='logLR')
        # annotate peaks
        p.broadpeak_annotation()


    ###################
    ## create bigWig ##
    ###################
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

    #########################
    ## transposon analysis ##
    #########################
    if ext_tre_bam_files is None or ext_ctl_bam_files is None:
        logging.info('transposon analysis skipped')
    else:
        if isinstance(ctl_bam_files, list) and isinstance(tre_bam_files, list):
            # fetch the scale
            te_path = prj_path['transposon_analysis']
            for i in ext_tre_bam_files:
                i_index = ext_tre_bam_files.index(i)
                # genome mapping BAM
                ext_tre_bam = i[0]
                tre_bam = tre_bam_files[i_index]
                if i_index >= len(ext_ctl_bam_files):
                    i_index = 0
                ext_ctl_bam = ext_ctl_bam_files[i_index][0]
                ctl_bam = ctl_bam_files[i_index]
                # fetch the normalize scale
                tre_bam_prefix = file_prefix(tre_bam)[0]
                tre_bam_path = os.path.join(prj_path['macs2_output'], tre_bam_prefix)
                p = Macs2(ip=tre_bam, control=ctl_bam, genome=args.g, output=tre_bam_path, 
                    prefix=tre_bam_prefix)
                d = p.get_effect_size() # ip_scale, ip_depth, input_scale, input_depth

                # bam to bigWig            
                te_sub_path = os.path.join(te_path, tre_bam_prefix)
                bam2bigwig2(ext_tre_bam, te_sub_path, scale=d['ip_scale'], 
                    overwrite=args.overwrite)
                bam2bigwig2(ext_ctl_bam, te_sub_path, scale=d['input_scale'], 
                    overwrite=args.overwrite)

                # save scale to file            
                s = os.path.join(te_sub_path, 'scale.lib')
                with open(s, 'wt') as fo:
                    fo.write(json.dumps(d))

                # create coverage plots
                pdf_out = os.path.join(te_sub_path, tre_bam_prefix + '.track_view.pdf')
                fasize  = 'abc.fa'
                bigwig2track_single(ext_tre_bam, ext_ctl_bam, fasize, 'P5', pdf_out)


    # return [ctl_bam_files, tre_bam_files]











def main():
    chipseq_genome()
    # chipseq_te()


if __name__ == '__main__':
    main()


# EOF
