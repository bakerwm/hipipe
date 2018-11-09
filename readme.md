# hipipe

A group of scripts for HiSeq data, including demx, trimming, alignment etc.


## design

### demx

1. compare two strings, mismatches

2. p7 and inline barcode, separate, or at the same time

3. PE or SE


### trim

1. `cutadapt` trim adapters

2. cut extra base at either ends

3. remove untrim reads

4. remove dup reads

5. support SE, PE reads (trimmomatic)


### align

1. determine aligner: bowtie, bowtie2, star, ...

2. parameters: unique-only, mismatch, n-mapped, 

3. index: spike-in, mt-trRNA, genome, TE, extra sequence ...

4. mappig stat


### rnaseq


### chipseq


### smRNAseq


### goldclip
























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

+ 3. view

Create bigWig for merged samples.


## ChIPseq-pipeline

+ 
























## Installation


## Usage

```
$ python <path_to_hipipe>/hipipe-trim.py
$ python <path_to_hipipe>/hipipe-align.py
$ python <path_to_hipipe>/hipipe-demx.py
```


See `foo.py` for example:

```
from trimmer import Trimmer

fq = 'demo/name_rep1.fq'
adapter3 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
path_out = 'abc'
rm_untrim = True
trim_times = 3
rm_dup = True
cut_after_trim = '5'
adapter_sliding = True
overwrite = True

a = Trimmer(fq, adapter3, path_out, 
            adapter_sliding=adapter_sliding,
            rm_dup=rm_dup, cut_after_trim=cut_after_trim,
            rm_untrim=rm_untrim, trim_times=trim_times, 
            overwrite=overwrite).run()


from alignment import Alignment

fqs = ['demo/name_rep1.fq', 'demo/name_rep2.fq']
path_out = 'abc'
smp_name = 'abc'
genome = 'dm3'
aligner = 'bowtie2'
align_to_rRNA = True
unique_only = True
overwrite = False

Alignment(fqs, path_out, smp_name, genome, 
          overwrite=overwrite, 
          align_to_rRNA=align_to_rRNA,
          aligner=aligner, 
          unique_only=unique_only).run()
```


