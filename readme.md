# hipipe

A group of scripts for HiSeq data, including demx, trimming, alignment etc.

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


