#!/usr/bin/env python3

import os, sys
from trimmer import *

fq_in = sys.argv[1]
fq_out = sys.argv[2]
fa = sys.argv[3]
fq = sys.argv[4]
# Fastx(fq_in, fq_out, overwrite=True).trimmer(cut_left=7, cut_right=7, min_len=21, cut_to_length=20)
# Fastx(fq_in, fq_out, overwrite=True).fa2fq(fa, fq)
# Fastx(fq_in, fq_out, overwrite=True).sample(n=4)

adapter3=['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC']
path_out='bbbbbb'
args = args_init(trim=True) # save as dictionary
fq1_files = args.pop('fq1', None) # remove  'fq1'
args.pop('path_out')
args.pop('len_min')
args.pop('adapter3')
args['rm_untrim'] = True
args['cut_to_length'] = 16
# args['overwrite'] = True
Cutadapt(fq_in, adapter3, path_out, len_min=15, **args).run()