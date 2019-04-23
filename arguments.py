
"""
Prepare arguments for goldclip pipeline
"""
import os
import pathlib

def args_init(args=None, demx=False, trim=False, align=False, call_peak=False, bam2bw=False):
        """Inititate the arguments, assign the default values to arg
        positional arg: smp, genome
        """
        if isinstance(args, dict):
            pass
        elif args is None:
            args = {} # init dictionary
        else:
            raise Exception('unknown argument: args=%s' % args)

        args['fq1'] = args.get('fq1', None)
        args['fq2'] = args.get('fq2', None)
        args['path_out'] = args.get('path_out', str(pathlib.Path.cwd()))
        if args['path_out'] is None:
            args['path_out'] = str(pathlib.Path.cwd())

        ## optional
        genome_path = os.path.join(str(pathlib.Path.home()), 'data', 'genome')
        args['genome_path'] = args.get('genome_path', None)
        if args['genome_path'] is None:
            args['genome_path'] = genome_path
        args['overwrite'] = args.get('overwrite', False)
        args['threads'] = args.get('threads', 8)

        ## demx
        if demx:
            args['demx_type'] = args.get('demx_type', 'p7') # p7, barcode, both
            args['n_mismatch'] = args.get('n_mismatch', 0)
            args['bc_n_left'] = args.get('bc_n_left', 3)
            args['bc_n_right'] = args.get('bc_n_right', 2)
            args['bc_in_read'] = args.get('bc_in_read', 1)
            args['cut'] = args.get('cut', False)

        ## trimming
        if trim:
            args['len_min']  = args.get('len_min', 15)
            args['adapter3'] = args.get('adapter3', ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'])
            args['keep_name'] = args.get('keep_name', True)
            args['gzip'] = args.get('gzip', True) # output fastq

            args['qual_min'] = args.get('qual_min', 20)
            args['error_rate'] = args.get('error_rate', 0.1)
            args['overlap'] = args.get('overlap', 3)
            args['percent'] = args.get('percent', 80)
            args['trim_times'] = args.get('trim_times', 1)

            args['rm_untrim'] = args.get('rm_untrim', False)
            args['save_untrim'] = args.get('save_untrim', False)
            args['save_too_short'] = args.get('save_too_short', False)
            args['save_too_long'] = args.get('save_too_long', False)

            args['adapter_sliding'] = args.get('adapter_sliding', False)
            args['cut_after_trim'] = args.get('cut_after_trim', '0') # NSR
            args['rmdup'] = args.get('rmdup', False)
            args['cut_after_rmdup'] = args.get('cut_after_rmdup', '0')
            args['cut_to_length'] = args.get('cut_to_length', 0)
            args['adapter5'] = args.get('adapter5', None)

            ## PE trimming options
            args['AD3'] = args.get('AD3', ['AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'])
            args['AD5'] = args.get('AD5', None)
            args['not_trim_adapter'] = args.get('not_trim_adapter', False)
            args['keep_temp_files'] = args.get('keep_temp_files', False)

            ## deprecated
            # args['rm_dup'] = args.get('rm_dup', True) # Deprecated since v0.3
            # args['cut_before_trim'] = args.get('cut_before_trim', '0') # Deprecated since v0.3
            # args['gzipped'] = args.get('gzipped', True) # Deprecated since v0.3
            # args['double_trim'] = args.get('double_trim', False) # Deprecated since v0.3


        ## alignment
        if align:
            args['genome'] = args.get('genome', None)
            args['spikein'] = args.get('spikein', None)
            #args['index_ext'] = args.get('index_ext', None)
            args['extra_index'] = args.get('extra_index', None)
            args['unique_only'] = args.get('unique_only', True) # unique map
            args['aligner'] = args.get('aligner', 'bowtie') # bowtie alignment
            args['te_index'] = args.get('te_index', None) #
            # args['align_to_te'] = args.get('align_to_te', False) #
            args['align_by_order'] = args.get('align_by_order', True) # align reads to multiple index by order
            args['n_map'] = args.get('n_map', 0)
            args['align_to_rRNA'] = args.get('align_to_rRNA', True)
            args['repeat_masked_genome'] = args.get('repeat_masked_genome', False)
            args['merge_rep'] = args.get('merge_rep', True)
            args['small_genome'] = args.get('small_genome', False)
            args['simple_name'] = args.get('simple_name', False)

            # check-point
            if args['spikein'] == args['genome']:
                args['spikein'] = None

        ## peak-calling
        if call_peak:
            args['peak_caller'] = args.get('peak_caller', 'pyicoclip')

            ## rtstop-calling
            args['threshold'] = args.get('threshold', 1)
            args['intersect'] = args.get('intersect', 0)
            args['threads'] = args.get('threads', 8)

        ## bam2bw
        if bam2bw:
            args['filterRNAstrand'] = args.get('filterRNAstrand', None)
            args['samFlagExclude'] = args.get('samFlagExclude', None)
            args['samFlagInclude'] = args.get('samFlagInclude', None)
            args['binsize'] = args.get('binsize', 10)

        return args


def args_default(lib='rnaseq'):
    """pre-defined arguments for regular hiseq libraries"""

    ## default parameters
    args_hiseq = {
        'rnaseq' : {
            'len_min': 20,
            'adapter3': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'],
            'cut_after_trim': '7,-7',
            'rmdup': False 
        },
        'chipseq': {
            'len_min': 20,
            'adapter3': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC']
        },
        'iclip': {
            'len_min': 15,
            'adapter3': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'],
            'rmdup': True,
            'cut_after_rmdup': '9',
            'adapter_sliding': True,
            'trim_times': 4,
            'rm_untrim': True
        },
        'eclip' : {
            'len_min': 15,
            'adapter3': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'],
            'cut_after_trim': '-7',
            'rmdup': True,
            'cut_after_rmdup': '10',
            'adapter_sliding': True,
            'trim_times': 4,
            'rm_untrim': True
        },
        'clipnsr' : {
            'len_min': 15,
            'adapter3': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'],
            'cut_after_trim': '7,-7',
            'rmdup': True,
            'adapter_sliding': True,
            'trim_times': 4,
            'rm_untrim': True
        },
        'atacseq': {
            'len_min': 20,
            'adapter3': ['CTGTCTCTTATACACATCT'],
            'adapter_sliding': True,
            'trim_times': 4
        },
        'smrna': {
            'len_min': 18,
            'adapter3': ['TGGAATTCTCGGGTGCCAAGG']
        }
    }

    ## library-type
    if lib is None:
        args_lib = {}
    elif lib in args_hiseq:
        args_lib = args_hiseq[args['library_type']] # return dict
    else:
        logging.info('unknown lib : %s' % lib)
        args_lib = {}
        # raise Exception('illegal argument: --library-type %s' % args['library_type'])
    return args_lib


default_arguments = {
    'trim': {
        'adapter3' : 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
        'adapter5' : None,
        'len_min' : 15,
        'read12' : 1,
        'qual_min' : 20,
        'err_rate' : 0.1,
        'AD3' : 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
        'AD5' : None,
        'overlap' : 3,
        'threads': 8,
        'overwrite' : False,
        'rm_untrim' : False,
        'keep_name' : True,
        'adapter_sliding' : False,
        'trim_times': 1,
        'double_trim': False,
        'rm_dup' : False,
        'cut_before_trim' : 0,
        'cut_after_trim': 0,
        'trim_to_length' : 0
    },

    'align': {
        'spikein' : None,
        'index_ext' : None,
        'threads': 8,
        'unique_only' : False,
        'n_map' : 0,
        'aligner' : 'STAR',
        'align_to_rRNA' : True,
        'repeat_masked_genome' : False,
        'merge_rep' : True,
        'overwrite' : False
    },

    'peak': {
        'peak_caller' : 'pyicoclip',
        'threads' : 8,
        'overwrite' : False,
    }, 

    'rtstop': {
        'threshold' : 1, # threshold
        'intersect' : 0, # intersect, 0, 1
        'overwrite' : False
    },

    'report': {
        'group': 'homer',
        'window': 10000,
        'threads': 8
    }
}


# decrepted
class Argument(object):
    """Pre-defined arguments for goldclip analysis
    triming
    alignment
    peak-calling
    rtstop-calling
    report
    """

    def __init__(self, mode=1):
        """Get the default values for arguments
        all 5 modules
        return dict
        """
        self.mode = mode

    def trim(self):
        ## illumina TruSeq
        args_trim_basic = {
            'adapter3' : 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
            'adapter5' : None,
            'len_min' : 15,
            'read12' : 1,
            'qual_min' : 20,
            'err_rate' : 0.1,
            'AD3' : 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
            'AD5' : None,
            'overlap' : 3,
            'threads': 8,
            'overwrite' : False,
            'rm_untrim' : False,
            'keep_name' : True,
            'adapter_sliding' : False,
            'rm_dup' : False,
            'cut_before_trim' : 0,
            'trim_to_length' : 0,
        }

        ##
        cut_args = {
            1 : { 'cut_after_trim' : '7,-7' }, # NSR
            2 : { 'cut_after_trim' : '10,-7' }, # read1: N{10}-----{7-nt}
            3 : { 'cut_after_trim' : '9'} # read1: NNN{bc-4-nt}NN
        }

        args_trim_cut = cut_args[self.mode]

        args_trim = {**args_trim_basic, **args_trim_cut}

        return args_trim


    def align(self):
        # required
        # -i, -o, -g, -n, 
        args_align = {
            'spikein' : None,
            'index_ext' : None,
            'threads': 8,
            'unique_only' : False,
            'n_map' : 0,
            'aligner' : 'STAR',
            'align_to_rRNA' : True,
            'repeat_masked_genome' : False,
            'merge_rep' : True,
            'overwrite' : False
        }
        return args_align


    def peak(self):
        args_peak = {
            'peak_caller' : 'pyicoclip',
            'threads' : 8,
            'overwrite' : False,
        }
        return args_peak


    def rtstop(self):
        args_rtstop = {
            'threshold' : 1, # threshold
            'intersect' : 0, # intersect, 0, 1
            'overwrite' : False
        }
        return args_rtstop


    def all(self):
        args_all = {**self.trim(), **self.align(), **self.peak(), **self.rtstop()}
        return args_all
