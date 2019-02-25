
"""
Prepare arguments for goldclip pipeline
"""
import pathlib

def args_init(args=None, demx=False, trim=True, align=True, call_peak=False):
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

        ## optional
        args['genome_path'] = args.get('genome_path', None)
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
            args['adapter3'] = args.get('adapter3', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC')
            args['len_min']  = args.get('len_min', 15)
            args['adapter5'] = args.get('adapter5', '')
            args['qual_min'] = args.get('qual_min', 20)
            args['error_rate'] = args.get('error_rate', 0.1)
            args['overlap'] = args.get('overlap', 3)
            args['percent'] = args.get('percent', 80)
            args['rm_untrim'] = args.get('rm_untrim', False)
            args['keep_name'] = args.get('keep_name', True)
            args['adapter_sliding'] = args.get('adapter_sliding', False)
            args['trim_times'] = args.get('trim_times', 1)
            args['double_trim'] = args.get('double_trim', False)
            args['rm_dup'] = args.get('rm_dup', True)
            args['cut_before_trim'] = args.get('cut_before_trim', '0')
            args['cut_after_trim'] = args.get('cut_after_trim', '7,-7') # NSR
            args['trim_to_length'] = args.get('trim_to_length', 0)
            args['gzipped'] = args.get('gzipped', False) # output fastq

            ## PE trimming options
            args['AD3'] = args.get('AD3', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
            args['AD5'] = args.get('AD5', None)

        ## alignment
        if align:
            args['spikein'] = args.get('spikein', None)
            args['index_ext'] = args.get('index_ext', None)
            args['unique_only'] = args.get('unique_only', True) # unique map
            args['aligner'] = args.get('aligner', 'bowtie') # bowtie alignment
            args['align-to-rRNA'] = args.get('align-to-rRNA', True)
            args['n_map'] = args.get('n_map', 0)
            args['align_to_rRNA'] = args.get('align_to_rRNA', True)
            args['repeat_masked_genome'] = args.get('repeat_masked_genome', False)
            args['merge_rep'] = args.get('merge_rep', True)

        ## peak-calling
        if call_peak:
            args['peak_caller'] = args.get('peak_caller', 'pyicoclip')

            ## rtstop-calling
            args['threshold'] = args.get('threshold', 1)
            args['intersect'] = args.get('intersect', 0)
            args['threads'] = args.get('threads', 8)

        return args


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
