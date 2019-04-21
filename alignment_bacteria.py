#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check if contaminations exists in input files

1. bowtie method (not recommended)

2. kraken2/centrifuge method


## NCBI Taxonomy classification


## standard output of kraken2
Like Kraken 1, Kraken 2 offers two formats of sample-wide results. 
Kraken 2's standard sample report format is tab-delimited with one line 
per taxon. The fields of the output, from left-to-right, are as follows:

1. Percentage of fragments covered by the clade rooted at this taxon
2. Number of fragments covered by the clade rooted at this taxon
3. Number of fragments assigned directly to this taxon
4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
   (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
   Taxa that are not at any of these 10 ranks have a rank code that is
   formed by using the rank code of the closest ancestor rank with
   a number indicating the distance from that rank.  E.g., "G2" is a
   rank code indicating a taxon is between genus and species and the
   grandparent taxon is at the genus rank.
5. NCBI taxonomic ID number
6. Indented scientific name

source: https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#standard-kraken-output-format

"""

import os
import sys
import subprocess
import logging
import pandas as pd


logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

class Kraken2(object):
    """Run taxonomic classification using program kraken2. (version 2.0.7-beta)
    mission:

    1. Construct database
    download NCBI data (human, bacteria, archaea, and virus)
    require >130 GB disk space to process the database,
    final size is ~100 GB.

    2. Run program

    3. Return report
    """

    def __init__(self, input, outdir, kraken2_db, kraken2_exe=None, 
        threads=16, overwrite=False):
        """retrieve the input fastq file, outdir director,
        check kraken2 from your PATH
        specify db_path
        """
        ## check arguments
        assert isinstance(input, str)
        assert isinstance(kraken2_db, str)
        assert isinstance(threads, int)

        if not os.path.exists(outdir):
            try:
                os.makedirs(outdir)
            except IOError:
                    logging.error('failed to create directories: %s' % outdir)

        if kraken2_exe is None:
            kraken2_exe = self.which('kraken2')
        if kraken2_exe is None:
            self.install_kraken2()
            sys.exit('kraken2 not found')

        ## output files
        if input.endswith('.gz'):
            args = '--gzip-compressed'
            tmp = os.path.splitext(os.path.basename(input))[0]
            prefix = os.path.splitext(tmp)[0]
        elif input.endswith('.bz2'):
            args = '--bizp2-compressed'
            tmp = os.path.splitext(os.path.basename(input))[0]
            prefix = os.path.splitext(tmp)[0]
        else:
            args = ''
            prefix = os.path.splitext(os.path.basename(input))[0]

        kraken2_log = os.path.join(outdir, prefix + '.kraken2.log')
        kraken2_output = os.path.join(outdir, prefix + '.kraken2.out')
        kraken2_report = os.path.join(outdir, prefix + '.kraken2.report')

        ## arguments
        self.input = input
        self.outdir = outdir
        self.kraken2_exe = kraken2_exe
        self.kraken2_db = kraken2_db
        self.threads = threads
        self.overwrite = overwrite
        self.args = args
        self.prefix = prefix
        self.kraken2_log = kraken2_log
        self.kraken2_output = kraken2_output
        self.kraken2_report = kraken2_report


    def which(self, program):
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file
        return None


    def install_kraken2(self):
        """Install kraken2 
        instructions
        """
        install_notes = """
        ----------------------------------------------------------------------
        Here is brief installation instruction for Kraken2,
        See official website for more details:
        https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual
        ----------------------------------------------------------------------

        1. Requirements

        100 GB disk space, 29 GB RAM

        2. Download source code from github: 

        $ cd ~/biosoft/
        $ git clone https://github.com/DerrickWood/kraken2.git
        $ cd kraken2
        $ ./install_kraken2.sh <path_to_kraken2>
        
        # <path_to_kraken2> if the directory where you want to install
        # I choose "./bin"
        # copy main Kraken2 scripts to a directory in your PATH
        
        $ cp <path_to_kraken2>/kraken2{,-build,-inspect} <PATH_directory>

        $ which kraken2

        your installation is correct, if above command returns the path of kraken2.

        """

        logging.info('Please go throught the following instructions \n%s' % install_notes)


    def install_std_db(self):
        """Install Kraken2 database

        1. Standard Database

        kraken2-build --standard --db <db_name>
        kraken2-build --standard --threads 16 --db <db_name>

        2. Custom Database

            2.1 install a taxonomy

            $ kraken2-build --download-taxonomy --db $DBNAME

        """
        db_notes = """
        Standard database is fine for most research
        (Kraken2 also support custom database, see official website for details)

        $ kraken2-build --standard --threads 16 --db <db_name>

        The program requires ~140 GB of disk space during creation, and hours to download
        NCBI datasets depends on your network.

        !!!                                                                !!!
        !!! The program stopped while processing the third library         !!!
        !!! you need to type "y"                                           !!!
        !!! to confirm the file replacement, x to "assembly_summary.txt"   !!!
        !!!                                                                !!!
        mv: replace 'assembly_summary.txt', overriding mode 0444 (r--r--r--)? y

        The full log looks like this:

        Step 1/2: Performing rsync file transfer of requested files
        Rsync file transfer complete.
        Step 2/2: Assigning taxonomic IDs to sequences
        Processed 297 projects (449 sequences, 755.01 Mbp)... done.
        All files processed, cleaning up extra sequence files... done, library complete.
        Masking low-complexity regions of downloaded library... done.
        Step 1/2: Performing rsync file transfer of requested files
        Rsync file transfer complete.
        Step 2/2: Assigning taxonomic IDs to sequences
        Processed 15151 projects (32846 sequences, 60.75 Gbp)... done.
        All files processed, cleaning up extra sequence files... done, library complete.
        Masking low-complexity regions of downloaded library... done.
        Step 1/2: Performing rsync file transfer of requested files
        Rsync file transfer complete.
        Step 2/2: Assigning taxonomic IDs to sequences
        Processed 8609 projects (11060 sequences, 278.49 Mbp)... done.
        All files processed, cleaning up extra sequence files... done, library complete.
        Masking low-complexity regions of downloaded library... done.
        mv: replace 'assembly_summary.txt', overriding mode 0444 (r--r--r--)? y
        Step 1/2: Performing rsync file transfer of requested files
        Rsync file transfer complete.
        Step 2/2: Assigning taxonomic IDs to sequences
        Processed 1 project (594 sequences, 3.26 Gbp)... done.
        All files processed, cleaning up extra sequence files... done, library complete.
        Downloading UniVec_Core data from server... done.
        Adding taxonomy ID of 28384 to all sequences... done.
        Masking low-complexity regions of downloaded library... done.
        Creating sequence ID to taxonomy ID map (step 1)...
        Sequence ID to taxonomy ID map complete. [0.302s]
        Estimating required capacity (step 2)...
        Estimated hash table requirement: 38219976408 bytes
        Capacity estimation complete. [9m22.571s]
        Building database files (step 3)...
        Taxonomy parsed and converted.
        CHT created with 15 bits reserved for taxid.
        Completed processing of 48085 sequences, 65038298666 bp
        Writing data to disk...  complete.
        Database files completed. [50m55.800s]
        Database construction complete. [Total: 1h0m19.486s]
        """

        logging.info('Install standard Kraken2 database: \n %s' % db_notes)


    def db_checker(self):
        """Check whether the path is a valid kraken2 database
        
        $ kraken2-inspect --db <db_path> --skip-counts

        Standard database

        # Database options: nucleotide db, k = 35, l = 31
        # Spaced mask = 11111111111111111111111111111111110011001100110011001100110011
        # Toggle mask = 1110001101111110001010001100010000100111000110110101101000101101
        # Total taxonomy nodes: 23303
        # Table size: 6689453365
        # Table capacity: 9554994102
        # Min clear hash value = 0


        minikraken2_v1_8GB database

        # Database options: nucleotide db, k = 35, l = 31
        # Spaced mask = 11111111111111111111111111111111111111001100110011001100110011
        # Toggle mask = 1110001101111110001010001100010000100111000110110101101000101101
        # Total taxonomy nodes: 21101
        # Table size: 1399914077
        # Table capacity: 2000000000
        # Min clear hash value = 13727554816041021440

        error

        kraken2-inspect: database ("minikraken2_v1_8GB/") does not contain necessary file taxo.k2d
        """
        kraken2_inspect = self.which('kraken2-inspect')
        db_prefix = os.path.basename(self.kraken2_db.rstrip('/'))
        cmd = '%s --threads %d --db %s --skip-counts' % (kraken2_inspect, self.threads, self.kraken2_db)
        log = os.path.join(self.outdir, db_prefix + '.checker.log')
        logging.info('Checking Kraken2 database: %s' % self.kraken2_db)
        ## run program
        with open(log, 'wt') as fo:
            subprocess.run(cmd, shell=True, stdout=fo, stderr=fo)
        ## outdir, header
        db_content = ['kraken2-inspect log', ]
        with open(log, 'rt') as fi:
            header = next(fi)
            db_content.append(header.rstrip())
            for line in fi:
                db_content.append(line.rstrip())

        if header.startswith('kraken2-inspect'):
            sys.exit(header)
        else:
            logging.info('\n'.join(db_content) + '\n')


    def run(self):
        """Run Kraken2 program
        output
        """
        self.install_kraken2()

        ## cmd
        cmd = '%s %s --threads %d --use-names --db %s --output %s --report %s %s' % (self.kraken2_exe, self.args,
            self.threads, self.kraken2_db, self.kraken2_output, self.kraken2_report, self.input)
        
        # sys.exit(cmd)
        ## run
        logging.info('Running Kraken2')
        if os.path.exists(self.kraken2_output) and os.path.exists(self.kraken2_report) and self.overwrite is False:
            logging.info('Kraken2 output exists, skipped ...')
        else:
            with open(self.kraken2_log, 'wt') as fo:
                subprocess.run(cmd, shell=True, stdout=fo, stderr=fo)
    
        ## check output
        log_content = ['kraken2 log', ]
        with open(self.kraken2_log, 'rt') as fi:
            header = next(fi).rstrip()
            log_content.append(header.rstrip())
            for line in fi:
                log_content.append(line.rstrip())

        if header == 'Loading database information... done.':
            logging.info('\n'.join(log_content) + '\n')
        else:
            sys.exit('Kraken2 failed, see log : %s' % log)


    def stat(self):
        """Overall statistics of Kraken2
        the log file
        : total 
        : classified
        : unclassified
        """
        kraken2_log = self.kraken2_log
        if os.path.exists(kraken2_log):
            pass
        else:
            self.run()
            # sys.exit('Kraken2 output not found, please use Kraken2().run() first')
        ##
        num_total = 1
        with open(kraken2_log, 'rt') as fi:
            for line in fi:
                num = line.strip().split(' ')[0]
                if 'processed' in line:
                    num_total = num
                elif 'sequences classified' in line:
                    num_c = num
                elif 'unclassified' in line:
                    num_u = num
                else:
                    pass
        num_c_pct = int(num_c) / int(num_total) * 100
        # num_u_pct = 1 #'%.2f' % float(num_u) / float(num_total)
        ## to screen
        stat_log = '%s : %s : %.2f%% of %s sequences classified: ' % (self.prefix, 
            num_c, float(num_c_pct), num_total)
        logging.info(stat_log)
        ## output
        return [num_c, num_total]


    def report(self, topN=10):
        """Inspect the output of kraken2 report file
        pandas table
        """
        assert isinstance(topN, int)
        ## run 
        if os.path.exists(self.kraken2_report):
            pass
        else:
            self.run()
        ## report
        df = pd.read_csv(self.kraken2_report, '\t', 
            names=['pct', 'reads_in_clade', 'reads_in_tax', 'code', 'taxid', 'name'])
        ## sort by reads_in_tax
        df = df.sort_values(['reads_in_tax'], ascending=False).reset_index()
        ## remove white spaces
        df['name'] = df['name'].str.strip()
        ## add prefix
        df['sample'] = self.prefix

        ## total reads
        num_c, num_total = self.stat()

        ## adt pct
        df['reads_hit'] = num_c
        df['reads_total'] = num_total
        df['hit_pct'] = df['reads_in_tax'] / int(num_c) * 100

        ## sub-sample
        df_rpt = df.loc[0:topN, ['sample', 'name', 'reads_in_tax', 'hit_pct', 'reads_hit', 'reads_total']]
        with pd.option_context('display.expand_frame_repr', False, 'display.max_colwidth', 100):
            print(df_rpt)
        return df_rpt



# fq = '/home/wangming/work/yu_2019/projects/20190312_zp_goldclip/results/goldclip_v2/06.unmap_reads/demo.fq'
# out = '/home/wangming/work/yu_2019/projects/20190312_zp_goldclip/results/goldclip_v2/06.unmap_reads/demo'
# db = '/home/wangming/data/custom_db/kraken2'

# Kraken2(fq, out, db).install_kraken2()
# Kraken2(fq, out, db).install_std_db()
# Kraken2(fq, out, db).db_checker()
# Kraken2(fq, out, db).run()
# Kraken2(fq, out, db).report()
# Kraken2(fq, out, db).stat()