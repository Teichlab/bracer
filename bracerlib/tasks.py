from __future__ import print_function

import csv

import six
import matplotlib as mpl

from bracerlib import io, core
from bracerlib.io import check_binary

mpl.use('pdf')
import re
import seaborn as sns
from matplotlib import pyplot as plt
from bracerlib import base_dir
from bracerlib import bracer_func
from configparser import ConfigParser, NoOptionError
import argparse
import sys
import os
import subprocess
import glob
import shutil
from collections import defaultdict, Counter
from time import sleep
import warnings
import pickle
from prettytable import PrettyTable
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import itertools
import pdb
import numpy as np
from numpy import percentile, array
from matplotlib.colors import hex2color, rgb2hex
import random
import copy
import colorsys

class TracerTask(object):

    base_parser = argparse.ArgumentParser(add_help=False, 
                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    base_parser.add_argument('--ncores', '-p', metavar="<CORES>", 
                help='Number of processor cores to use', type=int, default=1)
    base_parser.add_argument('--config_file', '-c', metavar="<CONFIG_FILE>", 
                help='Config file to use', default=None)
    base_parser.add_argument('--resource_dir', metavar="<RESOURCE_DIR>",
                help='Root directory for resources', default=None)

    config = None

    def run(self):
        pass

    def get_binary(self, name):
        tool_key = name.lower() + '_path'
        user_path = None
        if self.config.has_option('tool_locations', tool_key):
            user_path = self.resolve_relative_path(
                self.config.get('tool_locations', tool_key))
        return check_binary(name, user_path)

    def read_config(self, config_file):
        # First look for environmental variable
        if not config_file:
            config_file = os.environ.get('BRACER_CONF', None)
            if config_file is not None:
                config_file = os.path.expanduser(config_file)
                if not os.path.isfile(config_file):
                    config_file = None
        # Then check the default location
        if not config_file:
            config_file = os.path.expanduser('~/.bracerrc')
            if not os.path.isfile(config_file):
                print("Config file not found at ~/.bracerrc."
                    " Using default bracer.conf in repo...")
                config_file = os.path.join(base_dir, 'bracer.conf')
        bracer_func.check_config_file(config_file)
        config = ConfigParser()
        config.read(config_file)

        return config

    def resolve_relative_path(self, path):
        if not path.startswith("/"):
            base_directory = os.path.abspath(os.path.dirname(__file__))
            full_path = os.path.normpath(
                "/{}/../{}".format(base_directory, path))
        else:
            full_path = path

        return full_path

    def print_cell_summary(self, cell, output_file, loci):
        out_file = open(output_file, 'w')
        out_file.write(
            '------------------\n{name}\n------------------\n'.format(
                                                    name=cell.name))
        
        # Summarise the productive/total recombinants
        for l in loci:
            out_file.write('BCR_{locus} recombinants: {summary}\n'.format(
                locus=l, summary=cell.summarise_productivity(l)))
        
        out_file.write('\n\n')
        
        for l in loci:
            out_file.write("#BCR_{locus}#\n".format(locus=l))
            rs = cell.recombinants["BCR"][l]
            if rs is None:
                out_file.write("No BCR_{} recombinants found\n\n".format(l))
            else:
                for r in rs:
                    out_file.write(r.get_summary())
                    out_file.write("\n\n")
        
        out_file.close()

    def die_with_empty_cell(self, cell_name, output_dir, species):
        print("##No recombinants found##")
        cell = core.Cell(cell_name, None, is_empty=True, species=species, 
                                                        loci=self.loci)
        
        # Save cell in a pickle
        unfiltered_summary = \
            "{}/unfiltered_BCR_seqs/unfiltered_BCRs.txt".format(self.output_dir)
        pickle_file = \
            "{}/unfiltered_BCR_seqs/{}.pkl".format(self.output_dir, cell.name)
        filtered_summary = \
            "{}/filtered_BCR_seqs/filtered_BCRs.txt".format(self.output_dir)
        filtered_pickle = \
            "{}/filtered_BCR_seqs/{}.pkl".format(self.output_dir, cell.name)

        self.print_cell_summary(cell, unfiltered_summary, self.loci)
            

        with open(pickle_file, 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)

        cell.filter_recombinants()
        self.print_cell_summary(cell, filtered_summary, self.loci)
                                                                            
        with open(filtered_pickle, 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)
        
        exit(0)
    
    def get_species_root(self, species, root=None):
        if root is None:
            resources_root = os.path.join(base_dir, 'resources', species)
        else:
            resources_root = os.path.join(root, species)
        assert os.path.isdir(resources_root), "Species not found in resources"
        return (resources_root)


    def get_rscript_path(self):
        rscript_path = os.path.join(base_dir, 'bracerlib/lineage.R')
        return(rscript_path)
        
    def get_available_species(self, root=None):
        if root is None:
            resources_dir = os.path.join(base_dir, 'resources')
        else:
            resources_dir = root
        species_dirs = next(os.walk(resources_dir))[1]
        return(species_dirs)


class Assembler(TracerTask):

    def __init__(self, **kwargs):
        if not kwargs:
            
            # get list of all available species in resources
            
            parser = argparse.ArgumentParser(
                description='Reconstruct BCR sequences from RNAseq reads for a '
                            'single cell', parents=[self.base_parser], 
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
            parser.add_argument('--resume_with_existing_files', '-r',
                                help='look for existing intermediate files and '
                                'use those instead of starting from scratch', 
                                action="store_true")
            parser.add_argument('--assembled_file', help='Fasta file containing '
                                'already assembled sequences for the cell. '
                                'Providing this file skips the alignment and '
                                'assembly steps', default=False)
            parser.add_argument('--species', '-s', 
                                help='Species to use for reconstruction', 
                                choices=self.get_available_species(), 
                                default='Hsap')
            parser.add_argument('--loci',
                                help='Space-separated list of loci to reconstruct '
                                'for receptor', default=['H','K', 'L'], nargs = '+')
            parser.add_argument('--single_end', help='set this if your sequencing '
                                'data are single-end reads', action="store_true")
            parser.add_argument('--fragment_length',
                                help='Estimated average fragment length in the '
                                'sequencing library.'
                                ' Used for Kallisto quantification. REQUIRED for '
                                'single-end data.', default=False)
            parser.add_argument('--fragment_sd',
                                help='Estimated standard deviation of average '
                                'fragment length in the sequencing library.'
                                ' Used for Kallisto quantification. REQUIRED '
                                'for single-end data.', default=False)
            parser.add_argument('--max_junc_len',
                                help='Maximum permitted length of junction '
                                'string in recombinant identifier. '
                                'Used to filter out artefacts.', default=100)
            parser.add_argument('--no_trimming',
                                help='Do not trim raw reads to remove adapter sequences '
                                'and low quality reads.', action = "store_true")
            parser.add_argument('--keep_trimmed_reads',
                                help='Do not delete the output files from the trimming step.',
                                action = "store_true")
            parser.add_argument('cell_name', metavar="<CELL_NAME>", 
                                help='name of cell for file labels')
            parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                help='directory for output as <output_dir>/<cell_name>')
            parser.add_argument('fastq1', metavar="<FASTQ1>",
                                help='first fastq file', nargs='?')                        
            parser.add_argument('fastq2', metavar="<FASTQ2>",
                                help='second fastq file', nargs='?')

            args = parser.parse_args(sys.argv[2:])
            
            resource_dir = args.resource_dir
            self.cell_name = args.cell_name
            self.fastq1 = args.fastq1
            self.single_end = args.single_end
            self.fastq2 = args.fastq2
            self.ncores = str(args.ncores)
            self.assembled_file = args.assembled_file
            self.species = args.species
            self.resume_with_existing_files = args.resume_with_existing_files
            self.fragment_length = args.fragment_length
            self.fragment_sd = args.fragment_sd
            self.output_dir = args.output_dir
            self.loci = args.loci
            self.max_junc_len = args.max_junc_len
            self.no_trimming = args.no_trimming
            self.keep_trimmed_reads = args.keep_trimmed_reads
            config_file = args.config_file
            

        else:
            resource_dir = kwargs.get('resource_dir')
            self.cell_name = kwargs.get('cell_name')
            self.fastq1 = kwargs.get('fastq1')
            self.fastq2 = kwargs.get('fastq2')
            self.ncores = kwargs.get('ncores')
            self.assembled_file = kwargs.get('assembled_file')
            self.species = kwargs.get('species')
            self.resume_with_existing_files = kwargs.get(
                            'resume_with_existing_files')
            self.output_dir = kwargs.get('output_dir')
            self.single_end = kwargs.get('single_end')
            self.fragment_length = kwargs.get('fragment_length')
            self.fragment_sd = kwargs.get('fragment_sd')
            self.loci = kwargs.get('loci')
            self.max_junc_len = kwargs.get('max_junc_len')
            self.no_trimming = kwargs.get('no_trimming')
            self.keep_trimmed_reads = kwargs.get('keep_trimmed_reads')
            config_file = kwargs.get('config_file')

        self.trimmed_fastq1 = None
        self.trimmed_fastq2 = None
        self.config = self.read_config(config_file)
        self.species_root = self.get_species_root(self.species,
                                            root=resource_dir)
        
        if not self.assembled_file:
            self.assembled_file = None
 
        # Check that FASTA file containing assembled sequences exists
        # if running with --assembled_file
        if self.assembled_file:
             if not os.path.isfile(self.assembled_file):
                   raise OSError('2', 'Fasta file containing assembled '
                                'sequences not found. Please provide FASTA '
                                'file or run Assemble without --assembled_file', 
                                self.assembled_file)
        else:
            assert self.fastq1, \
                ('No FASTQ specified. Either provide FASTQ or provide '
                'assembled sequences in FASTA format with the '
                '--assembled_file option.')

         # Check FASTQ files exist
        if not self.assembled_file:
            if not os.path.isfile(self.fastq1):
                raise OSError('2', 'FASTQ file not found', self.fastq1)
            if not self.single_end and self.fastq2:
                if not os.path.isfile(self.fastq2):
                    raise OSError('2', 'FASTQ file not found', self.fastq2)

        # Check the FASTQ config is correct
        if not self.assembled_file:
            if not self.single_end:
                assert self.fastq2, \
                    ('Only one fastq file specified. ' 
                    'Either set --single_end or provide second fastq.')
            else:
                self.fastq2 = None
                if self.fastq2:
                    print('Two fastq files given with --single-end option.')
                    print('Ignoring second file.')
                assert self.fragment_length and self.fragment_sd, \
                    ('Must specify estimated average fragment length '
                    '(--fragment_length) and standard deviation (--fragment_sd) '
                    'for use with single-end data')
                assert self.fragment_length, \
                    'Must specify estimated average fragment length \
                    (--fragment_length) for use with single-end data'
                assert self.fragment_sd, \
                    'Must specify estimated fragment length standard deviation \
                    (--fragment_sd) for use with single-end data'

        
    def run(self, **kwargs):

        # Set-up output directories
        root_output_dir = os.path.abspath(self.output_dir)
        io.makeOutputDir(root_output_dir)
        self.output_dir = root_output_dir + "/" + self.cell_name

        io.makeOutputDir(self.output_dir)

        data_dirs = ['aligned_reads', 'Trinity_output', 'IgBLAST_output', 
                    'BLAST_output', 'unfiltered_BCR_seqs', 
                    'expression_quantification', 'filtered_BCR_seqs']
        if not self.no_trimming:
            data_dirs.append('trimmed_reads')

        for d in data_dirs:
            io.makeOutputDir("{}/{}".format(self.output_dir, d))

        # Perform BraCeR's core functions
        if not self.assembled_file:
            if not self.no_trimming:
                self.trim_reads()
            self.align()
            self.de_novo_assemble()
        
        self.blast()
        cell = self.ig_blast()

        if self.fastq1:
            self.quantify(cell)
        
        
        unfiltered_fasta_filename = \
            "{}/unfiltered_BCR_seqs/{}_BCRseqs.fa".format(self.output_dir, 
                                                            self.cell_name)
        unfiltered_cell_summary_file = \
            "{}/unfiltered_BCR_seqs/unfiltered_BCRs.txt".format(self.output_dir)

        unfiltered_pickle = \
            "{}/unfiltered_BCR_seqs/{}.pkl".format(self.output_dir, 
                                                    self.cell_name)
        filtered_fasta_filename = \
            "{}/filtered_BCR_seqs/{}_BCRseqs.fa".format(self.output_dir, 
                                                        self.cell_name)
        filtered_cell_summary_file = \
            "{}/filtered_BCR_seqs/filtered_BCRs.txt".format(self.output_dir)

        filtered_pickle = \
            "{}/filtered_BCR_seqs/{}.pkl".format(self.output_dir, 
                                                    self.cell_name) 
                                                                        
        unfiltered_fasta_file = open(unfiltered_fasta_filename, 'w')
        unfiltered_fasta_file.write(cell.get_fasta_string())
        unfiltered_fasta_file.close()

        self.print_cell_summary(cell, unfiltered_cell_summary_file, self.loci)
        
        # Assign isotype and bgcolor
        ranked_recs = cell.rank_recombinants()
        isotype = cell.determine_isotype(ranked_recs)
        bgcolor = cell.assign_bgcolor(isotype)
        cell.bgcolor = bgcolor
        cell.isotype = isotype

        # Create database dictionary
        cell.databasedict = cell.get_database_for_locus(self.loci)

        # Save cell in a pickle
        with open(unfiltered_pickle, 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)


        # Filter recombinants
        print("##Filtering by read count##")
        cell.filter_recombinants()
        filtered_fasta_file = open(filtered_fasta_filename, 'w')
        filtered_fasta_file.write(cell.get_fasta_string())
        filtered_fasta_file.close()
        self.print_cell_summary(cell, filtered_cell_summary_file, self.loci)

        # Save cell in a pickle
        with open(filtered_pickle, 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)


    def trim_reads(self):
        trim_galore = self.get_binary('trim_galore')
        cutadapt = self.get_binary('cutadapt')

        self.trimmed_fastq1, self.trimmed_fastq2 = bracer_func.run_trim_galore(
                    trim_galore, cutadapt, self.output_dir, self.cell_name, 
                    self.fastq1, self.fastq2, self.resume_with_existing_files, 
                    self.single_end)


    def align(self):
        bowtie2 = self.get_binary('bowtie2')

        synthetic_genome_path = os.path.join(self.species_root, 
                                'combinatorial_recombinomes')
        bowtie2_build = self.get_binary('bowtie2-build')
        # Align with Bowtie 2
        bracer_func.bowtie2_alignment(bowtie2, self.ncores, self.loci, 
                    self.output_dir, self.cell_name, synthetic_genome_path, 
                    self.fastq1, self.fastq2, self.resume_with_existing_files, 
                    self.single_end, bowtie2_build, self.trimmed_fastq1,
                    self.trimmed_fastq2)
        print()

    def de_novo_assemble(self):

        try:
            trinity = self.get_binary('trinity')
        except OSError:
            trinity = self.get_binary('Trinity')

        # Trinity version
        if not self.config.has_option('trinity_options', 'trinity_version'):
            try:
                subprocess.check_output([trinity, '--version'])
            except subprocess.CalledProcessError as err:
                if re.search('v2', err.output.decode('utf-8')):
                    self.config.set('trinity_options', 'trinity_version', '2')
                else:
                    self.config.set('trinity_options', 'trinity_version', '1')

        if self.config.has_option('trinity_options', 'trinity_grid_conf'):
            trinity_grid_conf = self.resolve_relative_path(
                self.config.get('trinity_options', 'trinity_grid_conf'))
        else:
            trinity_grid_conf = False

        # Is Trinity version compatible with --no_normalize_reads argument
        no_normalise = False
        trinity_version = self.config.get('trinity_options', 'trinity_version')
        
        try:
            subprocess.check_output([trinity, '--version'])
        except subprocess.CalledProcessError as err:
            first_line = err.output.decode('utf-8').split("\n")[0]
            if re.search('v2.3', first_line) or re.search('v2.4', first_line):
                no_normalise = True
            else:
                no_normalise = False

        # De novo assembly with Trinity
        trinity_JM = self.config.get('trinity_options', 'max_jellyfish_memory')


        successful_files = bracer_func.assemble_with_trinity(trinity, 
            self.loci, self.output_dir, self.cell_name, self.ncores, 
            trinity_grid_conf, trinity_JM, trinity_version, 
            self.resume_with_existing_files, self.single_end, 
            self.species, no_normalise)

        if len(successful_files) == 0:
            print("No successful Trinity assemblies")
            self.die_with_empty_cell(self.cell_name, self.output_dir, 
                                                        self.species)
        print()


    def ig_blast(self):
        igblastn = self.get_binary('igblastn')

        # Reference data locations

        igblast_index_location = os.path.join(self.species_root, 'igblast_dbs')
        imgt_seq_location = os.path.join(self.species_root, 'raw_seqs')
        igblast_seqtype = 'Ig'

        # IgBlast of assembled contigs
        bracer_func.run_IgBlast(igblastn, self.loci, self.output_dir, 
                    self.cell_name, igblast_index_location, igblast_seqtype, 
                    self.species, self.resume_with_existing_files, 
                    self.assembled_file)
        print()
        
        # Run IgBlast of assembled contigs using IMGT-gapped references
        gapped_index_location = os.path.join(self.species_root,
                                        'imgt_gapped_resources/igblast_dbs')
        bracer_func.run_IgBlast_IMGT_gaps_for_cell(igblastn, self.loci, 
                    self.output_dir, self.cell_name, igblast_index_location, 
                    gapped_index_location, igblast_seqtype, self.species, 
                    self.assembled_file)
        self.create_changeo_db()


        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cell = io.parse_IgBLAST(self.loci, self.output_dir, self.cell_name, 
                    imgt_seq_location, self.species, self.assembled_file, 
                    self.max_junc_len)
            if cell.is_empty:
                self.die_with_empty_cell(self.cell_name, self.output_dir, 
                                                            self.species)

        return cell


    def create_changeo_db(self):
        """Creates Change-O database from IgBlast result files after alignment
        to imgt-gapped sequences for CDR3 detection"""
        try:
            MakeDb =  self.config.get('tool_locations', 'changeo_path') + "/MakeDb.py"
            if not os.path.is_file(MakeDb):
                MakeDb = "MakeDb.py"
        except:
            MakeDb = "MakeDb.py"

        gapped_seq_location = os.path.join(self.species_root,
                            'imgt_gapped_resources/raw_seqs')
        ungapped_seq_location = os.path.join(self.species_root, 'raw_seqs')

        # Run MakeDb from Change-O toolkit
        for locus in self.loci:
            bracer_func.run_MakeDb_for_cell(MakeDb, locus, self.output_dir, 
                        self.species, gapped_seq_location, 
                        ungapped_seq_location, self.cell_name)


    def blast(self):
        blastn = self.get_binary('blastn')

        # Reference data locations
        blast_index_location = os.path.join(self.species_root, 'igblast_dbs')
        imgt_seq_location = os.path.join(self.species_root, 'raw_seqs')

        # BLAST assembled contigs
        bracer_func.run_Blast(blastn, self.loci, self.output_dir, 
                    self.cell_name, blast_index_location, self.species, 
                    self.resume_with_existing_files, self.assembled_file)
        print()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            isotype = io.parse_BLAST(self.loci, self.output_dir, 
                        self.cell_name, self.species, self.assembled_file)

        return isotype


    def quantify(self, cell):
        kallisto = self.get_binary('kallisto')
        
        if not self.config.has_option('kallisto_transcriptomes', self.species):
            raise OSError("No transcriptome reference specified for {species}. "
                        "Please specify location in config file."
                        .format(species=self.species))
        else:
            kallisto_base_transcriptome = self.resolve_relative_path(
                self.config.get('kallisto_transcriptomes', self.species))
            if not os.path.isfile(kallisto_base_transcriptome):
                raise OSError('2', 'Transcriptome reference not found', 
                                kallisto_base_transcriptome)

        #Quantification with kallisto
        bracer_func.quantify_with_kallisto(
                kallisto, cell, self.output_dir, self.cell_name, 
                kallisto_base_transcriptome, self.fastq1, self.fastq2, 
                self.ncores, self.resume_with_existing_files, self.single_end, 
                self.fragment_length, self.fragment_sd, self.trimmed_fastq1,
                self.trimmed_fastq2, self.keep_trimmed_reads)
        print()

        counts = bracer_func.load_kallisto_counts(
            "{}/expression_quantification/abundance.tsv".format(
                                                self.output_dir))
        
        for receptor, locus_dict in six.iteritems(cell.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None:
                    for rec in recombinants:
                        tpm = counts["BCR"][locus][rec.contig_name]
                        rec.TPM = tpm
                    

class Summariser(TracerTask):

    def __init__(self, **kwargs):

        if not kwargs:
            parser = argparse.ArgumentParser(description="Summarise set of cells \
                                             with reconstructed BCR sequences",
                                             parents=[self.base_parser])
            parser.add_argument('--species', '-s',
                                help='Species to use for reconstruction',
                                choices=self.get_available_species(), 
                                default='Hsap')
            parser.add_argument('--loci', help="Space-separated list of loci to \
                                summarise for receptor", default=['H','K', 'L'], 
                                nargs = '+')
            parser.add_argument('--use_unfiltered', '-u', 
                                help='Use unfiltered recombinants', 
                                action="store_true")
            parser.add_argument('--graph_format', '-f', metavar="<GRAPH_FORMAT>", 
                                help='Graphviz output format [pdf]', default='svg')
            parser.add_argument('--no_networks', help='Skip attempts to draw network \
                                graphs', action = "store_true")
            parser.add_argument('--IGH_networks', help='Base clonality only \
                                on heavy chain', action = "store_true")
            parser.add_argument('--dist', metavar="<DISTANCE>",
                                help='Distance value to use for clonal inference. ' \
                                'Heavily mutated datasets may require a higher ' \
                                'distance value, whereas datasets enriched ' \
                                'in naive B cells may require a lower value.', 
                                default='0.2')
            parser.add_argument('--include_multiplets', help='Do not exclude \
                                cells containing more than two recombinants \
                                for a locus from downstream analyses, \
                                including networks and clonotype analysis', 
                                action = "store_true")
            parser.add_argument('--infer_lineage', 
                                help='Construct lineage trees for clone groups \
                                shown in clonal network', action = "store_true")
            parser.add_argument('dir', metavar="<DIR>", help='Directory containing \
                                subdirectories for each cell to be summarised')
            args = parser.parse_args(sys.argv[2:])

            resource_dir = args.resource_dir
            self.root_dir = os.path.abspath(args.dir)
            self.graph_format = args.graph_format
            self.use_unfiltered = args.use_unfiltered
            self.draw_graphs = not args.no_networks
            self.IGH_networks = args.IGH_networks
            self.no_multiplets = not args.include_multiplets
            self.loci = args.loci
            self.species = args.species
            config_file = args.config_file
            self.dist = args.dist
            self.infer_lineage = args.infer_lineage
        else:
            resource_dir = kwargs.get('resource_dir')
            self.use_unfiltered = kwargs.get('use_unfiltered')
            self.root_dir = os.path.abspath(kwargs.get('root_dir'))
            self.draw_graphs = not (kwargs.get('no_networks'))
            self.IGH_networks = kwargs.get('IGH_networks')
            self.no_multiplets = kwargs.get('no_multiplets')
            self.graph_format = kwargs.get('graph_format')
            self.loci = kwargs.get('loci')
            self.species = kwargs.get('species')
            config_file = kwargs.get('config_file')
            self.dist = kwargs.get('dist')
            self.infer_lineage = kwargs.get('infer_lineage')

        # Read config file
        self.config = self.read_config(config_file)
        self.species_dir = self.get_species_root(self.species,
                                            root=resource_dir)

        
        
    def run(self):

        if self.draw_graphs:
            dot = self.get_binary("dot")
            neato = self.get_binary("neato")

        else:
            dot = ""
            neato = ""

        if self.infer_lineage:
            dnapars = self.get_binary("dnapars")
            rscript = self.get_binary("rscript")
            R = self.get_binary("R")
                    
        
        cells = {}
        empty_cells = []
        subdirectories = next(os.walk(self.root_dir))[1]

        if self.use_unfiltered:
            pkl_dir = "unfiltered_BCR_seqs"
            outdir = "{}/unfiltered_BCR_summary".format(self.root_dir)

        else:
            pkl_dir = "filtered_BCR_seqs"
            outdir = "{}/filtered_BCR_summary".format(self.root_dir)

        io.makeOutputDir(outdir)

        if self.infer_lineage:
            lineage_dir = "{}/lineage_trees".format(outdir)
            io.makeOutputDir(lineage_dir)

        outfile = open("{}/BCR_summary.txt".format(outdir), 'w')
        length_filename_root = "{}/reconstructed_lengths_BCR".format(outdir)
        isotype_filename_root = "{}/isotypes_BCR".format(outdir)

        for d in subdirectories:
            cell_pkl = "{root_dir}/{d}/{pkl_dir}/{d}.pkl".format(
                        pkl_dir=pkl_dir, d=d, root_dir=self.root_dir)
            if os.path.isfile(cell_pkl):
                with open(cell_pkl, 'rb') as pkl:
                    cl = pickle.load(pkl)
                cells[d] = cl
                if cl.is_empty or cl.missing_loci_of_interest(self.loci):
                    empty_cells.append(d)
                               
        
        # Count cells with productive chains for each locus and for each 
        # possible pair and write to summary file
        cell_recovery_count = self.write_reconstruction_statistics(outfile, 
                                                        self.loci, cells)

        # Make all recombinant table and write to summary file
        t = self.make_all_recombinant_table(self.loci, cells)
        outfile.write(t.get_string())
        outfile.write("\n")

        # Name potential multiplets (cells with more
        # than two recombinants for a locus
        multiplets = []
        if self.use_unfiltered:
            multiplets = self.detect_multiplets(self.no_multiplets, outfile, 
                                                        cells, self.loci)
        else:
            #Identify multiplets if using filtered recombinants
            unf_cells = {}
            unf_pkl_dir = "unfiltered_BCR_seqs"
            unf_outdir = "{}/unfiltered_BCR_summary".format(self.root_dir)
            for d in subdirectories:
                unf_pkl = "{root_dir}/{d}/{pkl_dir}/{d}.pkl".format(
                    pkl_dir=unf_pkl_dir, d=d, root_dir=self.root_dir)
                if os.path.isfile(unf_pkl):
                    with open(unf_pkl, 'rb') as unf_pkl:
                        cl = pickle.load(unf_pkl)
                    unf_cells[d] = cl
            multiplets = self.detect_multiplets(self.no_multiplets, outfile,
                                                        unf_cells, self.loci)


        # Create Change-O db file for filtered multiplets if no_multuplets is True
        if self.no_multiplets and len(multiplets) > 0:
            self.create_multiplet_database_file(outdir, self.loci, cells, 
                                                            multiplets)

        # Delete likely multiplets from downstream analyses unless --include_multiplets 
        if self.no_multiplets and len(multiplets) > 0:
            for cell_name in multiplets:
                del cells[cell_name]


        # Make full length statistics table and plot proportion of sequences 
        # that are full-length
        (full_length_counter, full_length_prod_counter, productive_counter, 
            total_counter) = self.count_full_length_sequences(self.loci, cells)
        
        (full_length_statistics_table, all_dict, prod_dict) = \
            self.make_full_length_statistics_table(self.loci, total_counter, 
            productive_counter, full_length_counter, full_length_prod_counter)
        
        outfile.write(full_length_statistics_table)
        outfile.write("\n")
        self.plot_full_length_sequences_proportions(all_dict, prod_dict, 
                                                      self.loci, outdir)
         
        # Report cells with two productive chains from same locus
        outstring = self.two_productive_chains_per_locus(self.loci, cells)
      
        if len(outstring) > 0:
            outfile.write(outstring)
      
        # Report cells with productive kappa and lambda chain
        outstring = self.report_kappa_lambda_cells(self.loci, cells)

        if len(outstring) > 0:
            outfile.write(outstring)

        # Make initial clonal assignments for B cells using Change-O 
        # DefineClones bygroup
        try:
            DefineClones = (self.config.get('tool_locations', 'changeo_path') + 
                                                            "/DefineClones.py")
            if not os.path.isfile(DefineClones):
                DefineClones = "DefineClones.py"
        except:
            DefineClones = "DefineClones.py"

        for locus in self.loci:
            self.make_changeo_input(outdir, locus, cells)
            bracer_func.run_DefineClones(DefineClones, locus, outdir, 
                                            self.species, self.dist)
            print()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")

        # Read Change-O result files and define clone groups for each locus
        (clones, cell_clones, cell_contig_clones) = self.read_changeo_results(
                                                            self.loci, outdir)

        # Get H clone groups consisting of 2 or more cells
        multiple_clones_H = dict()
        for clone, cells_info in six.iteritems(clones["H"]):
            if len(cells_info.keys()) > 1:
                multiple_clones_H[clone] = cells_info

        cells_with_clonal_H = []
        for clone, cells_info in six.iteritems(multiple_clones_H):
            for cell_name, contig_name in six.iteritems(multiple_clones_H[clone]):
                if not cell_name in cells_with_clonal_H:
                    cells_with_clonal_H.append(cell_name)

        # Create Change-O compatible database file containing all sequences 
        #with info from BraCeR (not IMGT-gapped)
        self.create_database_file(outdir, self.loci, cells, cell_contig_clones)

        # Create database file in the Change-O format with IMGT-gaps
        self.create_IMGT_gapped_db(outdir, cells, self.loci, cell_contig_clones)

        # Make isotype usage table and plot isotype distributions
        isotype_counter = self.count_isotype_usage(cells)
        (header, outstring) = self.make_isotype_table(cell_recovery_count, 
                                                          isotype_counter)
        outfile.write(header)
        outfile.write(outstring)
        outfile.write("\n")
        self.plot_isotype_distributions(isotype_counter, outdir)

        # Plot lengths of reconstructed sequences
        self.plot_length_distributions(self.loci, length_filename_root, cells)

        # Delete cells with no reconstructed sequences                
        for cell_name in empty_cells:
            del cells[cell_name]
                
        # Make clonotype networks
        network_colours = io.read_colour_file(os.path.join(
                    self.species_dir, "colours.csv"), receptor_name="BCR")

        component_groups, G = bracer_func.draw_network_from_cells(cells, outdir, 
                           self.graph_format, dot, neato, self.draw_graphs, 
                           self.loci, network_colours, cell_contig_clones, 
                           cells_with_clonal_H, self.no_multiplets, 
                           self.IGH_networks)
        
        # Print component groups to the summary
        outfile.write(
            "\n#Clonotype groups#\n"
            "This is a text representation of the groups shown in clonotype_"
             "network_with_identifiers.pdf.\n"
            "It does not exclude cells that only share heavy chain and not "
            "light chain if Summarise is run with --IGH_networks.\n\n")

        for g in component_groups:
            outfile.write(", ".join(g))
            outfile.write("\n")
            
        outfile.write("\n\n")

        # Reconstruct lineages - optional
        """These steps use MakeDb and CreateGermlines of the Change-O toolkit.
        Change-O reference: Gupta NT*, Vander Heiden JA*, Uduman M, Gadala-Maria D
        Yaari G, Kleinstein SH. Change-O: a toolkit for analyzing large-scale B cell
        immunoglobulin repertoire sequencing data. Bioinformatics 2015;
        doi: 10.1093/bioinformatics/btv359"""

        if self.infer_lineage and len(component_groups) > 0:
            # Create IgBlast input files
            if not self.IGH_networks and not self.loci ==["H"]:
                # Find the clonal sequences that are shared in each component group
                clone_groups = self.get_clone_groups_for_component_group(self.loci, 
                                        cells, component_groups, cell_contig_clones)
                for locus in self.loci:
                    self.make_clone_igblast_input_for_concatenation(outdir, locus, 
                            cells, component_groups, clone_groups, cell_contig_clones)
            else:
                for locus in self.loci:
                    self.make_clone_igblast_input(outdir, locus, cells)

            for locus in self.loci:
                # Align clonal sequences using IgBlast and IMGT-gapped references
                self.IgBlast_germline_reconstruction(outdir, locus)

                # Create Change-O database from IMGT-gapped IgBlast results
                self.create_changeo_db(outdir, locus)

                # Modify database
                self.modify_changeo_db(outdir, locus, cells, cell_contig_clones)

                # Reconstruct germline sequences with CreateGermlines
                self.create_germline_sequences(outdir, locus, cells)

            if not self.IGH_networks and not self.loci == ["H"]:
                self.concatenate_lineage_databases(outdir, self.loci, 
                                        component_groups, clone_groups)

            # Run lineage reconstruction with Alakazam
            self.create_lineage_trees(outdir, self.species)


        # Print name of empty cells to the summary
        outfile.write("#Cells with no reconstructed sequences#\n")
        if len(empty_cells) == 0:
                    outfile.write("None")
        elif len(empty_cells) == 1:
            for i in empty_cells:
                outfile.write("".join(i))
        else:
            outfile.write(", ".join(empty_cells))
        
        outfile.write("\n\n")
        
        # Plot clonotype sizes
        x_range, clonotype_sizes = self.plot_clonotype_sizes(outdir, cells, 
                                                            self.loci, G)

        # Write clonotype sizes to text file
        with open("{}/clonotype_sizes.txt".format(outdir), 'w') as f:
            data = zip(x_range, clonotype_sizes)
            f.write("clonotype_size\tclonotype_count\n")
            for t in data:
                f.write("{}\t{}\n".format(t[0], t[1]))

        outfile.close()


    def get_clone_groups_for_component_group(self, loci, cells, 
                        component_groups, cell_contig_clones):
        """Returns dictionary with number of chains belonging to each 
        locus-specific clone group (determined by Change-O) in each component 
        group in the clonal network to find the locus-specific clone groups 
        of the most commonly shared chains within each component group"""
        
        counter = dict()
        i = 0
        for g in component_groups:
            i += 1
            counter[i] = dict()
            for l in loci:
                counter[i][l] = dict()
                for cell_name in g:
                    cell = cells[cell_name]
                    for rec in cell.recombinants["BCR"][l]:
                        contig_name = rec.contig_name + "_" + rec.identifier
                        if contig_name in cell_contig_clones[l][cell_name].keys():
                            clone = cell_contig_clones[l][cell_name][contig_name]
                            clone = clone + "_" + l
                            if not clone in counter[i][l].keys():
                                counter[i][l][clone] = 1
                            else:
                                counter[i][l][clone] += 1
                            
        # Find most highly shared heavy and light chain
        top_chain = dict()
        for i in counter.keys():
            top_chain[i] = dict()
            heavy_count = 0
            light_count = 0

            for l in loci:
                for clone, count in six.iteritems(counter[i][l]):
                    if l =="H":
                        if count > heavy_count:
                            heavy_count = count
                            top_chain[i]["heavy"] = clone
                    else:
                        if count > light_count:
                            light_count = count
                            top_chain[i]["light"] = clone
                    
        return (top_chain)
        

    def write_reconstruction_statistics(self, outfile, loci, cells):
        """Writes statistics for reconstruction of productive chains and 
        paired chains to summary file"""
        
        cell_recovery_count = dict()

        # Count cells with productive chains for each locus and possible pairs
        for l in loci:
            cell_recovery_count[l] = 0
        
        possible_pairs = ["".join(x) for x in itertools.combinations(loci, 2)]

        for p in possible_pairs:
            cell_recovery_count[p] = 0
        
        for cell_name, cell in six.iteritems(cells):
            prod_counts = dict()
            for l in loci:
                prod_counts[l] = cell.count_productive_recombinants(l)
                if prod_counts[l] > 0:
                    cell_recovery_count[l] += 1
            for p in possible_pairs:
                if prod_counts[p[0]] > 0 and prod_counts[p[1]] > 0:
                    cell_recovery_count[p] += 1

        total_cells = len(cells)

        for l in loci:
            count = cell_recovery_count[l]
            pc = round((count/float(total_cells))*100, 1)
            outfile.write(
                "BCR_{locus} reconstruction:\t{count} / {total} ({pc}%)\n".format(
                                locus=l, count=count, total=total_cells, pc=pc))
        outfile.write("\n")
        
        
        for p in possible_pairs:
            count = cell_recovery_count[p]
            pc = round((count/float(total_cells))*100, 1)
            if pc == "KL":
                outstring = ("Productive reconstruction of K and L:\t" +
                            "{count} / {total} ({pc}%)\n".format(p=p, 
                            count=count, total=total_cells, pc=pc))
                outfile.write(outstring)
            else:
                 outfile.write(
                        "Paired {p} productive reconstruction:\t{count} / {total} ({pc}%)\n".format(
                                p=p, count=count, total=total_cells, pc=pc))
            
        outfile.write("\n")

        return(cell_recovery_count)


    def make_all_recombinant_table(self, loci, cells):
        """Creates table with number of reconstructed recombinants per locus"""
        
        all_counters = defaultdict(Counter)
        prod_counters = defaultdict(Counter)
        for cell in cells.values():
            for l in loci:
                all_counters[l].update({cell.count_total_recombinants(l): 1})
                prod_counters[l].update({cell.count_productive_recombinants(l): 1})

        all_recombinant_counts = []

        for locus in all_counters:
            all_recombinant_counts = (all_recombinant_counts +
                    list(all_counters[locus].keys()))

        max_recombinant_count = max(all_recombinant_counts)

        table_header = ['', '0 recombinants', '1 recombinant', '2 recombinants']
        recomb_range = range(0, 3)
        if max_recombinant_count > 2:
            extra_header = [str(x) + " recombinants" for x in range(
                                        3, max_recombinant_count + 1)]
            table_header = table_header + extra_header
            recomb_range = range(0, max_recombinant_count + 1)
        t = PrettyTable(table_header)
        t.padding_width = 1
        t.align = "l"

        for counter_name in ['all_counters', 'prod_counters']:
            counter_type = counter_name.split("_")[0]
            counter_set = eval(counter_name)
            for l in loci:
                counter = counter_set[l]
                count_array = [counter[x] for x in recomb_range]
                total_with_at_least_one = sum(count_array[1:])
                if total_with_at_least_one > 0:
                    percentages = [''] + [" (" + str(round((float(x) \
                        / total_with_at_least_one) * 100)) + "%)" \
                        for x in count_array[1:]]
                else:
                    percentages = [''] + [" (N/A%)" for x in count_array[1:]]
                row = []
                for i in recomb_range:
                    row.append(str(count_array[i]) + percentages[i])
                label = '{} {}'.format(counter_type, l)
                t.add_row([label] + row)
        
        return(t)


    def detect_multiplets(self, no_multiplets, outfile, cells, loci):
        """Names cells with more than two recombinants for a locus"""
        
        multiplets = []
        outfile.write("\n\n##Cells with more than two recombinants for a locus##\n")
        if no_multiplets:
            outfile.write("\nThe following cells are likely multiplets or "
                "contaminated as they contain more than two recombinants for "
                "a locus, and were excluded from downstream analyses.\n")
        else:
             outfile.write("\nThe following cells are likely multiplets or "
             "contaminated as they contain more than two recombinants for "
             "a locus. Consider removing them from downstream analysis.\n")
        found_multi = False
        for cell in cells.values():
            if cell.has_excess_recombinants:
                outfile.write("*{}*\n".format(cell.name))
                multiplets.append(cell.name)
                for l in loci:
                    count = cell.count_total_recombinants(l)
                    outfile.write("BCR_{l}:\t{count}\n".format(l=l, count=count))
                outfile.write("\n")
                found_multi = True
        if not found_multi:
            outfile.write("None\n\n")

        return(multiplets)


    def get_quartiles(self, locus):
        
        fasta = os.path.join(self.species_dir, 'combinatorial_recombinomes',
                'BCR_{locus}.fa'.format(locus=locus))
        
        # Remove the start N padding and end C sequence from the lengths
        constant_fasta = os.path.join(self.species_dir, 'raw_seqs', 
                'BCR_{locus}_C.fa'.format(locus=locus))
        C_len = len(next(SeqIO.parse(constant_fasta, "fasta")))

        seq = str(next(SeqIO.parse(fasta, "fasta")).seq)
        if seq.startswith("N"):
            r = re.compile(r"^N+")
            N_len = r.search(seq).end()
        else:
            N_len = 0
        
        adj = C_len+N_len                                                               
                                                                            
        lengths = array([len(rec)-adj for rec in SeqIO.parse(fasta, "fasta")])
        quartiles = (percentile(lengths, 25), percentile(lengths, 75))

        return(quartiles)


    def plot_length_distributions(self, loci, length_filename_root, cells):
        """Plots length distributions of reconstructed sequences"""
        
        lengths = defaultdict(list)
        for cell in cells.values():
            for l in loci:
                lengths[l].extend(cell.get_trinity_lengths(l))

        quartiles = dict()
        for l in loci:
            q = self.get_quartiles(l)
            quartiles[l] = q

        for l in loci:
            q = quartiles[l]
            lns = lengths[l]
            if len(lns) > 1:
                plt.figure()
                plt.axvline(q[0], linestyle="--", color='k')
                plt.axvline(q[1], linestyle="--", color='k')
                sns.distplot(lns)
                sns.despine()
                plt.xlabel("BCR_{} reconstructed length (bp)".format(l))
                plt.ylabel("Density")
                plt.savefig("{}_{}.pdf".format(length_filename_root, l))
            if len(lns) > 0:
                with open("{}_{}.txt".format(length_filename_root,l), 'w') as f:
                    for l in sorted(lns):
                        f.write("{}\n".format(l))
            
        
    def two_productive_chains_per_locus(self, loci, cells):
        """Report cells with two (or more) productive chains from same locus"""

        outstring = ""
        for l in loci:
            double_rec = []
            cells_string = ""
            for cell_name, cell in six.iteritems(cells):
                prod_counts = dict()
                prod_counts[l] = cell.count_productive_recombinants(l)
                if prod_counts[l] > 1:
                    double_rec.append(cell_name)
            if len(double_rec) > 0 :
                cells_string = ', '.join(double_rec)
                locus_string = ("Cells with two productive {} chains: \n".format(l) +
                                                        cells_string + "\n\n")
                outstring += locus_string
        if len(outstring) != 0:
            outstring = ("\n\n##Cells with two (or more) productive "
                        "recombinants for a locus##\n\n" + outstring)

        return (outstring)


    def report_kappa_lambda_cells(self, loci, cells):
        """Report cells with both productive kappa and lambda chain"""
        outstring = ""
        if ("K" and "L") in loci:
            
            kappa_lambda = []
            cells_string = ""
            for cell_name, cell in six.iteritems(cells):
                store_cell = True
                for l in ["K", "L"]:
                    prod_counts = dict()
                    prod_counts[l] = cell.count_productive_recombinants(l)
                    if prod_counts[l] == 0:
                        store_cell = False
                if store_cell == True:
                    kappa_lambda.append(cell_name)
            if len(kappa_lambda) > 0 :
                cells_string = ', '.join(kappa_lambda)
                outstring = ("\n\n##Cells with productive K and L chain##\n\n" +
                        cells_string + "\n\n")

        return(outstring)


    def create_database_file(self, outdir, loci, cells, cell_contig_clones):
        """Creates a tab-delimited Change-O database containing all 
        reconstructed sequences and info from BraCeR using ungapped IMGT
        reference sequences, although CDR3 sequence and producivity info
        comes from alignment to IMGT-gapped references and parsing with
        Change-O if info is available."""

        changeo_string = ("CELL\tSEQUENCE_ID\tCONTIG_NAME\tLOCUS\tFUNCTIONAL\t"
                        "IN_FRAME\tSTOP\tINDELS\tV_CALL\tTOP_V_ALLELE\tD_CALL"
                        "\tTOP_D_ALLELE\tJ_CALL\tTOP_J_ALLELE\tSEQUENCE_INPUT"
                        "\tSEQUENCE_VDJ\tV_SEQ_START\tV_SEQ_LENGTH\t"
                        "D_SEQ_START\tD_SEQ_LENGTH\tJ_SEQ_START\tJ_SEQ_LENGTH"
                        "\tCDR3\tCDR3_LENGTH\tJUNCTION\tJUNCTION_LENGTH\tTPM\t"
                        "C_CALL\tISOTYPE\tCLONE\n")

        database_file = "{}/changeodb.tab".format(outdir)

        with open(database_file, 'w') as output:
            output.write(changeo_string)
            for cell_name, cell in six.iteritems(cells):
                databasedict = cell.databasedict
                for locus in loci:
                    recs = cell.recombinants["BCR"][locus]
                    if not recs is None:
                        for rec in recs:
                            C_gene = rec.C_gene
                            if C_gene == None:
                                C_gene = "None"
                            if "*" in C_gene:
                                C_gene = C_gene.split("*")[0]
                            isotype = C_gene
                            full_name = rec.contig_name + "_" + rec.identifier
                            clone = None
                            if locus in cell_contig_clones.keys():
                                if cell_name in cell_contig_clones[locus].keys():
                                    if full_name in cell_contig_clones[locus][cell_name].keys():
                                        clone = cell_contig_clones[locus][cell_name][full_name]
                        
                            string = (databasedict[locus][rec.contig_name] + 
                                "\t{}\t{}\n".format(isotype, clone))
                            output.write(string)


    def create_multiplet_database_file(self, outdir, loci, cells, multiplets):
        """Creates a tab-delimited Change-O database containing all potential
        multiplets filtered out unless --include_multiplets option is given"""
        
        changeo_string = ("CELL\tSEQUENCE_ID\tCONTIG_NAME\tLOCUS\tFUNCTIONAL\t"
                        "IN_FRAME\tSTOP\tINDELS\tV_CALL\tTOP_V_ALLELE\tD_CALL"
                        "\tTOP_D_ALLELE\tJ_CALL\tTOP_J_ALLELE\tSEQUENCE_INPUT"
                        "\tSEQUENCE_VDJ\tV_SEQ_START\tV_SEQ_LENGTH\tD_SEQ_START"
                        "\tD_SEQ_LENGTH\tJ_SEQ_START\tJ_SEQ_LENGTH\tCDR3\t"
                        "CDR3_LENGTH\tJUNCTION\tJUNCTION_LENGTH\tTPM\tC_CALL\n")

        database_file = "{}/filtered_multiplets_changeodb.tab".format(outdir)

        with open(database_file, 'w') as output:
            output.write(changeo_string)
            sorted_cell_names = sorted(multiplets)
            for cell_name in sorted_cell_names:
                cell = cells[cell_name]
                databasedict = cell.databasedict
                for locus in loci:
                    recs = cell.recombinants["BCR"][locus]
                    for rec in recs:
                        string = databasedict[locus][rec.contig_name] + "\n"
                        output.write(string)


    def make_changeo_input(self, outdir, locus, cells):
        """Creates tab-delimited database file for each locus compatible with 
        Change-O"""
        
        changeo_string = ("SEQUENCE_ID\tV_CALL\tD_CALL\tJ_CALL\tSEQUENCE_VDJ\t" 
                          "JUNCTION_LENGTH\tJUNCTION\tCELL\tISOTYPE\n")

        changeo_input = "{}/changeo_input_{}.tab".format(outdir, locus)

        with open(changeo_input, 'w') as output:
            empty = True
            for cell_name, cell in six.iteritems(cells):
                productive = cell.count_productive_recombinants(locus)

                if productive > 0:
                    if empty == True:
                        output.write(changeo_string)
                        empty = False
                    write_string = (cell.changeodict[locus].rstrip() + 
                                "\t{}\t{}\n".format(cell.name, cell.isotype))
                    output.write(write_string)
        

    def read_changeo_results(self, loci, outdir):
        """Reads Change-O database files after running DefineClones and 
        finds clone groups for each locus"""

        clones = dict()
        cell_clones = dict()
        cell_contig_clones = dict()
        cell_list = []
        contig_list = []

        for l in loci:
            clones[l] = dict()
            cell_clones[l] = dict()
            cell_contig_clones[l] = dict()
            changeo_result = "{}/changeo_input_{}_clone-pass.tab".format(
                                                                outdir, l)
            if not os.path.exists(changeo_result) \
                or not os.path.getsize(changeo_result) > 0:
                continue
            else:
                with open(changeo_result, 'r') as input_file:
                    for line in input_file:
                        if not line.startswith("SEQUENCE_ID"):
                            
                            fields = line.split("\t")
                            clone = fields[len(fields)-1].rstrip()
                            cell = fields[0].split("_TRINITY")[0]
                            contig_name = fields[0].split("{}_".format(cell))[1]
                            if not cell in cell_clones[l].keys():
                                cell_clones[l][cell] = clone
                            if not clone in clones[l].keys():
                                clones[l][clone] = dict()
                                contig_list = [contig_name]
                                clones[l][clone][cell] = []

                            if not cell in cell_contig_clones[l].keys():
                                cell_contig_clones[l][cell] = dict()

                            if not cell in clones[l][clone].keys():
                                contig_list = [contig_name]
                                clones[l][clone][cell] = []
                            
                            clones[l][clone][cell].append(contig_name)
                            cell_contig_clones[l][cell][contig_name] = clone

        return (clones, cell_clones, cell_contig_clones)


    def make_clone_igblast_input_for_concatenation(self, outdir, locus, cells, 
                        component_groups, clone_groups, cell_contig_clones):
        """Creates input file for each locus containing sequences belonging to
        a clone group shared in a component group to use as input for IgBlast 
        in order to obtain IMGT-gapped sequences for germline sequence 
        assignment by Change-O"""

        igblast_input = "{}/igblast_input_{}.fa".format(outdir, locus)
        changeo_result = "{}/changeo_input_{}_clone-pass.tab".format(outdir, locus)
        if not os.path.exists(changeo_result) \
            or not os.path.getsize(changeo_result) > 0:
            pass

        else:
            with open(changeo_result, 'r') as input_file:
                with open(igblast_input, "w") as output:
                    for line in input_file:
                        if len(line) > 0 and not line.startswith("SEQUENCE_ID"):
                            write = False
                            fields = line.split("\t")
                            name = fields[0]
                            cell = fields[0].split("_TRINITY")[0]
                            contig_name = fields[0].split("{}_".format(cell))[1]
                            seq = fields[4] + "\n"
                            header = ">" + name + "\n"
                            clone = cell_contig_clones[locus][cell][contig_name] + "_" + locus
                            i = 0
                            
                            for g in component_groups:
                                i += 1
                                top_heavy = clone_groups[i]["heavy"]
                                top_light = clone_groups[i]["light"]
                                for cell_name in g:
                                    if cell == cell_name:
                                        if clone == top_heavy or clone == top_light:
                                            output.write(header)
                                            output.write(seq)



    def make_clone_igblast_input(self, outdir, locus, cells):
        """Creates input file for each locus containing sequences belonging to 
        a clone group to use as input for IgBlast in order to obtain 
        IMGT-gapped sequences for germline sequence assignment by Change-O"""

        igblast_input = "{}/igblast_input_{}.fa".format(outdir, locus)

        changeo_result = "{}/changeo_input_{}_clone-pass.tab".format(
                                                        outdir, locus)
        if not os.path.exists(changeo_result) \
            or not os.path.getsize(changeo_result) > 0:
            pass
        
        else:
            with open(changeo_result, 'r') as input_file:
                with open(igblast_input, "w") as output:
                    for line in input_file:
                        if len(line) > 0 and not line.startswith("SEQUENCE_ID"):
                            fields = line.split("\t")
                            name = fields[0]
                            cell = fields[0].split("_TRINITY")[0]
                            contig_name = fields[0].split("{}_".format(cell))[1]
                            seq = fields[4] + "\n"
                            header = ">" + name + "\n"
                            output.write(header)
                            output.write(seq)


    def IgBlast_germline_reconstruction(self, outdir, locus):
        
        igblastn = self.get_binary('igblastn')

        # Reference data locations
        ungapped_igblast_index_location = os.path.join(self.species_dir,
                                            'igblast_dbs')
        gapped_igblast_index_location = os.path.join(self.species_dir, 
                                    'imgt_gapped_resources/igblast_dbs')
        igblast_seqtype = 'Ig'

        # IgBlast of sequences
        bracer_func.run_IgBlast_for_lineage_reconstruction(igblastn, locus, 
                outdir, ungapped_igblast_index_location, 
                gapped_igblast_index_location, igblast_seqtype, 
                self.species)
        

    def create_changeo_db(self, outdir, locus):
        """Creates Change-O database from IgBlast result files after alignment
        to imgt-gapped sequences for germline reconstruction"""

        try:
            MakeDb =  self.config.get('tool_locations', 'changeo_path') + "/MakeDb.py"
            if not os.path.is_file(MakeDb):
                MakeDb = "MakeDb.py"
        except:
            MakeDb = "MakeDb.py"

        gapped_seq_location = os.path.join(self.species_dir,
                            'imgt_gapped_resources/raw_seqs')
        ungapped_seq_location = os.path.join(self.species_dir, 'raw_seqs')
        igblast_seqtype = 'Ig'

        # Run MakeDb from Change-O toolkit
        bracer_func.run_MakeDb(MakeDb, locus, outdir, self.species, 
                        gapped_seq_location, ungapped_seq_location)
        

    def modify_changeo_db(self, outdir, locus, cells, cell_contig_clones):
        """Adds ISOTYPE, CLONE and CELL columns to the Change-O database file 
        before germline reconstruction"""
        
        input_file = "{}/igblast_{}_db-pass.tab".format(outdir, locus)
        
        if os.path.isfile(input_file):
            output_file = "{}/igblast_{}_db-modified.tab".format(outdir, locus)
            with open(input_file, "r") as input:
                with open(output_file, "w") as output:
                    for line in input:
                        if line.startswith("SEQUENCE_ID"):
                            header = line.rstrip()
                            header += "\tISOTYPE\tCLONE\tCELL\n"
                            output.write(header)
                        elif (len(line) > 0 and not line == "\n"):
                            fields = line.split("\t")
                            name = fields[0]
                            cell_name = fields[0].split("_TRINITY")[0]
                            full_contig_name = fields[0].split("{}_".format(cell_name))[1]
                            contig_name = full_contig_name.split("_IG")[0]

                            for cell in cells.values():
                                if cell.name == cell_name:
                                    recs = cell.recombinants["BCR"][locus]
                                    for rec in recs:
                                        if contig_name == rec.contig_name:
                                            C_gene = rec.C_gene
                                            if C_gene == None or C_gene == "None":
                                                C_gene = "Unknown"
                                            if "*" in C_gene:
                                                C_gene = C_gene.split("*")[0]
                                            isotype = C_gene
                                            clone = cell_contig_clones[locus][cell.name][full_contig_name]
                                            line = (line.rstrip() + 
                                                "\t{}\t{}\t{}\n".format(
                                                    isotype, clone, cell_name))
                                            output.write(line)

        # Remove old Change-O database file
        if os.path.isfile(input_file):
            if os.path.isfile(output_file):
                os.remove(input_file)


    def create_germline_sequences(self, outdir, locus, cells):
        """Runs CreateGermlines from Change-O toolkit to infer germline 
        sequences from clonal sequences"""

        try:
            CreateGermlines =  (self.config.get('tool_locations', 
                                    'changeo_path') + "/CreateGermlines.py")
            if not os.path.isfile(CreateGermlines):
                CreateGermlines = "CreateGermlines.py"
        except:
            CreateGermlines = "CreateGermlines.py"

        gapped_seq_location = os.path.join(self.species_dir,
                            'imgt_gapped_resources/raw_seqs')
        ungapped_seq_location = os.path.join(self.species_dir, 'raw_seqs')

        bracer_func.run_CreateGermlines(CreateGermlines, locus, outdir, 
                self.species, gapped_seq_location, ungapped_seq_location)


    def concatenate_lineage_databases(self, outdir, loci, component_groups, top_chain):
        """Concatenates the most commonly shared heavy and light chain
        database entries for each cell in each component group for
        lineage reconstruction"""

        output_file = "{}/concatenated_lineage_input.tab".format(outdir)
        # Read in info for each cell from locus-specific files
        cell_info = dict()
        for locus in loci:
            input_file = "{}/igblast_{}_db-modified_germ-pass.tab".format(outdir, locus)
            cell_info[locus] = dict()
            if os.path.isfile(input_file):
                with open(input_file, "r") as input:
                    for line in input:
                        if not line.startswith("SEQUENCE_ID"):
                            info = line.split("\t")
                            cell_name = line.split("_TRINITY")[0]
                            cell_info[locus][cell_name] = dict()
                            cell_info[locus][cell_name]["sequence_imgt"] = info[11]
                            cell_info[locus][cell_name]["germline_imgt"] = info[48]
                            cell_info[locus][cell_name]["V_call"] = info[7]
                            cell_info[locus][cell_name]["J_call"] = info[9]
                            cell_info[locus][cell_name]["junction_length"] = info[28]
                            cell_info[locus][cell_name]["isotype"] = info[45]


        with open(output_file, "w") as output:
        
            header = ("SEQUENCE_ID\tSEQUENCE_IMGT\tCLONE\tGERMLINE_IMGT_D_MASK" + 
                    "\tV_CALL\tJ_CALL\tJUNCTION_LENGTH\tISOTYPE\tCELL\n")
                
            output.write(header)
            clone = 0
            for g in component_groups:
                # Set clone=component_group number
                clone += 1
                light_locus = top_chain[clone]["light"].split("_")[1] 
                for cell_name in g:
                    outstring = ""
                    
                    sequence_id = cell_name
                    sequence_imgt = cell_info["H"][cell_name]["sequence_imgt"] + cell_info[light_locus][cell_name]["sequence_imgt"]
                    germline_imgt = cell_info["H"][cell_name]["germline_imgt"] + cell_info[light_locus][cell_name]["germline_imgt"]
                    V_call = cell_info["H"][cell_name]["V_call"] + "," + cell_info[light_locus][cell_name]["V_call"]
                    J_call = cell_info["H"][cell_name]["J_call"] + "," + cell_info[light_locus][cell_name]["J_call"]
                    junction_length = cell_info["H"][cell_name]["junction_length"]
                    isotype = cell_info["H"][cell_name]["isotype"]
                    outstring = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                        sequence_id, sequence_imgt, clone, 
                                        germline_imgt, V_call, J_call, 
                                        junction_length, isotype, cell_name)
                    output.write(outstring)



    def create_lineage_trees(self, outdir, species):
        """Executes R script with Alakazam commands for lineage reconstruction"""
        try:
            Rscript = self.get_binary('rscript')
        except:
            Rscript = self.get_binary('Rscript')

        dnapars = self.get_binary('dnapars')
        script = self.get_rscript_path()
        # Look for alternative species names
        if ("Hsap" or "human" or "hsap") in species:
            species = "Hsap"
        elif ("Mmus" or "mouse" or "mmus") in species:
            species = "Mmus"
        command = Rscript + " " + script + " {} {} {}".format(outdir, dnapars, species) 
        subprocess.check_call(command, shell=True)


    def create_IMGT_gapped_db(self, outdir, cells, loci, cell_contig_clones):
        """Creates IMGT-gapped Change-O database from each recombinant's
        IMGT-gapped db-string"""
        
        db_file = "{}/IMGT_gapped_db.tab".format(outdir)
        header = "CELL\tSEQUENCE_ID\tLOCUS\t" 

        header += ("TRINITY_STRING\tSEQUENCE_INPUT\tFUNCTIONAL\tIN_FRAME\tSTOP\t"
            "MUTATED_INVARIANT\tINDELS\tV_CALL\tD_CALL\tJ_CALL\tSEQUENCE_VDJ"
            "\tSEQUENCE_IMGT\tV_SEQ_START\tV_SEQ_LENGTH\tV_GERM_START_VDJ\t"
            "V_GERM_LENGTH_VDJ\tV_GERM_START_IMGT\tV_GERM_LENGTH_IMGT\t"
            "NP1_LENGTH\tD_SEQ_START\tD_SEQ_LENGTH\tD_GERM_START\t"
            "D_GERM_LENGTH\tNP2_LENGTH\tJ_SEQ_START\tJ_SEQ_LENGTH\t"
            "J_GERM_START\tJ_GERM_LENGTH\tJUNCTION_LENGTH\tJUNCTION\t"
            "FWR1_IMGT\tFWR2_IMGT\tFWR3_IMGT\tFWR4_IMGT\tCDR1_IMGT\t"
            "CDR2_IMGT\tCDR3_IMGT\tV_SCORE\tV_IDENTITY\tV_EVALUE\tV_BTOP\t"
            "J_SCORE\tJ_IDENTITY\tJ_EVALUE\tJ_BTOP")
        header += "\tISOTYPE\tC_GENE\tCLONE\n"

        with open(db_file, "w") as output:
            output.write(header)
            for cell in cells.values():
                for locus in loci:
                    recs = cell.recombinants["BCR"][locus]
                    if recs:
                        for rec in recs:
                            try:
                                db_string = rec.gapped_db_string
                            except:
                                db_string = ""
                            sequence_id = "{}_{}_{}".format(cell.name, 
                                                rec.contig_name, locus)
                            full_contig_name = "{}_{}".format(rec.contig_name, rec.identifier)
                            C_gene = rec.C_gene
                            if C_gene == None or C_gene == "None":
                                isotype = "Unknown"
                            elif "*" in C_gene:
                                isotype = C_gene.split("*")[0]
                            else:
                                isotype = None
                            try:
                                clone = cell_contig_clones[locus][cell.name][full_contig_name]
                            except:
                                clone = None
                            outline = "{}\t{}\t{}\t".format(cell.name, sequence_id, locus)
                            if len(db_string) > 0:
                                outline += db_string.strip()
                            else:
                                outline += 44*"\t"
                            outline += "\t{}\t{}\t{}\n".format(isotype, C_gene, clone)
                            output.write(outline)
        


    def count_full_length_sequences(self, loci, cells):
        """Count all full length sequences and productive full length sequences"""

        full_length_counter = dict()
        full_length_prod_counter = dict()
        productive_counter = dict()
        total_counter = dict()
        
        for l in self.loci:
            (full_length_counter[l], full_length_prod_counter[l], 
            productive_counter[l], total_counter[l]) = (0, 0, 0, 0)

        for cell_name, cell in six.iteritems(cells):
            for l in self.loci:
                full_length_count = cell.count_full_length_recombinants(l)
                full_length_prod_count = cell.count_productive_full_length_recombinants(l)
                productive = cell.count_productive_recombinants(l)
                total = cell.count_total_recombinants(l)
                if full_length_count > 0:
                    full_length_counter[l] += full_length_count
                if full_length_prod_count > 0:
                    full_length_prod_counter[l] +=full_length_prod_count
                if productive > 0:
                    productive_counter[l] += productive
                if total > 0:
                    total_counter[l] += total

        return (full_length_counter, full_length_prod_counter, 
                productive_counter, total_counter)


    def plot_full_length_sequences_proportions(self, all_dict, prod_dict, 
                                                            loci, outdir):
        D_all = all_dict
        D_prod = prod_dict
        values_all = []
        values_prod = []
        for l in self.loci:
            values_all.append(int(D_all[l]))
            values_prod.append(int(D_prod[l]))
        n_groups = len(self.loci)
        values_all = tuple(values_all)
        values_prod = tuple(values_prod)
        fig, ax = plt.subplots()
        x_ticks = tuple(self.loci)

        index = np.arange(n_groups)
        bar_width = 0.25
        opacity = 0.2

        rects1 = plt.bar(index, values_all, bar_width, color='#000000', 
                                                           label='All')
        rects2 = plt.bar(index + bar_width, values_prod, bar_width, 
                                     color='#cccccc', label='Prod')

        plt.xlabel('Locus')
        plt.ylabel('Percentage')
        plt.title('Percentage full-length of all or productive sequences')
        plt.xticks(index + bar_width, x_ticks)
        plt.legend()

        plt.savefig("{}/full_length_seqs.pdf".format(outdir))


    def make_full_length_statistics_table(self, loci, total_counter, 
        productive_counter, full_length_counter, full_length_prod_counter):
        """ Make full length statistics table"""

        header = ("##Proportion of full-length sequences of all recovered " +
                  "sequences##\n\n\t")
        all_outstring = "all\t"
        prod_outstring = "prod\t"
        all_dict = dict()
        prod_dict = dict()
        for l in self.loci:
            header = header + "{}\t".format(l)
            all_seqs = total_counter[l]
            prod_seqs = productive_counter[l]
            full_length_seqs = full_length_counter[l]
            prod_full_length_seqs = full_length_prod_counter[l]
            if prod_seqs > 0:
                prod_percent = float(prod_full_length_seqs)/int(prod_seqs)*100
                prod_percent = format(prod_percent, '.0f')
                prod_dict[l] = prod_percent
            else:
                prod_percent = "N/A"
                prod_dict[l] = "0"
            if all_seqs > 0:
                all_percent = float(full_length_seqs)/int(all_seqs)*100
                all_percent = format(all_percent, '.0f')
                all_dict[l] = all_percent
            else:
                all_percent = "N/A"
                all_dict[l] = "0"

            all_string = "{}/{} ({}%)\t".format(full_length_seqs, all_seqs, 
                                                                all_percent)
            prod_string = "{}/{} ({}%)\t".format(prod_full_length_seqs, 
                                                prod_seqs, prod_percent)
            all_outstring += all_string
            prod_outstring += prod_string

        outstring = "{header}\n{all_outstring}\n{prod_outstring}\n".format(
                    header=header, all_outstring=all_outstring, 
                    prod_outstring=prod_outstring)
        full_length_statistics_table = outstring

        return (full_length_statistics_table, all_dict, prod_dict)


    def count_isotype_usage(self, cells):
        prod_counters = defaultdict(Counter)
        cell_isotypes = []
        isotype_counter = dict()
        for cell in cells.values():
            productive = cell.count_productive_recombinants("H")
            if productive > 0:
                isotype = cell.isotype
                if isotype == None:
                    isotype = "Unknown"
                cell_isotypes.append(isotype)
        for isotype in cell_isotypes:
            if not isotype in isotype_counter:
                isotype_counter[isotype] = 1
            else:
                isotype_counter[isotype] += 1
        return (isotype_counter)


    def make_isotype_table(self, cell_recovery_count, isotype_counter):
        prod_H = cell_recovery_count["H"]
        header = ("##Isotype of cells with productive heavy chain##\n\n"
                 + "Isotype\tcells\t% of cells\n")
        outstring = ""
        for isotype, number in six.iteritems(isotype_counter):
            number = str(number)
            percent = float(number)/int(prod_H)*100
            percent = format(percent, '.2f')
            string = "{isotype}\t{number}\t{percent}\n".format(isotype=isotype, 
                                                number=number, percent=percent)
            outstring = outstring + string
        return (header, outstring)


    def plot_isotype_distributions(self, isotype_counter, outdir):
        isotypes = []
        isotype_counts = []
        D = isotype_counter
        for isotype, count in six.iteritems(isotype_counter):
            if isotype == "None":
                isotype = "Unknown"
            isotypes.append(isotype)
            isotype_counts.append(count)
        if len(isotypes) > 1:
            plt.figure()
            w = 0.85
            plt.bar(range(len(D)), D.values(), width=w, color='black', 
                                                        align='center')
            plt.xticks(range(len(D)), list(D.keys()))
            plt.xlabel("Isotype")
            plt.ylabel("Cell count")
            plt.savefig("{}/isotype_distribution.pdf".format(outdir))


    def plot_clonotype_sizes(self, outdir, cells, loci, G):

        plt.figure()
        clonotype_sizes = bracer_func.get_component_groups_sizes(
                                                cells, loci, G)
        w = 0.85
        x_range = range(1, len(clonotype_sizes) + 1)
        plt.bar(x_range, height=clonotype_sizes, width=w, color='black',
                                                        align='center')
        plt.gca().set_xticks(x_range)
        plt.xlabel("Clonotype size")
        plt.ylabel("Clonotype count")
        plt.savefig("{}/clonotype_sizes.pdf".format(outdir))

        return(x_range, clonotype_sizes)


class Tester(TracerTask):

    def __init__(self, **kwargs):
        if not kwargs:
            parser = argparse.ArgumentParser(description="Test BraCeR installation "
                              "with small dataset", parents=[self.base_parser])
            parser.add_argument('--graph_format', '-f', metavar="<GRAPH_FORMAT>", 
                                help='graphviz output format [pdf]', default='pdf')
            parser.add_argument('--no_networks', help='Skip attempts to draw '
                                'network graphs', action="store_true")
            parser.add_argument('--resume_with_existing_files', '-r',
                                help='Look for existing intermediate files '
                                'and use those instead of starting from scratch',
                                action="store_true")
            parser.add_argument('--infer_lineage', help='Construct lineage trees '
                                'for clone groups shown in clonal network', 
                                action ="store_true")
            parser.add_argument('--output', '-o', 
                                help='Directory for output data of test')
            #parser.add_argument('--no_trimming', help='Do not trim reads',
                                #action="store_true")
            
            args = parser.parse_args(sys.argv[2:])

            self.resource_dir = args.resource_dir
            self.output_dir = args.output
            self.ncores = args.ncores
            self.config_file = args.config_file
            self.graph_format = args.graph_format
            self.no_networks = args.no_networks
            self.resume = args.resume_with_existing_files
            self.infer_lineage = args.infer_lineage
            #self.no_trimming = args.no_trimming
        else:
            self.resource_dir = kwargs.get('resource_dir')
            self.output_dir = kwargs.get('output')
            self.ncores = kwargs.get('ncores')
            self.config_file = kwargs.get('config_file')
            self.graph_format = kwargs.get('graph_format', 'svg')
            self.no_networks = kwargs.get('no_networks')
            self.resume = kwargs.get('resume_with_existing_files')
            self.infer_lineage = kwargs.get('infer_lineage')
            #self.no_trimming = kwargs.get('no_trimming')

        self.trimmed_fastq1 = None
        self.trimmed_fastq2 = None
        self.keep_trimmed_reads = False

    def run(self):
        test_dir = os.path.join(base_dir, 'test_data')
        test_names = ['cell1']
        if self.output_dir:
            out_dir = os.path.join(self.output_dir, 'results')
        else:
            out_dir = os.path.join(test_dir, 'results')

        for name in test_names:
            f1 = "{}/{}_1.fastq".format(test_dir, name)
            f2 = "{}/{}_2.fastq".format(test_dir, name)

            Assembler(resource_dir=self.resource_dir, ncores=str(self.ncores), 
                      config_file=self.config_file, 
                      resume_with_existing_files=self.resume, species='Hsap', 
                      fastq1=f1, fastq2=f2, cell_name=name, output_dir=out_dir,
                      single_end=False, fragment_length=False, fragment_sd=False, 
                      loci=['H', 'K', 'L'], max_junc_len=100, 
                      no_trimming=True, trimmed_fastq1=self.trimmed_fastq1,
                      trimmed_fastq2=self.trimmed_fastq2).run()

        Summariser(resource_dir=self.resource_dir, config_file=self.config_file, 
                   use_unfiltered=False, graph_format=self.graph_format, 
                   no_networks=self.no_networks, root_dir=out_dir, 
                   loci=['H', 'K', 'L'], no_multiplets=False, IGH_networks=False, 
                   infer_lineage=self.infer_lineage, species='Hsap').run()


class Builder(TracerTask):

    """ Build Combinatorial Recombinomes and reference databases for a given species """

    def __init__(self, **kwargs):
        self.leader_padding = 20
        if not kwargs:
            parser = argparse.ArgumentParser(description="Build resources from sequences", 
                                                            parents=[self.base_parser])
            parser.add_argument('--force_overwrite', '-f', 
                                help='Force overwrite of existing resources',
                                action='store_true')
            parser.add_argument('species', metavar="<SPECIES>", help='species (eg Mmus)')
            parser.add_argument('locus_name', metavar="<LOCUS_NAME>", help='Name of locus (H, K or L)')
            parser.add_argument('N_padding', metavar="<N_PADDING>", 
                                help='Number of ambiguous N nucleotides '
                                'between V and J', type=int)
            parser.add_argument('colour', metavar="<COLOUR>", default = 'random', 
                                help='Colour for productive recombinants. '
                                'Specify as HTML (eg E41A1C) or use "random', 
                                type = self.check_colour)
            parser.add_argument('V_seqs', metavar="<V_SEQS>", 
                                help='FASTA file containing V gene sequences')
            parser.add_argument('J_seqs', metavar="<J_SEQS>", 
                                help='FASTA file containing J gene sequences')
            parser.add_argument('C_seqs', metavar="<C_SEQS>", 
                                help='FASTA file containing C gene sequence(s) '
                                'for creation of recombinomes')
            parser.add_argument('D_seqs', metavar="<D_SEQS>", nargs='?', default=False,
                                help='FASTA file containing D gene sequences (optional)')
            parser.add_argument('--C_db', metavar="<ALT_C_SEQS>", nargs='?',
                                help='Specify alternative FASTA file (if other '
                                'than the one used to make recombinomes) '
                                'containing all C gene sequences for creation of '
                                'BLAST database to correctly identify isotype (optional)',
                                default=False)
            parser.add_argument('--V_gapped', metavar="<GAPPED_V_SEQS>", nargs='?',
                                help='FASTA file containing IMGT-gapped V '
                                'reference sequences (optional, but highly recommended). '
                                'Required for accurate CDR3 detection, lineage '
                                'reconstruction and creation of IMGT-gapped '
                                'tab-delimited databases', default=False)
            parser.add_argument('--igblast_aux', metavar="<IGBLAST_AUXILARY_FILE>",
                                nargs='?', help='IgBlast auxiliary file for species. '
                                'See IgBLAST documentation for details.',
                                default=False)
            
            args = parser.parse_args(sys.argv[2:])
            resource_dir = args.resource_dir
            self.ncores = args.ncores
            self.force_overwrite = args.force_overwrite
            self.species = args.species
            self.locus_name = args.locus_name
            self.N_padding = args.N_padding
            
            self.raw_seq_files = {}
            self.raw_seq_files['V'] = args.V_seqs
            self.raw_seq_files['J'] = args.J_seqs
            self.raw_seq_files['C'] = args.C_seqs
            self.prod_colour = args.colour
            self.gapped = False
            if args.D_seqs:
                self.raw_seq_files['D'] = args.D_seqs
            if args.C_db:
                self.raw_seq_files['c'] = args.C_db
            
            self.igblast_aux = args.igblast_aux

            self.gapped_raw_seq_files = {}
            if args.V_gapped:
                self.gapped = True
                self.gapped_raw_seq_files['V'] = args.V_gapped
            
            config_file = args.config_file
            
        else:
            resource_dir = kwargs.get('resource_dir')
            self.ncores = kwargs.get('ncores')
            self.force_overwrite = kwargs.get('force_overwrite')
            self.species = kwargs.get('species')
            self.locus_name = kwargs.get('locus_name')
            self.N_padding = kwargs.get('N_padding')
            self.raw_seq_files = {}
            self.raw_seq_files['V'] = kwargs.get('V_seqs')
            self.raw_seq_files['J'] = kwargs.get('J_seqs')
            self.raw_seq_files['C'] = kwargs.get('C_seqs')
            self.prod_colour = kwargs.get('colour')
            self.gapped = False
            if kwargs.get('D_seqs'):
                self.raw_seq_files['D'] = kwargs.get('D_seqs')
            if kwargs.get('C_db'):
                self.raw_seq_files['c'] = kwargs.get('C_db')
            if kwargs.get('igblast_aux'):
                self.igblast_aux = kwargs.get('igblast_aux')
            self.gapped_raw_seq_files = {}
            if kwargs.get('V_gapped'):
                self.gapped = True
                self.gapped_raw_seq_files['V'] = kwargs.get('V_gapped')
            
            config_file = kwargs.get('config_file')

        self.config = self.read_config(config_file)

        if resource_dir is None:
            self.species_dir = os.path.join(base_dir, 'resources', self.species)
        else:
            self.species_dir = os.path.join(resource_dir, self.species)
        
        
        

    def run(self):

        self.init_dirs()
        
        self.calculate_colours(self.prod_colour)
        VDJC_files, gapped_V_file  = self.copy_raw_files()
        recombinome_fasta = self.make_recombinomes(VDJC_files)
        self.make_bowtie2_index(recombinome_fasta)
        missing_dbs = self.make_igblast_db(VDJC_files)
        for s in missing_dbs:
            print("\nIMPORTANT: there is no IgBLAST database for BCR_{segment}\n".format(
                    segment=s))
            print("Run build with {segment} segments for BCR before using bracer assemble\n".format(
                    segment=s))

        if self.gapped:
            missing_gapped_dbs = self.make_igblast_db_gapped(gapped_V_file, VDJC_files)
            for s in missing_gapped_dbs:
                print("\nIMPORTANT: there is no IgBLAST database for IMGT-gapped BCR_{}_{}\n".format(
                                            self.locus, s))
                print("Run build with --V_gapped and IMGT-gapped BCR_V segments before using bracer assemble\n")
    
    def check_colour(self, c):
        if c == 'random':
            return(c)
        else:
            try:
                if not c.startswith("#"):
                    c = "#" + c
                hex2color(c)
                return(c)
            except ValueError:
                msg = "{c} is not a valid html colour. Specify as xxXXxx".format(c=c)
                raise argparse.ArgumentTypeError(msg)
                
    def calculate_colours(self, c):
        c = c.lower()
        pal = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']
        allowed_pal = copy.copy(pal)
        colour_file = os.path.join(self.species_dir, "colours.csv")
        if os.path.exists(colour_file):
            colour_map, used_colours = io.read_colour_file(colour_file, return_used_list=True, 
                                                            receptor_name = "BCR")     
            if len(used_colours) < len(pal):
                for uc in used_colours:
                    if uc in pal:
                        allowed_pal.remove(uc)
        else:
            used_colours = None
            colour_map = {}
        if c == 'random':
            prod_colour = random.choice(allowed_pal)

        else:
            if used_colours is not None:
                if ("BCR" in colour_map and 
                   self.locus_name in colour_map["BCR"]
                   and not colour_map["BCR"][self.locus_name][0] == c):
                        if c in used_colours and c not in allowed_pal:
                            msg = "{c} already in use. Please specify a different colour.".format(c=c)
                            raise argparse.ArgumentTypeError(msg)
            
            prod_colour = c
        
        prod_rgb = hex2color(prod_colour)
        h, s, v = colorsys.rgb_to_hsv(*prod_rgb)
        
        nonprod_rgb = colorsys.hsv_to_rgb(h, s*0.5, self.adj_v(v))
        nonprod_colour = str(rgb2hex(nonprod_rgb))
        
        d1 = {self.locus_name : (prod_colour, nonprod_colour)}
        
        if "BCR" in colour_map:
            colour_map["BCR"].update(d1)
        else:
            colour_map["BCR"] = d1
        
        io.write_colour_file(colour_file, colour_map)
        
                
        
    def adj_v(self, v):
        new_v = v * 1.3
        if new_v > 1:
            new_v = 1
        return(new_v)

    def check_duplicate(self, new_path, segment=None, descriptor="Resource"):
        error_string = "{descriptor} already exists for BCR_{locus}".format(
            descriptor=descriptor, locus=self.locus_name)
        if segment:
            error_string += "_" + segment
        error_string += ". Use --force_overwrite to replace existing file."
        if os.path.isfile(new_path):
            assert self.force_overwrite, error_string

    def init_dirs(self):

        # Set up output directories
        subdirs = ['igblast_dbs', 'combinatorial_recombinomes', 'raw_seqs']

        if self.gapped:
            subdirs = subdirs + ['imgt_gapped_resources', 
                'imgt_gapped_resources/igblast_dbs', 
                'imgt_gapped_resources/raw_seqs']

        io.makeOutputDir(self.species_dir)
        for d in subdirs:
            io.makeOutputDir(os.path.join(self.species_dir, d))
            

    def copy_raw_files(self):

        """ Move user-specified files to internal resources file system """
        
        gene_segs = 'VJC'
        
        VDJC_files = {}
        
        if 'D' in self.raw_seq_files:
            gene_segs += 'D'
        if 'c' in self.raw_seq_files:
            gene_segs += 'c'

        for s in gene_segs:
            seg = s
            if seg == "c":
                seg = "C_all"
            fn = "BCR_{locus}_{seg}.fa".format(locus=self.locus_name, seg=seg)
            out_file = os.path.join(self.species_dir, 'raw_seqs', fn)
            VDJC_files[s] = out_file
            self.check_duplicate(out_file, segment=s,
                                 descriptor="Sequence File")

            # Clean up sequence names if long IMGT headers are provided
            self.clean_IMGT_fasta_files_for_IgBlast_dbs(self.raw_seq_files[s], out_file)

        # Copy IMGT-gapped V segment sequences if provided
        gapped_V_file = None
        if self.gapped:
            fn = "BCR_{locus}_V.fa".format(locus=self.locus_name)
            out_file = os.path.join(self.species_dir, 'imgt_gapped_resources/raw_seqs', fn)
            gapped_V_file = out_file
            self.check_duplicate(out_file, segment='V',
                                        descriptor="Sequence File")
            shutil.copy(self.gapped_raw_seq_files['V'], out_file)


        # Copy IgBlast auxiliary file for species if provided
        if self.igblast_aux:
            fn = "{}_gl.aux".format(self.species)
            out_file = os.path.join(self.species_dir, 
                        'imgt_gapped_resources/igblast_dbs', fn)
            if not os.path.isfile(out_file):
                shutil.copy(self.igblast_aux, out_file)

        return (VDJC_files, gapped_V_file)

     
    def load_segment_seqs(self, filename):
        seqs = {}
        with open(filename, 'rU') as fn:
            for record in SeqIO.parse(fn, 'fasta'):
                seqs[record.id] = str(record.seq)

        return seqs
        
    def make_recombinomes(self, VDJC_files):
        
        out_fasta = os.path.join(
            self.species_dir, 'combinatorial_recombinomes',
            'BCR_{locus}.fa'.format(locus=self.locus_name))

        self.check_duplicate(out_fasta, descriptor="Combinatorial recombinome")
        
        seqs = {}
        
        for s in 'VJC':
            in_file = VDJC_files[s]
            seqs[s] = self.load_segment_seqs(in_file)


        const_seq = list(seqs['C'].values())[0].upper()
        N_junction_string = "N" * self.N_padding
        N_leader_string = "N" * self.leader_padding
        
        seqs_to_write = []

        # Compile sequences to write
        if len(seqs['C']) == 1:
            for V_name, V_seq in six.iteritems(seqs['V']):
                for J_name, J_seq in six.iteritems(seqs['J']):
                    chr_name = ">chr={V_name}_{J_name}".format(J_name=J_name,
                                                           V_name=V_name)
                    seq = (N_leader_string + V_seq.lower() + N_junction_string + 
                      J_seq.lower() + const_seq)
                    seqs_to_write.append("{chr_name}\n{seq}\n".format(
                                            seq=seq, chr_name=chr_name))
        
        elif len(seqs['C']) > 1:
            for V_name, V_seq in six.iteritems(seqs['V']):
                for J_name, J_seq in six.iteritems(seqs['J']):
                    for C_name, C_seq in six.iteritems(seqs['C']):
                        chr_name = ">chr={V_name}_{J_name}_{C_name}_C_region".format(
                                        J_name=J_name, V_name=V_name, C_name=C_name)
                        seq = N_leader_string + V_seq.lower() + N_junction_string + J_seq.lower() + C_seq.upper()
                        seqs_to_write.append("{chr_name}\n{seq}\n".format(seq=seq, chr_name=chr_name))
                    
        with open(out_fasta, 'w') as f:
            for seq in seqs_to_write:
                f.write(seq)

        return out_fasta

    def make_bowtie2_index(self, recombinome_fasta):
        
        bowtie2_build = self.get_binary('bowtie2-build')
        index_base = os.path.join(
            self.species_dir, 'combinatorial_recombinomes',
            'BCR_{locus}'.format(locus=self.locus_name))

        self.check_duplicate(index_base + ".1.bt2", descriptor="Bowtie2 index")
        
        command = [bowtie2_build, '-q', recombinome_fasta, index_base]
        try:
            subprocess.check_call(command)
        except subprocess.CalledProcessError:
            print("bowtie2-build failed")
    
    def clean_IMGT_fasta_files_for_IgBlast_dbs(self, in_file, out_file):
        """Clean IMGT germline fasta files for IgBLAST database build.
        This function is modified from 
        https://bitbucket.org/kleinstein/immcantation/src/tip/scripts/clean_imgtdb.py"""
        # Load sequences into memory and process them
        name_set = set()
        seq_list = list()
        for rec in SeqIO.parse(in_file, 'fasta'):
            if '|' in rec.description:
                name = rec.description.split('|')[1]
            else:
                name = rec.description
            if name not in name_set:
                name_set.add(name)
                seq = SeqRecord(rec.seq.ungap('.').upper(), id=name, name=name, description=name)
                seq_list.append(seq)

        # Overwrite file
        with open(out_file, 'w') as out_handle:
            writer = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
            writer.write_file(seq_list)


    def make_igblast_db(self, VDJC_files):
        
        igblast_dir = os.path.join(self.species_dir, 'igblast_dbs')
        
        makeblastdb = self.get_binary('makeblastdb')
        missing_dbs = []
      
        gene_segs = 'VDJ'

        # Use alternative C sequence file if provided
        if 'c' in self.raw_seq_files:
            gene_segs += 'c'
	
	# If not alternative C sequence file is provided, use the same as for making recombinomes
        else:
            gene_segs += 'C'

        for s in gene_segs:
            seg = s
            if s == "c":
                seg = "C_all"
            fn = "BCR_{segment}.fa".format(segment=seg)
            fasta_file = os.path.join(igblast_dir, fn)

            # Create file if it doesn't already exist
            open(fasta_file, 'a').close()

            #pdb.set_trace()
            if s in VDJC_files:
                with open(fasta_file) as e:
                    existing_seqs = SeqIO.to_dict(SeqIO.parse(e, "fasta"))
                with open(VDJC_files[s]) as n:
                    new_seqs = SeqIO.to_dict(SeqIO.parse(n, "fasta"))
                
                non_overwritten_seqs = []
                
                for seq_name, seq in six.iteritems(new_seqs):
                    if seq_name in existing_seqs:
                        if not self.force_overwrite:
                            non_overwritten_seqs.append(seq_name)
                        else:
                            existing_seqs.update({seq_name: seq})
                    else:
                        existing_seqs.update({seq_name: seq})
                with open(fasta_file, 'w') as f:
                    SeqIO.write(existing_seqs.values(), f, "fasta")
                
                if len(existing_seqs) == 0:
                    missing_dbs.append(s)
            
                if len(non_overwritten_seqs) > 0:
                    print('The follwing IgBLAST DB sequences for '
                          'BCR_{segment} already found in {file}.'.format(
                                            segment=seg, file=fasta_file))
                    print('These sequences were not overwritten. '
                          'Use --force_overwrite to replace with new ones')
                    for seq in non_overwritten_seqs:
                        print(seq)
                
                command = [makeblastdb, '-parse_seqids', '-dbtype', 'nucl',
                           '-in', fasta_file]
                try:
                    subprocess.check_call(command)
                except subprocess.CalledProcessError:
                    print("makeblastdb failed for BCR_{segment}".format(
                                                        segment=seg))

        return missing_dbs

    def make_igblast_db_gapped(self, gapped_V_file, VDJC_files):
        """Create locus-specific IgBlast databases from IMGT-gapped sequences"""
        
        igblast_dir = os.path.join(self.species_dir, 'imgt_gapped_resources/igblast_dbs')
        makeblastdb = self.get_binary('makeblastdb')
        missing_gapped_dbs = []
        
        gene_segs = "VDJ"
        for s in gene_segs:
            if s=="V":
                # Clean up IMGT-gapped file to remove gaps
                cleaned_file = os.path.join(self.species_dir, 
                        'imgt_gapped_resources/raw_seqs/BCR_{}_V_cleaned.fa'.format(
                                                                self.locus_name))
                self.clean_IMGT_fasta_files_for_IgBlast_dbs(gapped_V_file, cleaned_file)
            elif s in VDJC_files:
                cleaned_file = VDJC_files[s]

            fn = "BCR_{}_{}.fa".format(self.locus_name, s)
            fasta_file = os.path.join(igblast_dir, fn)
            # Create file if it doesn't already exist
            open(fasta_file, 'a').close()

            with open(fasta_file) as e:
                existing_seqs = SeqIO.to_dict(SeqIO.parse(e, "fasta"))
            with open(cleaned_file) as n:
                new_seqs = SeqIO.to_dict(SeqIO.parse(n, "fasta"))
            non_overwritten_seqs = []
            for seq_name, seq in six.iteritems(new_seqs):
                if seq_name in existing_seqs:
                    if not self.force_overwrite:
                        non_overwritten_seqs.append(seq_name)
                    else:
                        existing_seqs.update({seq_name: seq})
                else:
                    existing_seqs.update({seq_name: seq})
            with open(fasta_file, 'w') as f:
                SeqIO.write(existing_seqs.values(), f, "fasta")

            if len(existing_seqs) == 0:
                missing_gapped_dbs.append(s)

            if len(non_overwritten_seqs) > 0:
                print('The follwing IgBLAST DB sequences for '
                    'BCR_{locus}_{s} already found in {file}.'.format(
                                locus=self.locus, s=s, file=fasta_file))
                print('These sequences were not overwritten. '
                    'Use --force_overwrite to replace with new ones')
                for seq in non_overwritten_seqs:
                    print(seq)

            command = [makeblastdb, '-parse_seqids', '-dbtype', 'nucl',
                                                    '-in', fasta_file]

            try:
                subprocess.check_call(command)
            except subprocess.CalledProcessError:
                print("makeblastdb failed for BCR_{}_{} (IMGT-gapped)".format(self.locus, s))

        return missing_gapped_dbs
