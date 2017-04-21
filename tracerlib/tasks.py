from __future__ import print_function

import csv

import six
import matplotlib as mpl

from tracerlib import io, core
from tracerlib.io import check_binary

mpl.use('pdf')
import re
import seaborn as sns
from matplotlib import pyplot as plt
from tracerlib import base_dir
from tracerlib import tracer_func
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
import itertools
import pdb
import numpy as np
from numpy import percentile, array
from matplotlib.colors import hex2color, rgb2hex
import random
import copy
import colorsys

class TracerTask(object):

    base_parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    base_parser.add_argument('--ncores', '-p', metavar="<CORES>", help='number of processor cores to use', type=int,
                             default=1)
    base_parser.add_argument('--config_file', '-c', metavar="<CONFIG_FILE>", help='config file to use',
                             default='~/.tracerrc')

    assembler = "trinity"
    config = None
    def run(self):
        pass

    def get_binary(self, name):
        tool_key = name.lower() + '_path'
        user_path = None
        if self.config.has_option('tool_locations', tool_key):
            user_path = self.resolve_relative_path(self.config.get('tool_locations', tool_key))
        return check_binary(name, user_path)

    def read_config(self, config_file):
        # Read config file
        if not config_file:
            config_file = '~/.tracerrc'
        config_file = os.path.expanduser(config_file)
        if not os.path.isfile(config_file):
            print("Config file not found at ~/.tracerrc. Using default tracer.conf in repo...")
            config_file = os.path.join(base_dir, 'tracer.conf')
        tracer_func.check_config_file(config_file)
        config = ConfigParser()
        config.read(config_file)
        return config

    def resolve_relative_path(self, path):
        if not path.startswith("/"):
            base_directory = os.path.abspath(os.path.dirname(__file__))
            full_path = os.path.normpath("/{}/../{}".format(base_directory, path))
        else:
            full_path = path
        return full_path

    def print_cell_summary(self, cell, output_file, receptor_name, loci):
        out_file = open(output_file, 'w')
        out_file.write('------------------\n{name}\n------------------\n'.format(name=cell.name))
        
        # summarise the productive/total recombinants
        for l in loci:
            out_file.write('{receptor}_{locus} recombinants: {summary}\n'.format(receptor=receptor_name, locus=l,
                                                            summary=cell.summarise_productivity(receptor_name, l)))
        
        out_file.write('\n\n')
        
        for l in loci:
            out_file.write("#{receptor}_{locus}#\n".format(receptor=receptor_name, locus=l))
            rs = cell.recombinants[receptor_name][l]
            if rs is None:
                out_file.write("No {receptor}_{locus} recombinants found\n\n".format(receptor=receptor_name, locus=l))
            else:
                for r in rs:
                    out_file.write(r.get_summary())
                    out_file.write("\n\n")
        
        out_file.close()

    def die_with_empty_cell(self, cell_name, output_dir, species, assembler):
        print("##No recombinants found##")
        cell = core.Cell(cell_name, None, is_empty=True, species=species, receptor=self.receptor_name, loci=self.loci)
        
        # Save cell in a pickle

        unfiltered_summary = "{output_dir}/unfiltered_{receptor}_seqs/unfiltered_{receptor}s.txt".format(
                                                  output_dir=self.output_dir, receptor=self.receptor_name)
        pickle_file = "{output_dir}/unfiltered_{receptor}_seqs/{cell_name}.pkl".format(output_dir=self.output_dir,
                                                          cell_name=cell.name, receptor=self.receptor_name)
        filtered_summary = "{output_dir}/filtered_{receptor}_seqs/filtered_{receptor}s.txt".format(
                                                         output_dir=self.output_dir, receptor=self.receptor_name)
                                                                            
        filtered_pickle = "{output_dir}/filtered_{receptor}_seqs/{cell_name}.pkl".format(output_dir=self.output_dir,
                                                                  cell_name=cell.name, receptor=self.receptor_name)


        self.print_cell_summary(cell, unfiltered_summary, self.receptor_name, self.loci)
                                                                        

        with open(pickle_file, 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)

        cell.filter_recombinants()
        self.print_cell_summary(cell, filtered_summary, self.receptor_name, self.loci)
                                                                            
        with open(filtered_pickle, 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)
        
        exit(0)
    
    def get_resources_root(self, species):
        resources_dir = os.path.join(base_dir, 'resources')
        resources_root = os.path.join(resources_dir, species)
        return(resources_root)
    
    def get_available_species(self):
        resources_dir = os.path.join(base_dir, 'resources')
        species_dirs = next(os.walk(resources_dir))[1]
        return(species_dirs)


class Assembler(TracerTask):

    def __init__(self, **kwargs):
        if not kwargs:
            
            # get list of all available species in resources
            
            parser = argparse.ArgumentParser(
                description="Reconstruct TCR sequences from RNAseq reads for a single cell", 
                parents=[self.base_parser], formatter_class=argparse.ArgumentDefaultsHelpFormatter)
            parser.add_argument('--resume_with_existing_files', '-r',
                                help='look for existing intermediate files and use those instead of starting from scratch',
                                action="store_true")
            parser.add_argument('--assembled_file',
                                help='Fasta file containing already assembled sequences for the cell. Providing this \
                                file skips the alignment and assembly steps', default=False)
            parser.add_argument('--species', '-s',
                                help='Species to use for reconstruction', choices=self.get_available_species(), default='Mmus')
            parser.add_argument('--receptor_name',
                                help="Name of receptor to reconstruct", default='TCR')
            parser.add_argument('--loci',
                                help="Space-separated list of loci to reconstruct for receptor", 
                                default=['A','B'], nargs = '+')
            parser.add_argument('--seq_method', '-m',
                                help='Method for constructing sequence to assess productivity, \
                                quantify expression and for output reporting. See README for details.',
                                choices=['imgt', 'assembly'], default='imgt')
            parser.add_argument('--single_end', help='set this if your sequencing data are single-end reads',
                                action="store_true")
            parser.add_argument('--fragment_length',
                                help='Estimated average fragment length in the sequencing library.'
                                     ' Used for Kallisto quantification. REQUIRED for single-end data.',
                                default=False)
            parser.add_argument('--fragment_sd',
                                help='Estimated standard deviation of average fragment length in the sequencing library.'
                                     ' Used for Kallisto quantification. REQUIRED for single-end data.',
                                default=False)
            parser.add_argument('--max_junc_len',
                                help="Maximum permitted length of junction string in recombinant identifier. "
                                     "Used to filter out artefacts. May need to be longer for TCRdelta.",
                                default=50)

            parser.add_argument('fastq1', metavar="<FASTQ1>", help='first fastq file', nargs='?')
            parser.add_argument('fastq2', metavar="<FASTQ2>", help='second fastq file', nargs='?')
            parser.add_argument('cell_name', metavar="<CELL_NAME>", help='name of cell for file labels')
            parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                help='directory for output as <output_dir>/<cell_name>')

            args = parser.parse_args(sys.argv[2:])

            self.cell_name = args.cell_name
            self.fastq1 = args.fastq1
            self.single_end = args.single_end
            self.fastq2 = args.fastq2
            self.ncores = str(args.ncores)
            self.assembled_file = args.assembled_file
            self.species = args.species
            self.seq_method = args.seq_method
            self.resume_with_existing_files = args.resume_with_existing_files
            self.fragment_length = args.fragment_length
            self.fragment_sd = args.fragment_sd
            self.output_dir = args.output_dir
            self.receptor_name = args.receptor_name
            self.loci = args.loci
            self.max_junc_len = args.max_junc_len
            config_file = args.config_file
            

        else:
            self.cell_name = kwargs.get('cell_name')
            self.fastq1 = kwargs.get('fastq1')
            self.fastq2 = kwargs.get('fastq2')
            self.ncores = kwargs.get('ncores')
            self.loci = kwargs.get('loci')
            self.max_junc_len = kwargs.get('max_junc_len')
            config_file = kwargs.get('config_file')
           

        self.config = self.read_config(config_file)
        
        self.assembler = "trinity"
        if not self.assembled_file:
            self.assembled_file = None

        #pdb.set_trace()

        # Check the fastq config is correct
        if not self.assembled_file:
            if not self.single_end:
                assert self.fastq2, "Only one fastq file specified. Either set --single_end or provide second fastq."
            else:
                self.fastq2 = None
                if self.fastq2:
                    print("Two fastq files given with --single-end option. Ignoring second file.")
                assert self.fragment_length and self.fragment_sd, \
                    'Must specify estimated average fragment length (--fragment_length)' \
                    ' and standard deviation (--fragment_sd) for use with single-end data'
                assert self.fragment_length, \
                    'Must specify estimated average fragment length (--fragment_length) for use with single-end data'
                assert self.fragment_sd, \
                    'Must specify estimated fragment length standard deviation (--fragment_sd) for use with single-end data'

        # Check FASTQ files exist
        if not self.assembled_file:
            if not os.path.isfile(self.fastq1):
                raise OSError('2', 'FASTQ file not found', self.fastq1)
            if not self.single_end and self.fastq2:
                if not os.path.isfile(self.fastq2):
                    raise OSError('2', 'FASTQ file not found', self.fastq2)

        # Check that fasta file containing already assembled sequences exists (if provided)
        if self.assembled_file:
            if not os.path.isfile(self.assembled_file):
                raise OSError('2', 'Fasta file containing assembled sequences not found. Please provide fasta \
                             file or run TraCeR without --assembled_file option.', self.assembled_file)
        
        
        
    def run(self, **kwargs):

        # Set-up output directories
        root_output_dir = os.path.abspath(self.output_dir)
        io.makeOutputDir(root_output_dir)
        self.output_dir = root_output_dir + "/" + self.cell_name

        io.makeOutputDir(self.output_dir)

        data_dirs = ['aligned_reads', 'Trinity_output', 'IgBLAST_output', 'BLAST_output',
                     'unfiltered_{receptor}_seqs'.format(receptor=self.receptor_name),'expression_quantification', 
                     'filtered_{receptor}_seqs'.format(receptor=self.receptor_name)]

        for d in data_dirs:
            io.makeOutputDir("{}/{}".format(self.output_dir, d))

        # Perform TraCeR's core functions
        if not self.assembled_file:
            self.align()
            self.de_novo_assemble()

        self.blast()
        cell = self.ig_blast()

        if not self.assembled_file:
            self.quantify(cell)
        
        
        unfiltered_fasta_filename = "{output_dir}/unfiltered_{receptor}_seqs/{cell_name}_{receptor}seqs.fa".format(output_dir=self.output_dir,
                                                                                        cell_name=self.cell_name, receptor=self.receptor_name)
        unfiltered_cell_summary_file =  "{output_dir}/unfiltered_{receptor}_seqs/unfiltered_{receptor}s.txt".format(output_dir=self.output_dir, 
                                                                                                                      receptor=self.receptor_name)

        unfiltered_pickle = "{output_dir}/unfiltered_{receptor}_seqs/{cell_name}.pkl".format(output_dir=self.output_dir, cell_name=cell.name,
                                                                                                                     receptor=self.receptor_name)


        filtered_fasta_filename = "{output_dir}/filtered_{receptor}_seqs/{cell_name}_{receptor}seqs.fa".format(output_dir=self.output_dir,
                                                                                          cell_name=self.cell_name, receptor=self.receptor_name)
        filtered_cell_summary_file =  "{output_dir}/filtered_{receptor}_seqs/filtered_{receptor}s.txt".format(output_dir=self.output_dir,
                                                                                                                    receptor=self.receptor_name)
        filtered_pickle = "{output_dir}/filtered_{receptor}_seqs/{cell_name}.pkl".format(output_dir=self.output_dir, cell_name=cell.name,
                                                                                                                    receptor=self.receptor_name)

                                                                        
        unfiltered_fasta_file = open(unfiltered_fasta_filename, 'w')
        unfiltered_fasta_file.write(cell.get_fasta_string())
        unfiltered_fasta_file.close()

        self.print_cell_summary(cell, unfiltered_cell_summary_file, self.receptor_name, self.loci)
        
        # Assign isotype and bgcolor (for B cells)
        if self.receptor_name == "BCR":
            ranked_recs = cell.rank_recombinants()
            isotype = cell.determine_isotype(ranked_recs)
            bgcolor = cell.assign_bgcolor(isotype)
            cell.bgcolor = bgcolor
            cell.isotype = isotype

        # Save cell in a pickle
        with open(unfiltered_pickle, 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)


        # Filter recombinants
        print("##Filtering by read count##")
        cell.filter_recombinants()
        filtered_fasta_file = open(filtered_fasta_filename, 'w')
        filtered_fasta_file.write(cell.get_fasta_string())
        filtered_fasta_file.close()
        self.print_cell_summary(cell, filtered_cell_summary_file, self.receptor_name, self.loci)

        # Save cell in a pickle
        
        with open(filtered_pickle, 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)


        
    def get_index_location(self, name):
        location = os.path.join(base_dir, 'resources', self.species, name)

        return location

    def align(self):
        bowtie2 = self.get_binary('bowtie2')

        synthetic_genome_path = self.get_index_location('combinatorial_recombinomes')
        # Align with bowtie
        tracer_func.bowtie2_alignment(
            bowtie2, self.ncores, self.receptor_name, self.loci, self.output_dir, self.cell_name, 
            synthetic_genome_path, self.fastq1, self.fastq2, self.resume_with_existing_files, self.single_end)
        print()

    def de_novo_assemble(self):

        trinity = self.get_binary('trinity')

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
            trinity_grid_conf = self.resolve_relative_path(self.config.get('trinity_options', 'trinity_grid_conf'))
        else:
            trinity_grid_conf = False

        # De novo assembly with trinity
        trinity_JM = self.config.get('trinity_options', 'max_jellyfish_memory')
        trinity_version = self.config.get('trinity_options', 'trinity_version')
        successful_files = tracer_func.assemble_with_trinity(
            trinity, self.receptor_name, self.loci, self.output_dir, self.cell_name, self.ncores, trinity_grid_conf, 
            trinity_JM, trinity_version, self.resume_with_existing_files, self.single_end, self.species)
        if len(successful_files) == 0:
            print("No successful Trinity assemblies")
            self.die_with_empty_cell(self.cell_name, self.output_dir, self.species, self.assembler)

        print()


    def ig_blast(self):
        igblastn = self.get_binary('igblastn')

        # Reference data locations
        igblast_index_location = self.get_index_location('igblast_dbs')
        imgt_seq_location = self.get_index_location('raw_seqs')

        igblast_seqtype = self.config.get('IgBlast_options', 'igblast_seqtype')


        # IgBlast of assembled contigs
        tracer_func.run_IgBlast(igblastn, self.receptor_name, self.loci, self.output_dir, self.cell_name, 
           igblast_index_location, igblast_seqtype, self.species, self.resume_with_existing_files, self.assembled_file)
        print()
        
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cell = io.parse_IgBLAST(self.receptor_name, self.loci, self.output_dir, self.cell_name, imgt_seq_location, 
                                    self.species, self.seq_method, self.assembler, self.assembled_file, self.max_junc_len)
            if cell.is_empty:
                self.die_with_empty_cell(self.cell_name, self.output_dir, self.species, self.assembler)

        return cell

    def blast(self):
        blastn = self.get_binary('blastn')


        # Reference data locations
        blast_index_location = self.get_index_location('igblast_dbs')
        imgt_seq_location = self.get_index_location('raw_seqs')

        # Blast of assembled contigs
        tracer_func.run_Blast(blastn, self.receptor_name, self.loci, self.output_dir, self.cell_name, blast_index_location,
                                self.species, self.resume_with_existing_files, self.assembler, self.assembled_file)
        print()


        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            isotype = io.parse_BLAST(self.receptor_name, self.loci, self.output_dir, self.cell_name, self.species, self.assembler)
        return isotype


    def quantify(self, cell):
        kallisto = self.get_binary('kallisto')
        
        if not self.config.has_option('kallisto_transcriptomes', self.species):
            raise OSError("No transcriptome reference specified for {species}. Please specify location in config file."
                          .format(species = self.species))
        else:
            kallisto_base_transcriptome = self.resolve_relative_path(self.config.get('kallisto_transcriptomes',
                                                                                 self.species))

        #Quantification with kallisto
        tracer_func.quantify_with_kallisto(
            kallisto, cell, self.output_dir, self.cell_name, kallisto_base_transcriptome, self.fastq1, self.fastq2,
            self.ncores, self.resume_with_existing_files, self.single_end, self.fragment_length, self.fragment_sd, self.receptor_name)
        print()

        counts = tracer_func.load_kallisto_counts("{}/expression_quantification/abundance.tsv".format(self.output_dir))
        
        for receptor, locus_dict in six.iteritems(cell.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None:
                    for rec in recombinants:
                        tpm = counts[receptor][locus][rec.contig_name]
                        rec.TPM = tpm
                    

class Summariser(TracerTask):

    def __init__(self, **kwargs):

        if not kwargs:
            parser = argparse.ArgumentParser(description="Summarise set of cells with reconstructed TCR/BCR sequences",
                                             parents=[self.base_parser])
            parser.add_argument('--species', '-s',
                                help='Species to use for reconstruction',
                                choices=self.get_available_species(), default='Mmus')
            parser.add_argument('--receptor_name',
                                help="Name of receptor to summarise", default='TCR')
            parser.add_argument('--loci',
                                help="Space-separated list of loci to summarise for receptor", 
                                default=['A','B'], nargs = '+')
            parser.add_argument('--use_unfiltered', '-u', help='use unfiltered recombinants', action="store_true")
            parser.add_argument('--keep_invariant', '-i', help='ignore invariant cells when constructing networks',
                                action="store_true")
            parser.add_argument('--graph_format', '-f', metavar="<GRAPH_FORMAT>", help='graphviz output format [pdf]',
                                default='pdf')
            parser.add_argument('--no_networks', help='skip attempts to draw network graphs', action = "store_true")
            parser.add_argument('--IGH_networks', help='base network graphs only on IGH chains', action = "store_true")
            parser.add_argument('dir', metavar="<DIR>",
                                help='directory containing subdirectories for each cell to be summarised')
            args = parser.parse_args(sys.argv[2:])

            self.root_dir = os.path.abspath(args.dir)
            self.graph_format = args.graph_format
            self.keep_invariant = args.keep_invariant
            self.use_unfiltered = args.use_unfiltered
            self.draw_graphs = not args.no_networks
            self.IGH_networks = args.IGH_networks
            self.receptor_name = args.receptor_name
            self.loci = args.loci
            self.species = args.species
            config_file = args.config_file
        else:
            self.use_unfiltered = kwargs.get('use_unfiltered')
            self.root_dir = os.path.abspath(kwargs.get('root_dir'))
            self.draw_graphs = not (kwargs.get('no_networks'))
            self.IGH_networks = kwargs.get('IGH_networks')
            self.graph_format = kwargs.get('graph_format')
            self.keep_invariant = kwargs.get('keep_invariant')
            self.receptor_name = kwargs.get('receptor_name')
            self.loci = kwargs.get('loci')
            self.species = kwargs.get('species')
            config_file = kwargs.get('config_file')

        # Read config file
        self.config = self.read_config(config_file)
        
        self.species_dir = self.get_resources_root(self.species)
        
        invariant_cells = self.resolve_relative_path(os.path.join('resources', self.species,
                                                                         'invariant_cells.json'))
        if os.path.isfile(invariant_cells):
            self.invariant_cells = io.parse_invariant_cells(invariant_cells)
        else:
            self.invariant_cells = None


        
    def run(self):

        if self.draw_graphs:
            dot = self.resolve_relative_path(self.config.get('tool_locations', 'dot_path'))
            neato = self.resolve_relative_path(self.config.get('tool_locations', 'neato_path'))

            # check that executables from config file can be used
            not_executable = []
            for name, x in six.iteritems({"dot": dot, "neato": neato}):
                if not io.is_exe(x):
                    not_executable.append((name, x))
            if len(not_executable) > 0:
                print()
                print("Could not execute the following required tools. Check your configuration file.")
                for t in not_executable:
                    print( t[0], t[1])
                print()
                exit(1)
        else:
            dot = ""
            neato = ""

        cells = {}
        empty_cells = []
        subdirectories = next(os.walk(self.root_dir))[1]

        if self.use_unfiltered:
            pkl_dir = "unfiltered_{}_seqs".format(self.receptor_name)
            outdir = "{}/unfiltered_{}_summary".format(self.root_dir, self.receptor_name)

        else:
            pkl_dir = "filtered_{}_seqs".format(self.receptor_name)
            outdir = "{}/filtered_{}_summary".format(self.root_dir, self.receptor_name)

        io.makeOutputDir(outdir)

        

        outfile = open("{}/{}_summary.txt".format(outdir, self.receptor_name), 'w')
        length_filename_root = "{}/reconstructed_lengths_{}".format(outdir, self.receptor_name)
        isotype_filename_root = "{}/isotypes_{}".format(outdir, self.receptor_name)
        cdr3_filename_root = "{}/cdr3_lengths_{}".format(outdir, self.receptor_name)

        for d in subdirectories:
            cell_pkl = "{root_dir}/{d}/{pkl_dir}/{d}.pkl".format(pkl_dir=pkl_dir, d=d, root_dir=self.root_dir)
            if os.path.isfile(cell_pkl):
                with open(cell_pkl, 'rb') as pkl:
                    cl = pickle.load(pkl)
                cells[d] = cl
                if cl.is_empty or cl.missing_loci_of_interest(self.receptor_name, self.loci):
                    empty_cells.append(d)
               
                
                
        
        cell_recovery_count = dict()
        # count cells with productive chains for each locus and for each possible pair
        for l in self.loci:
            cell_recovery_count[l] = 0

        
        possible_pairs = ["".join(x) for x in itertools.combinations(self.loci, 2)]
        
        for p in possible_pairs:
            cell_recovery_count[p] = 0
        
        for cell_name, cell in six.iteritems(cells):
            prod_counts = dict()
            for l in self.loci:
                prod_counts[l] = cell.count_productive_recombinants(self.receptor_name, l)
                if prod_counts[l] > 0:
                    cell_recovery_count[l] += 1
         
            for p in possible_pairs:
                if prod_counts[p[0]] > 0 and prod_counts[p[1]] > 0:
                    cell_recovery_count[p] += 1

        total_cells = len(cells)
       

        for l in self.loci:
            count = cell_recovery_count[l]
            pc = round((count/float(total_cells))*100, 1)
            outfile.write("{receptor}_{locus} reconstruction:\t{count} / {total} ({pc}%)\n".format(
                                                                                    receptor=self.receptor_name,
                                                                                    locus=l, count=count,
                                                                                    total=total_cells, pc=pc))
        outfile.write("\n")
        
        for p in possible_pairs:
            count = cell_recovery_count[p]
            pc = round((count/float(total_cells))*100, 1)
            if pc == "KL":
                outfile.write("Productive reconstruction of K and L:\t{count} / {total} ({pc}%)\n".format(
                                                                                    p=p,
                                                                                    count=count,
                                                                                    total=total_cells, pc=pc))
            else:
                outfile.write("Paired {p} productive reconstruction:\t{count} / {total} ({pc}%)\n".format(
                                                                                    p=p,
                                                                                    count=count,
                                                                                    total=total_cells, pc=pc))
        outfile.write("\n")
        
        
        all_counters = defaultdict(Counter)
        prod_counters = defaultdict(Counter)
        if self.receptor_name == "BCR":
            isotype_counters = defaultdict(Counter)
            possible_isotypes = ["IGHM", "IGHG1", "IGHG2A", "IGHG2B", "IGHG2C", "IGHG3", "IGHG4", "IGHA", "IGHA1", "IGHA2", "IGHE", "IGHD"]
        
        for cell in cells.values():
            for l in self.loci:
                all_counters[l].update({cell.count_total_recombinants(self.receptor_name, l): 1})
                prod_counters[l].update({cell.count_productive_recombinants(self.receptor_name, l): 1})
        
        all_recombinant_counts = []
        
        for locus in all_counters:
            all_recombinant_counts = all_recombinant_counts + \
                                     list(all_counters[locus].keys())
        max_recombinant_count = max(all_recombinant_counts)
        
        #max_recombinant_count = max(list(counters['all_alpha'].keys()) + list(counters['all_beta'].keys()))
        table_header = ['', '0 recombinants', '1 recombinant', '2 recombinants']
        recomb_range = range(0, 3)
        if max_recombinant_count > 2:
            extra_header = [str(x) + " recombinants" for x in range(3, max_recombinant_count + 1)]
            table_header = table_header + extra_header
            recomb_range = range(0, max_recombinant_count + 1)

        t = PrettyTable(table_header)
        t.padding_width = 1
        t.align = "l"

        
        #make all recombinant table
        for counter_name in ['all_counters', 'prod_counters']:
            counter_type = counter_name.split("_")[0]
            counter_set = eval(counter_name)
            for l in self.loci:
                counter = counter_set[l]
                count_array = [counter[x] for x in recomb_range]
                total_with_at_least_one = sum(count_array[1:])
                if total_with_at_least_one > 0:
                    percentages = [''] + [" (" + str(round((float(x) / total_with_at_least_one) * 100)) + "%)" for x in
                                          count_array[1:]]
                else:
                    percentages = [''] + [" (N/A%)" for x in count_array[1:]]
                row = []
                for i in recomb_range:
                    row.append(str(count_array[i]) + percentages[i])
                label = '{} {}'.format(counter_type, l)
                t.add_row([label] + row)
        
        
        
        outfile.write(t.get_string())
        outfile.write("\n")

        # If using unfiltered, name cells with more than two recombinants#
        if self.use_unfiltered:
            outfile.write("\n\n#Cells with more than two recombinants for a locus#\n")
            found_multi = False
            for cell in cells.values():
                if cell.has_excess_recombinants:
                    outfile.write("###{}###\n".format(cell.name))
                    for l in self.loci:
                        count = cell.count_total_recombinants(self.receptor_name, l)
                        outfile.write("{receptor}_{l}:\t{count}\n".format(receptor=self.receptor_name, l=l, count=count))
                    outfile.write("\n")
                    found_multi = True
            if not found_multi:
                outfile.write("None\n\n")
      
        # Make full length statistics table and plot proportion of sequences that are full-length
        (full_length_counter, full_length_prod_counter, productive_counter, total_counter) = self.count_full_length_sequences(self.loci, cells, self.receptor_name)
        (full_length_statistics_table, all_dict, prod_dict) = self.make_full_length_statistics_table(self.loci, total_counter, productive_counter, full_length_counter, full_length_prod_counter)
        outfile.write(full_length_statistics_table)
        outfile.write("\n")
        self.plot_full_length_sequences_proportions(all_dict, prod_dict, self.loci, outdir)
         
        #Report cells with two productive chains from same locus
        outstring = self.two_productive_chains_per_locus(self.loci, cells, self.receptor_name)
        if len(outstring) > 0:
            outfile.write(outstring)
      
        # B CELL SPECIFIC TASKS

        if self.receptor_name == "BCR":
            #Report cells with productive kappa and lambda chain
            outstring = self.report_kappa_lambda_cells(self.loci, cells, self.receptor_name)
            if len(outstring) > 0:
                outfile.write(outstring)

            # Make initial clonal assignments for B cells using ChangeO DefineClones bygroup
            for locus in self.loci:
                # Create input file for ChangeO
                self.make_changeo_input(outdir, locus, self.receptor_name, cells)
                # Run ChangeO
                changeo = self.get_binary('changeo')
                tracer_func.run_changeo(changeo, locus, outdir, self.species)
                print()
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")

            # Read ChangeO result files and define clone groups for each locus
            (clones, cell_clones, cell_contig_clones) = self.read_changeo_results(self.loci, outdir)
            print(cell_contig_clones)

            print("#####################   Cell contig clones   ############################\n")
            for locus in ["H", "K", "L"]:
                print("\n---{}---".format(locus))
                for cell_name, contig_info in six.iteritems(cell_contig_clones[locus]):
                    print(cell_name)
                    
                    for contig_name, clone in six.iteritems(cell_contig_clones[locus][cell_name]):
                        print(contig_name, clone)

            # Get H clone groups consisting of 2 or more cells
            multiple_clones_H = dict()
            for clone, cells_info in six.iteritems(clones["H"]):
                if len(cells_info.keys()) > 1:
                    multiple_clones_H[clone] = cells_info

            cells_with_clonal_H = []
            print("#####################   H clone groups (> 1 seq)   ############################\n")
            for clone, cells_info in six.iteritems(multiple_clones_H):
                print("\n---{}---".format(clone))
                for cell_name, contig_name in six.iteritems(multiple_clones_H[clone]):
                    print(cell_name, contig_name)
                    if not cell_name in cells_with_clonal_H:
                        cells_with_clonal_H.append(cell_name)
            print(cells_with_clonal_H)
               

            #print("Multiple_clones_H")
            #print(multiple_clones_H)

            # Get groups of clones sharing heavy and light chain
            """(paired_clone_groups, paired_clone_groups_contigs) = self.get_initial_clone_groups(self.loci, clones, multiple_clones_H, cell_clones, cell_contig_clones)
            print("Paired clone groups: ", paired_clone_groups)


            print("#####################   Paired clone groups   ############################\n")
            for clone, cell_names in six.iteritems(paired_clone_groups):
                print("\n---{}---".format(clone))
                
                if not cell_names is None:
                    for cell_name in cell_names:
                        print(cell_name)

            print("Paired clone groups contigs: ", paired_clone_groups_contigs)
            
            print("#####################   Paired clone groups contigs   ############################\n")
            for clone, locus_info in six.iteritems(paired_clone_groups_contigs):
                print("\n---{}---".format(clone))
                if locus in paired_clone_groups_contigs[clone].keys():
                    for locus, contigs in six.iteritems(paired_clone_groups_contigs[clone]):
                        print(locus)
                        if not contigs is None:
                            for contig in (paired_clone_groups_contigs[clone][locus]):
                                print(contig)
                        #for cell_name, info in six.iteritems(paired_clone_groups_contigs[clone][locus]):
                            #print(cell_name)
                            #for inf in info:
                                #print(inf)


            clonal_cells = []
            for clone, cell_list in six.iteritems(paired_clone_groups):
                for cell in cell_list:
                    clonal_cells.append(cell)
            #print("Clonal cells: ", clonal_cells)

            
            # Align clonal sequences with Muscle
            muscle = self.get_binary('muscle')
            for l in self.loci:
                # Create Muscle input file for each paired clone group for each locus
                self.make_muscle_clonal_alignment_input(outdir, l, self.receptor_name, cells, clonal_cells, paired_clone_groups, paired_clone_groups_contigs)
                # Run Muscle
                for clone, cell_list in six.iteritems(paired_clone_groups):
                    tracer_func.run_muscle(muscle, l, outdir, self.species, clone)
            
            # Create dictionary for sequence alignments for each clonal group
            (alignment_dict, first_cell_dict, differences_dict) = self.create_alignment_dict(paired_clone_groups, self.loci, outdir)
            #print("initial alignment dict", alignment_dict)
            #print("first cell dict", first_cell_dict)
            #print("differences_dict", differences_dict)

            # Modify alignment dictionary
            trimmed_alignment_dict = self.modify_alignment_dict(alignment_dict, first_cell_dict)

            (trimmed_alignment_dict, differences_dict) = self.get_differences_dict(trimmed_alignment_dict, first_cell_dict)
            #print(trimmed_alignment_dict)
            print()
            print("#####################   Differences_dict   ############################\n")
            for clone, locus_data in six.iteritems(differences_dict):
                print("\n---{}---".format(clone))
                 
                for locus, cell_data in six.iteritems(differences_dict[clone]):
                    print("\n**{}**".format(locus))
                    for cell_info, differences in six.iteritems(differences_dict[clone][locus]):
                        print(cell_info, differences)

            # Get distances between sequences in potential clonal groups
            matrix = self.load_distance_matrix(self.species)
            if matrix is not None:
                pass
            else:
                pass
                # Calculate Hamming distance?
            edit_distances = self.get_edit_distance(matrix, trimmed_alignment_dict, differences_dict)

            print()
            print("######################   Edit distances   ############################\n")
            for clone, locus_data in six.iteritems(edit_distances):
                print("\n---{}---".format(clone))

                for locus, cell_data in six.iteritems(edit_distances[clone]):
                    print("\n**{}**".format(locus))
                    for cell_info, differences in six.iteritems(edit_distances[clone][locus]):
                        print(cell_info, differences)

            # Normalise edit distances to sequence length
            n_edit_distances = self.get_normalised_edit_distance(edit_distances, trimmed_alignment_dict)

            print()
            print("######################   Normalised edit distances   ############################\n")
            for clone, locus_data in six.iteritems(n_edit_distances):
                print("\n---{}---".format(clone))

                for locus, cell_data in six.iteritems(n_edit_distances[clone]):
                    print("\n**{}**".format(locus))
                    for cell_info, differences in six.iteritems(n_edit_distances[clone][locus]):
                        print(cell_info, differences)

            
            # Get sum of normalised distances for all loci for each cell pair in each clone group
            s_edit_distances = self.get_sum_normalised_edit_distance(n_edit_distances)


            print()
            print("######################   Summarised normalised edit distances   ############################\n")
            for clone, cell_data in six.iteritems(s_edit_distances):
                print("\n---{}---".format(clone))

                for cell_info, distances in six.iteritems(s_edit_distances[clone]):
                    print(cell_info, distances)
    


            # Print output of initial clonal grouping
            outstring = "\n\n###Initial clonal groups determined by ChangeO###\n\n"
            clonal_cells = []
            for clone, cell_list in six.iteritems(paired_clone_groups):
                #print(clone, cell_list)
                string = ", ".join(cell_list) + "\n"
                outstring += string
                for cell in cell_list:
                    clonal_cells.append(cell)
            outfile.write(outstring)
            outfile.write("\n")"""
        
            # Make isotype usage table and plot isotype distributions
            #print("XXXXXXXXXXXXXXXXXXXXX PRINTING CELLS XXXXXXXXXXXXXXXXXXXXXXXXXXX")
            #print(cells)
            isotype_counter = self.count_isotype_usage(cells)
            (header, outstring) = self.make_isotype_table(cell_recovery_count, isotype_counter)
            outfile.write(header)
            outfile.write(outstring)
            outfile.write("\n")
            self.plot_isotype_distributions(isotype_counter, outdir)

        
        # reporting invariant cells
        invariant_cells = []
        
        if self.invariant_cells is not None:
             for ivc in self.invariant_cells:
                 ivc_loci = []
                 found_ivcs = {}
                 defining_locus = ivc.defining_locus
                 ivc_loci.append(defining_locus)
                 for cell in cells.values():
                     found_idents = {}
                     found_defining_locus, defining_id = ivc.check_for_match(cell, defining_locus)
                     if found_defining_locus:
                         found_idents[ivc.defining_locus] = defining_id
                         
                         for l in ivc.invariant_recombinants.keys():
                             if not l==defining_locus:
                                 ivc_loci.append(l)
                                 if l in cell.recombinants[ivc.receptor_type] and \
                                 cell.recombinants[ivc.receptor_type][l] is not None:
                                     found_other_locus, invar_id = ivc.check_for_match(cell, l)
                                     if found_other_locus:
                                         pass
                                     else:
                                         invar_id = "Invariant recombinant not found for {}_{}. {} found in total ({})".format(
                                             ivc.receptor_type, l, len(cell.recombinants[ivc.receptor_type][l]), 
                                             cell.getMainRecombinantIdentifiersForLocus(ivc.receptor_type, l))
                                         
                                 else:
                                     invar_id = "No sequences reconstructed for {}_{}".format(ivc.receptor_type, l)
                                 found_idents[l] = invar_id
                                     
                                 
                         found_ivcs[cell.name] = found_idents
                         invariant_cells.append(cell.name)
             
                 if len(found_ivcs) > 0:
                     outfile.write("\n#{} cells#\n".format(ivc.name))
                     
                     outfile.write("Expected: {}\n".format(ivc.expected_string))
                     outfile.write("Found {} possible cells.\n\n".format(len(found_ivcs)))
                     
                     sorted_names = sorted(list(found_ivcs.keys()))
                     for n in sorted_names:
                         outfile.write("### {} ###\n".format(n))
                         ivc_details = found_ivcs[n]
                         for l in ivc_loci:
                             outfile.write("{}_{}: {}\n".format(ivc.receptor_type, l, ivc_details[l]))
                     outfile.write("\n")
                 

        # plot lengths of reconstructed sequences
        lengths = defaultdict(list)
        for cell in cells.values():
            for l in self.loci:
                lengths[l].extend(cell.get_trinity_lengths(self.receptor_name, l))

        # plot length distributions
        quartiles = dict()
        for l in self.loci:
            q = self.get_quartiles(self.receptor_name, l)
            quartiles[l] = q
            
        for l in self.loci:
            q = quartiles[l]
            lns = lengths[l]
            if len(lns) > 1:
                plt.figure()
                plt.axvline(q[0], linestyle="--", color='k')
                plt.axvline(q[1], linestyle="--", color='k')
                sns.distplot(lns)
                sns.despine()
                plt.xlabel("{receptor}_{locus} reconstructed length (bp)".format(receptor=self.receptor_name,
                                                                                 locus=l))
                plt.ylabel("Density")
                plt.savefig("{}_{}.pdf".format(length_filename_root, l))
            if len(lns) > 0:
                with open("{}_{}.txt".format(length_filename_root,l), 'w') as f:
                        for l in sorted(lns):
                            f.write("{}\n".format(l))
                        

        for cell_name in empty_cells:
            del cells[cell_name]

        if not self.keep_invariant:
            for cell_name in invariant_cells:
                del cells[cell_name]
        

        # Write recombinant details
        with open("{}/recombinants.txt".format(outdir), 'w') as f:
            f.write("cell_name\tlocus\trecombinant_id\tproductive\treconstructed_length\n")
            sorted_cell_names = sorted(list(cells.keys()))
            for cell_name in sorted_cell_names:
                cell = cells[cell_name]
                for locus in self.loci:
                    recombinants = cell.recombinants[self.receptor_name][locus]
                    if recombinants is not None:
                        for r in recombinants:
                            f.write(
                                "{name}\t{locus}\t{ident}\t{productive}\t{length}\n".format(
                                    name=cell_name, locus=locus, ident=r.identifier,
                                    productive=r.productive, length=len(r.trinity_seq)))
                f.write("\n")
            f.write("\n\n")
            for cell_name in empty_cells:
                f.write("{cell_name}\tNo seqs found for {receptor}_{loci}\n".format(cell_name=cell_name,
                                                                                    receptor=self.receptor_name,
                                                                                    loci=self.loci))
                
        # make clonotype networks
        network_colours = io.read_colour_file(os.path.join(self.species_dir, "colours.csv"))
        if self.receptor_name == "BCR":
            component_groups = tracer_func.draw_network_from_cells(cells, outdir, self.graph_format, dot, neato,
                                                               self.draw_graphs, self.receptor_name, self.loci,
                                                               network_colours, cell_contig_clones, cells_with_clonal_H)
        
        # Print component groups to the summary#
        outfile.write(
            "\n#Clonotype groups#\n"
            "This is a text representation of the groups shown in clonotype_network_with_identifiers.pdf.\n"
            "It does not exclude cells that only share beta and not alpha.\n\n")
        for g in component_groups:
            outfile.write(", ".join(g))
            outfile.write("\n\n")
        
        """# plot clonotype sizes
        plt.figure()
        clonotype_sizes = tracer_func.get_component_groups_sizes(cells, self.receptor_name, self.loci)
        w = 0.85
        x_range = range(1, len(clonotype_sizes) + 1)
        plt.bar(x_range, height=clonotype_sizes, width=w, color='black', align='center')
        plt.gca().set_xticks(x_range)
        plt.xlabel("Clonotype size")
        plt.ylabel("Clonotype count")
        plt.savefig("{}/clonotype_sizes.pdf".format(outdir))
        
        # write clonotype sizes to text file
        with open("{}/clonotype_sizes.txt".format(outdir), 'w') as f:
            data = zip(x_range, clonotype_sizes)
            f.write("clonotype_size\tclonotype_count\n")
            for t in data:
                f.write("{}\t{}\n".format(t[0], t[1]))

        outfile.close()"""
    
    def get_quartiles(self, receptor, locus):
        
        fasta = os.path.join(self.species_dir, 'combinatorial_recombinomes','{receptor}_{locus}.fa'.format(
                                                                            receptor=receptor, locus=locus))
        
        # need to remove the start N padding and end C sequence from the lengths
        constant_fasta = os.path.join(self.species_dir, 'raw_seqs', '{receptor}_{locus}_C.fa'.format(
                                                                            receptor=receptor, locus=locus))
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
            
        
    def two_productive_chains_per_locus(self, loci, cells, receptor):
        """#Report cells with two productive chains from same locus"""
        outstring = ""
        for l in loci:
            double_rec = []
            cells_string = ""
            for cell_name, cell in six.iteritems(cells):
                prod_counts = dict()
                prod_counts[l] = cell.count_productive_recombinants(self.receptor_name, l)
                if prod_counts[l] == 2:
                    double_rec.append(cell_name)
            if len(double_rec) > 0 :
                cells_string = ', '.join(double_rec)
                locus_string = "Cells with two productive {} chains: \n".format(l) + cells_string + "\n\n"
                outstring += locus_string
        if len(outstring) != 0:
            outstring = "\n\n#Cells with two productive recombinants for a locus#\n\n" + outstring
        return (outstring)


    def report_kappa_lambda_cells(self, loci, cells, receptor):
        """Report cells with both productive kappa and lambda chain"""
        outstring = ""
        kappa_lambda = []
        cells_string = ""
        for cell_name, cell in six.iteritems(cells):
            store_cell = True
            for l in ["K", "L"]:
                prod_counts = dict()
                prod_counts[l] = cell.count_productive_recombinants(self.receptor_name, l)
                if prod_counts[l] == 0:
                    store_cell = False
            if store_cell == True:
                kappa_lambda.append(cell_name)
        if len(kappa_lambda) > 0 :
            cells_string = ', '.join(kappa_lambda)
            outstring = "\n\n#Cells with productive K and L chain#\n\n" + cells_string + "\n\n"
        return(outstring)


    def make_changeo_input(self, outdir, locus, receptor, cells):
        """Creates input file for each locus compatible with ChangeO"""
        
        changeo_string = "SEQUENCE_ID\tV_CALL\tD_CALL\tJ_CALL\tSEQUENCE_VDJ\tJUNCTION_LENGTH\tJUNCTION\n"
        changeo_input = "{}/changeo_input_{}.tab".format(outdir, locus)
        with open(changeo_input, 'w') as output:
            empty = True
            for cell_name, cell in six.iteritems(cells):
                productive = cell.count_productive_recombinants(receptor, locus)

                if productive > 0:
                    if empty == True:
                        output.write(changeo_string)
                        empty = False
                    output.write(cell.changeodict[locus])
        

    def read_changeo_results(self, loci, outdir):
        """Reads ChangeO result files and defines clone groups for each locus"""
        clones = dict()
        cell_clones = dict()
        cell_contig_clones = dict()
        cell_list = []
        contig_list = []

        for l in loci:
            clones[l] = dict()
            cell_clones[l] = dict()
            cell_contig_clones[l] = dict()
            
            changeo_result = "{}/changeo_input_{}_clone-pass.tab".format(outdir, l)
            if not os.path.exists(changeo_result):
                continue
            elif not os.path.getsize(changeo_result) > 0:
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


    """def get_initial_clone_groups(self, loci, clones, multiple_clones_H, cell_clones, cell_contig_clones):
        Get groups of B cell clones sharing clonally related heavy and light chains (from changeo output)
        paired_clone_groups = dict()
        paired_clone_groups_contigs = dict()
        # Create list of cells belonging to at least one heavy chain clone group with more than 1 members
        cell_list = []
        
        for clone, cells_info in six.iteritems(multiple_clones_H):    
            for cell_name, contig_name in six.iteritems(multiple_clones_H[clone]):
                if not cell_name in cell_list:
                    cell_list.append(cell_name)  

        for i in range(len(cell_list)):
            current_cell = cell_list[i]
            shared_H = False
            Hcontig2 = []
            Hcontig1 = []
            contig_list = cell_contig_clones["H"][current_cell].keys()
            for item in contig_list:
                contig_name = current_cell + "_" + item
                Hcontig2.append(contig_name)
            
            comparison_cells = cell_list[i + 1:]
            for comparison_cell in comparison_cells:

                comparison_contigs = cell_contig_clones["H"][comparison_cell].keys()
                current_contigs = cell_contig_clones["H"][current_cell].keys()
                for comparison_contig in comparison_contigs:
                    for current_contig in current_contigs:
                        if cell_contig_clones["H"][comparison_cell][comparison_contig] == cell_contig_clones["H"][current_cell][current_contig]:
                            shared_H = True
                
                    
                for l in self.loci:
                    if not comparison_cell in cell_contig_clones[l].keys():
                        cell_contig_clones[l][comparison_cell] = None

                    if not current_cell in cell_contig_clones[l].keys():
                        cell_contig_clones[l][current_cell] = None

                if shared_H == True:
                    clone = False
                    Kcontig1 = []
                    Kcontig2 = []
                    Lcontig1 = []
                    Lcontig2 = []
                    shared_K = False
                    shared_L = False

                    # Filter out cells with no light chain - make optional in future!
                    if ((cell_contig_clones["K"][comparison_cell] is None) or (cell_contig_clones["K"][current_cell] is None)) and \
                            ((cell_contig_clones["L"][comparison_cell] is None) or (cell_contig_clones["L"][current_cell] is None)):
                        clone = False

                    # Check for shared kappa chains
                    if (cell_contig_clones["K"][comparison_cell] is not None) and (cell_contig_clones["K"][current_cell] is not None):
                        comparison_contigs = cell_contig_clones["K"][comparison_cell].keys()
                        current_contigs = cell_contig_clones["K"][current_cell].keys()
                        for comparison_contig in comparison_contigs:
                            for current_contig in current_contigs:
                                if cell_contig_clones["K"][comparison_cell][comparison_contig] == cell_contig_clones["K"][current_cell][current_contig]:
                                    clone = True
                                    contig = current_cell + "_" + current_contig
                                    Kcontig1.append(contig)
                                    contig = comparison_cell + "_" + comparison_contig
                                    Kcontig2.append(contig)
                                    shared_K = True

                    # Check for shared lambda chains
                    if (cell_contig_clones["L"][comparison_cell] is not None) and (cell_contig_clones["L"][current_cell] is not None):
                        comparison_contigs = cell_contig_clones["L"][comparison_cell].keys()
                        current_contigs = cell_contig_clones["L"][current_cell].keys()
                        for comparison_contig in comparison_contigs:
                            for current_contig in current_contigs:
                                if cell_contig_clones["L"][comparison_cell][comparison_contig] == cell_contig_clones["L"][current_cell][current_contig]:
                                    clone = True
                                    contig = current_cell + "_" + current_contig
                                    Lcontig1.append(contig)
                                    contig = comparison_cell + "_" + comparison_contig
                                    Lcontig2.append(contig)
                                    shared_L = True

                    # Add to dictionary if cells are clonal

                    if clone == True and shared_H == True:
                        found = False
                        clones_so_far = len(paired_clone_groups.keys())
                        if clones_so_far == 0:
                            paired_clone_groups[1] = [current_cell, comparison_cell]
                            paired_clone_groups_contigs[1] = dict()
                            paired_clone_groups_contigs[1]["K"] = None
                            paired_clone_groups_contigs[1]["L"] = None
                            if shared_K:
                                paired_clone_groups_contigs[1]["K"] = [Kcontig1, Kcontig2]
                            if shared_L:
                                paired_clone_groups_contigs[1]["L"] = [Lcontig1, Lcontig2]
                            paired_clone_groups_contigs[1]["H"] = [Hcontig1, Hcontig2]

                        else:
                            for i in range(1, clones_so_far+1):
                                if (current_cell or comparison_cell) in paired_clone_groups[i]:
                                    if paired_clone_groups_contigs[i]["K"] is not None:
                                        if (Kcontig1 or Kcontig2) in paired_clone_groups_contigs[i]["K"]:
                                            found = True
                                            if not Kcontig1 in paired_clone_groups_contigs[i]["K"]:
                                                paired_clone_groups_contigs[i]["K"].append(Kcontig1)
                                            elif not  Kcontig2 in paired_clone_groups_contigs[i]["K"]:
                                                paired_clone_groups_contigs[i]["K"].append(Kcontig2)
                                    if paired_clone_groups_contigs[i]["L"] is not None:
                                        if (Lcontig1 or Lcontig2) in paired_clone_groups_contigs[i]["L"]:
                                            found = True
                                            if not Lcontig1 in paired_clone_groups_contigs[i]["L"]:
                                                paired_clone_groups_contigs[i]["L"].append(Lcontig1)
                                            elif not  Lcontig2 in paired_clone_groups_contigs[i]["L"]:
                                                paired_clone_groups_contigs[i]["L"].append(Lcontig2)
                                            found = True
                                    if paired_clone_groups_contigs[i]["H"] is not None:
                                        if (Hcontig1 or Hcontig2) in paired_clone_groups_contigs[i]["H"]:
                                            found = True
                                            if not Hcontig1 in paired_clone_groups_contigs[i]["H"]:
                                                paired_clone_groups_contigs[i]["H"].append(Hcontig1)
                                            elif not  Hcontig2 in paired_clone_groups_contigs[i]["H"]:
                                                paired_clone_groups_contigs[i]["H"].append(Hcontig2)
                                    if not current_cell in paired_clone_groups[i]:
                                        paired_clone_groups[i].append(current_cell)

                                    elif not comparison_cell in paired_clone_groups[i]:
                                        paired_clone_groups[i].append(comparison_cell)


                            if not found == True:
                                paired_clone_groups[clones_so_far + 1] = [current_cell, comparison_cell]
                                paired_clone_groups_contigs[clones_so_far + 1] = dict()
                                paired_clone_groups_contigs[clones_so_far + 1]["K"] = None
                                paired_clone_groups_contigs[clones_so_far + 1]["L"] = None
                                if shared_K:
                                    paired_clone_groups_contigs[clones_so_far + 1]["K"] = [Kcontig1, Kcontig2]
                                if shared_L:
                                    paired_clone_groups_contigs[clones_so_far + 1]["L"] = [Lcontig1, Lcontig2]
                                paired_clone_groups_contigs[clones_so_far + 1]["H"] = [Hcontig1, Hcontig2]
        return (paired_clone_groups, paired_clone_groups_contigs)"""


    def get_initial_clone_groups(self, loci, clones, multiple_clones_H, cell_clones, cell_contig_clones):
        """Get groups of B cell clones sharing clonally related heavy and light chains (from changeo output)"""
        paired_clone_groups = dict()
        paired_clone_groups_contigs = dict()
        for clone, cells in six.iteritems(multiple_clones_H):
            cell_list = []
            for cell, data in six.iteritems(multiple_clones_H[clone]):
                cell_list.append(cell)
            
            for i in range(len(cell_list)):
                current_cell = cell_list[i]

                Hcontig2 = []
                Hcontig1 = []
                contig_list = cells[current_cell]
                for item in contig_list:
                    contig_name = current_cell + "_" + item
                    Hcontig2.append(contig_name)

                comparison_cells = cell_list[i + 1:]
                for comparison_cell in comparison_cells:
                    #Hcontig1 = []
                   
                    contig_list = cells[comparison_cell]
                    for item in contig_list:
                        contig_name = comparison_cell + "_" + item
                        Hcontig1.append(contig_name)
                    for l in self.loci:
                        if not comparison_cell in cell_contig_clones[l].keys():
                            cell_contig_clones[l][comparison_cell] = None
                            
                        if not current_cell in cell_contig_clones[l].keys():
                            cell_contig_clones[l][current_cell] = None

                    clone = False
                    Kcontig1 = []
                    Kcontig2 = []
                    Lcontig1 = []
                    Lcontig2 = []
                    shared_K = False
                    shared_L = False
                    
                    # Filter out cells with no light chain - make optional in future!
                    if ((cell_contig_clones["K"][comparison_cell] is None) or (cell_contig_clones["K"][current_cell] is None)) and \
                            ((cell_contig_clones["L"][comparison_cell] is None) or (cell_contig_clones["L"][current_cell] is None)):
                        clone = False

                    # Check for shared kappa chains
                    if (cell_contig_clones["K"][comparison_cell] is not None) and (cell_contig_clones["K"][current_cell] is not None):
                        comparison_contigs = cell_contig_clones["K"][comparison_cell].keys()
                        current_contigs = cell_contig_clones["K"][current_cell].keys()
                        for comparison_contig in comparison_contigs:
                            for current_contig in current_contigs:
                                if cell_contig_clones["K"][comparison_cell][comparison_contig] == cell_contig_clones["K"][current_cell][current_contig]:
                                    clone = True
                                    contig = current_cell + "_" + current_contig
                                    Kcontig1.append(contig)
                                    contig = comparison_cell + "_" + comparison_contig
                                    Kcontig2.append(contig)
                                    shared_K = True
                     
                    # Check for shared lambda chains
                    if (cell_contig_clones["L"][comparison_cell] is not None) and (cell_contig_clones["L"][current_cell] is not None):
                        comparison_contigs = cell_contig_clones["L"][comparison_cell].keys()
                        current_contigs = cell_contig_clones["L"][current_cell].keys()
                        for comparison_contig in comparison_contigs:
                            for current_contig in current_contigs:
                                if cell_contig_clones["L"][comparison_cell][comparison_contig] == cell_contig_clones["L"][current_cell][current_contig]:
                                    clone = True
                                    contig = current_cell + "_" + current_contig
                                    Lcontig1.append(contig)
                                    contig = comparison_cell + "_" + comparison_contig
                                    Lcontig2.append(contig)
                                    shared_L = True
         
                    # Add to dictionary if cells are clonal 

                    if clone == True:
                        found = False
                        clones_so_far = len(paired_clone_groups.keys())
                        if clones_so_far == 0:
                            paired_clone_groups[1] = [current_cell, comparison_cell]
                            paired_clone_groups_contigs[1] = dict()
                            paired_clone_groups_contigs[1]["K"] = None
                            paired_clone_groups_contigs[1]["L"] = None
                            if shared_K:
                                paired_clone_groups_contigs[1]["K"] = [Kcontig1, Kcontig2]
                            if shared_L:
                                paired_clone_groups_contigs[1]["L"] = [Lcontig1, Lcontig2]
                            paired_clone_groups_contigs[1]["H"] = [Hcontig1, Hcontig2]
                        else:
                            for i in range(1, clones_so_far+1):
                                if (current_cell or comparison_cell) in paired_clone_groups[i]:
                                    if paired_clone_groups_contigs[i]["K"] is not None:
                                        if (Kcontig1 or Kcontig2) in paired_clone_groups_contigs[i]["K"]:
                                            found = True
                                            if not Kcontig1 in paired_clone_groups_contigs[i]["K"]:
                                                paired_clone_groups_contigs[i]["K"].append(Kcontig1)
                                            elif not  Kcontig2 in paired_clone_groups_contigs[i]["K"]:
                                                paired_clone_groups_contigs[i]["K"].append(Kcontig2)  
                                    if paired_clone_groups_contigs[i]["L"] is not None:                                 
                                        if (Lcontig1 or Lcontig2) in paired_clone_groups_contigs[i]["L"]:
                                            found = True
                                            if not Lcontig1 in paired_clone_groups_contigs[i]["L"]:
                                                paired_clone_groups_contigs[i]["L"].append(Lcontig1)
                                            elif not  Lcontig2 in paired_clone_groups_contigs[i]["L"]:
                                                paired_clone_groups_contigs[i]["L"].append(Lcontig2) 
                                            found = True
                                    if paired_clone_groups_contigs[i]["H"] is not None:
                                        if (Hcontig1 or Hcontig2) in paired_clone_groups_contigs[i]["H"]:
                                            found = True
                                            if not Hcontig1 in paired_clone_groups_contigs[i]["H"]:
                                                paired_clone_groups_contigs[i]["H"].append(Hcontig1)
                                            elif not  Hcontig2 in paired_clone_groups_contigs[i]["H"]:
                                                paired_clone_groups_contigs[i]["H"].append(Hcontig2)
                                    if not current_cell in paired_clone_groups[i]:
                                        paired_clone_groups[i].append(current_cell)
                                        
                                    elif not comparison_cell in paired_clone_groups[i]:
                                        paired_clone_groups[i].append(comparison_cell)
                                        
                                        
                            if not found == True:
                                paired_clone_groups[clones_so_far + 1] = [current_cell, comparison_cell]
                                paired_clone_groups_contigs[clones_so_far + 1] = dict()
                                paired_clone_groups_contigs[clones_so_far + 1]["K"] = None
                                paired_clone_groups_contigs[clones_so_far + 1]["L"] = None
                                if shared_K:
                                    paired_clone_groups_contigs[clones_so_far + 1]["K"] = [Kcontig1, Kcontig2]
                                if shared_L:
                                    paired_clone_groups_contigs[clones_so_far + 1]["L"] = [Lcontig1, Lcontig2]
                                paired_clone_groups_contigs[clones_so_far + 1]["H"] = [Hcontig1, Hcontig2]
        return (paired_clone_groups, paired_clone_groups_contigs)



    def make_muscle_clonal_alignment_input(self, outdir, locus, receptor, cells, clonal_cells, paired_clone_groups, paired_clone_groups_contigs):
        """Creates input file for each locus only containing sequences from cells belonging to paired clonal groups"""

        changeo_output = "{}/changeo_input_{}_clone-pass.tab".format(outdir, locus)
        if os.path.exists(changeo_output):
       
            for clone, locus_data in six.iteritems(paired_clone_groups_contigs):
                print()
                print(clone)
                muscle_input = "{}/muscle_input_{}_{}.fa".format(outdir, locus, clone)
                muscle_input2 = "{}/muscle_input_{}_{}b.fa".format(outdir, locus, clone)
                with open(muscle_input, 'w') as output:
                    if paired_clone_groups_contigs[clone][locus] is not None:
                        print()
                        print(locus)
                        print()
                        
                        with open(changeo_output, 'r') as input:
                            contig_list = []
                            contig_list2 = []
                            print(paired_clone_groups_contigs[clone][locus])
                            for item in paired_clone_groups_contigs[clone][locus]:
                                print(item)
                                for item2 in item:
                                    print("Item 2", item2)
                                    # Check if two clonal chains for locus
                                    if len(item2) == 2:
                                        print("Len item2 is 2")
                                        print("First item", item2[0], "Second item", item2[1])
                                        contig_list.append(item2[0])
                                        contig_list2.append(item2[1])
                                        print("appended items to contig lists")
                                        print("Contig list", contig_list, "Contig list 2", contig_list2)
                                    else:
                                        contig_list.append(item2)
                                        #print("Len item2 is not 2")
                                        print("Contig list", contig_list)
                            # Check if two clonal chains for locus (contig_list2 not empty)
                            if len(contig_list2) > 0:
                                fh = open(muscle_input2, "w")
                                print("Multiple clonal chains for locus detected")

                                
                            for line in input:
                                if not line.startswith("SEQUENCE_ID"):
                                    cell_name = line.split("_TRINITY")[0]
                                    contig_name = line.split("\t")[0].split("{}_".format(cell_name))[1]
                                    cell_contig_name = line.split("\t")[0]
                                    #print(cell_contig_name)
                                    

                                    if cell_contig_name in (contig_list or contig_list2):
                                        sequence = line.split()[4] + "\n"
                                        fasta_header = ">" + line.split()[0] + "\n"
                                        if cell_contig_name in contig_list:
                                            output.write(fasta_header)
                                            output.write(sequence)
                                        else:
                                            print("Writing to extra file")
                                            fh.write(fasta_header)
                                            fh.write(sequence)
                            
                            if os.path.exists(muscle_input2):
                                fh.close()

               


    def create_alignment_dict(self, paired_clone_groups, loci, outdir):
        alignment_dict = dict()
        differences_dict = dict()
        first_cell_dict = dict()
        for clone, cell_list in six.iteritems(paired_clone_groups):
            alignment_dict[clone] = dict()
            differences_dict[clone] = dict()
            first_cell_dict[clone] = dict()
            for l in loci:
                muscle_result_file = "{}/muscle_out_{}_{}.aln".format(outdir, l, clone)
                if not os.path.exists(muscle_result_file):
                    continue
                alignment_dict[clone][l] = dict()
                differences_dict[clone][l] = dict()
                first_cell = False
                with open(muscle_result_file, 'r') as infile:
                    count = 0
                    for line in infile:
                        cell_name = None
                        line = line.lstrip()
                        if line.startswith("MUSCLE"):
                            continue
                        if (line.startswith("\n") or len(line) == 0) and count == 0:
                            continue
                        if (line.startswith("\n") or len(line) == 0) and count > 0:
                            count += 1
                            if count > 2:
                                cell_name = "summary"
                                line = ""
                            
                        elif line.startswith("*"):
                            cell_name = "summary"
                            line = line.lstrip()
                            count = 0
                        else:
                            cell_name = line.split("_TRINITY")[0]
                            start = line.find("  ")
                            line = line[start:].lstrip()
                            count = 1
                            if first_cell == False:
                                first_cell = cell_name
                                first_cell_dict[clone][l] = first_cell
                                
                        if not cell_name in alignment_dict[clone][l].keys():
                            if not cell_name is None:
                                alignment_dict[clone][l][cell_name] = [line]
                        else:
                            if not cell_name is None:
                                alignment_dict[clone][l][cell_name].append(line)

        return (alignment_dict, first_cell_dict, differences_dict)


    def modify_alignment_dict(self, alignment_dict, first_cell_dict):
        trimmed_alignment_dict = dict()
        for clone, data in six.iteritems(alignment_dict):
            trimmed_alignment_dict[clone] = dict()
            for l, l_data in six.iteritems(alignment_dict[clone]):
                trimmed_alignment_dict[clone][l] = dict()
                alignment_string = alignment_dict[clone][l]
                
                num_lines = len(alignment_string["summary"])
                skip_lines = 0
                for line in alignment_string["summary"]:
                    if len(line) > 0:
                        break
                    else:
                        skip_lines += 1
                
                first_cell = first_cell_dict[clone][l]
                total_lines = len(alignment_string[first_cell])

                for cell_name, alignment in six.iteritems(alignment_dict[clone][l]):
                    new_alignment = ""
                     
                    if cell_name == "summary":
                        for i in range(skip_lines, total_lines):
                            difference = len(alignment_dict[clone][l][first_cell][i]) - len(alignment_dict[clone][l]["summary"][i])
                            if i == skip_lines:
                                add = alignment[i][:len(alignment[i]) -1]
                            else:
                                if difference > 0:
                                    add = difference*" " +  alignment[i][:len(alignment[i])-1]
                                
                                else:
                                    add = alignment[i][:len(alignment[i]) -1]
                            if len(new_alignment) == 0:
                                new_alignment = add
                            else:
                                new_alignment += add
                    else:
                        for i in range(skip_lines, total_lines):
                            difference = len(alignment_dict[clone][l][first_cell][i]) - len(alignment_dict[clone][l]["summary"][i])
                            if i == skip_lines:
                                add = alignment[i][difference:len(alignment[i]) -1]
                            else:
                                add = alignment[i][:len(alignment[i]) -1]  

                            if len(new_alignment) == 0:
                                new_alignment = add
                            else:
                                new_alignment += add            
                    
                    trimmed_alignment_dict[clone][l][cell_name] = new_alignment
 
        
        return(trimmed_alignment_dict)

    def get_differences_dict(self, trimmed_alignment_dict, first_cell_dict):
        differences_dict = dict()
        for clone, data in six.iteritems(trimmed_alignment_dict):
            differences_dict[clone] = dict()
            for l, l_data in six.iteritems(trimmed_alignment_dict[clone]):
                differences_dict[clone][l] = dict()
                # Identify polymorphic sites in sequence alignments
                summary_string = trimmed_alignment_dict[clone][l]["summary"]
                length = len(summary_string)
                polymorphic = []
                for i in range(0, length):
                    if summary_string[i] is not "*":
                        polymorphic.append(i)
    
                for i in range(0, len(polymorphic)):
                    for cell_name, alignment in six.iteritems(trimmed_alignment_dict[clone][l]):
                        if i == 0 and not cell_name == "summary":
                            differences_dict[clone][l][cell_name] = [alignment[polymorphic[0]]]
                        elif not cell_name == "summary":
                            differences_dict[clone][l][cell_name].append(alignment[polymorphic[i]])
                del(trimmed_alignment_dict[clone][l]["summary"])
        return(trimmed_alignment_dict, differences_dict)


    def load_distance_matrix(self, species):

        matrix = dict()
        if self.species == "Mmus":
            # Use M1N substitution distance model (from ChangeO)
                matrix["A"] = {"A":0, "C":2.86, "G":1, "T":2.14, "N":0, "-":0}
                matrix["C"] = {"A":2.86, "C":0, "G":2.14, "T":1, "N":0, "-":0}
                matrix["G"] = {"A":1, "C":2.14, "G":0, "T":2.86, "N":0, "-":0}
                matrix["T"] = {"A":2.14, "C":1, "G":2.86, "T":0, "N":0, "-":0}
                matrix["N"] = {"A":0, "C":0, "G":0, "T":0, "N":0, "-":0}
                matrix["-"] = {"A":0, "C":0, "G":0, "T":0, "N":0, "-":0}
            # Use hh_s1f substitution matrix from ChangeO
        elif species == "Hsap":
                matrix["A"] = {"A":0, "C":1.21, "G":0.64, "T":1.16, "N":0, "-":0}
                matrix["C"] = {"A":1.21, "C":0, "G":1.16, "T":0.64, "N":0, "-":0}
                matrix["G"] = {"A":0.64, "C":1.16, "G":0, "T":1.21, "N":0, "-":0}
                matrix["T"] = {"A":1.16, "C":0.64, "G":1.21, "T":0, "N":0, "-":0}
                matrix["N"] = {"A":0, "C":0, "G":0, "T":0, "N":0, "-":0}
                matrix["-"] = {"A":0, "C":0, "G":0, "T":0, "N":0, "-":0}

        else:
            matrix = None
            # Use Hamming distance instead of SHM distance model

        return(matrix)



    def get_edit_distance(self, matrix, alignment_dict, differences_dict):
        """Returns edit distance between sequences in each clone group for each locus"""
        edit_distances = dict()
        for clone, clone_data in six.iteritems(alignment_dict):
            edit_distances[clone] = dict()
            for l in list(alignment_dict[clone].keys()):
                equal = False
                cell_list = list(alignment_dict[clone][l].keys())
                edit_distances[clone][l] = dict()
                if differences_dict[clone][l] == {}:
                    equal = True

                for i in range(len(cell_list)- 1):
                    current_cell = cell_list[i]
                    
                    comparison_cells = cell_list[i + 1:]
                    for comparison_cell in comparison_cells:
                        pair = (current_cell, comparison_cell)
                        seq_differences = differences_dict[clone][l]
                        distance = 0
                        if equal == True:
                            distance = 0
                        
                        elif seq_differences[current_cell] == seq_differences[comparison_cell]:
                            distance = 0
                        elif not (current_cell) in seq_differences.keys():
                            distance = 0
                        else:
                            for n in range(len(seq_differences[current_cell])):
                                nt1 = seq_differences[current_cell][n]
                                nt2 = seq_differences[comparison_cell][n]
                                if nt1 != nt2:
                                    distance = float(distance) + float(matrix[nt1][nt2])
                        edit_distances[clone][l][pair] = distance
        return(edit_distances)

    def get_normalised_edit_distance(self, edit_distances, alignment_dict):
        n_edit_distance = dict()
        for clone, data in six.iteritems(edit_distances):
            n_edit_distance[clone] = dict()
            for l, l_data in six.iteritems(edit_distances[clone]):
                n_edit_distance[clone][l] = dict()
                for (current_cell, comparison_cell), distance in six.iteritems(edit_distances[clone][l]):
                    pair = (current_cell, comparison_cell)
                    
                    
                    seq_length = len(alignment_dict[clone][l][current_cell])
                    n_distance = distance / seq_length * 100
                    n_edit_distance[clone][l][pair] = n_distance
        return(n_edit_distance)

    def get_sum_normalised_edit_distance(self, n_edit_distances):
        s_edit_distance = dict()
        for clone, data in six.iteritems(n_edit_distances):
            s_edit_distance[clone] = dict()
            for l, l_data in six.iteritems(n_edit_distances[clone]):
                for (current_cell, comparison_cell), distance in six.iteritems(n_edit_distances[clone][l]):
                   
                    
                    pair1 = (current_cell, comparison_cell)
                    pair2 = (comparison_cell, current_cell)
                    if not pair1 in s_edit_distance[clone].keys():
                        if not pair2 in s_edit_distance[clone].keys():
                            s_edit_distance[clone][pair1] = distance
                        else:
                            s_edit_distance[clone][pair2] += distance
                    else:
                        s_edit_distance[clone][pair1] += distance
                        
        return(s_edit_distance)


    def count_full_length_sequences(self, loci, cells, receptor):
        """Count all full length sequences and productive full length sequences"""
        full_length_counter = dict()
        full_length_prod_counter = dict()
        productive_counter = dict()
        total_counter = dict()
        
        for l in self.loci:
            (full_length_counter[l], full_length_prod_counter[l], productive_counter[l], total_counter[l]) = (0, 0, 0, 0)

        for cell_name, cell in six.iteritems(cells):
            for l in self.loci:
                full_length_count = cell.count_full_length_recombinants(self.receptor_name, l)
                full_length_prod_count = cell.count_productive_full_length_recombinants(self.receptor_name, l)
                productive = cell.count_productive_recombinants(self.receptor_name, l)
                total = cell.count_total_recombinants(self.receptor_name, l)
                if full_length_count > 0:
                    full_length_counter[l] += full_length_count
                if full_length_prod_count > 0:
                    full_length_prod_counter[l] +=full_length_prod_count
                if productive > 0:
                    productive_counter[l] += productive
                if total > 0:
                    total_counter[l] += total
        return (full_length_counter, full_length_prod_counter, productive_counter, total_counter)


    def plot_full_length_sequences_proportions(self, all_dict, prod_dict, loci, outdir):
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

        rects1 = plt.bar(index, values_all, bar_width, color='#000000', label='All')
        rects2 = plt.bar(index + bar_width, values_prod, bar_width, color='#cccccc', label='Prod')

        plt.xlabel('Locus')
        plt.ylabel('Percentage')
        plt.title('Percentage full-length of all or productive sequences')
        plt.xticks(index + bar_width, x_ticks)
        plt.legend()

        #plt.tight_layout()
        plt.savefig("{}/full_length_seqs.pdf".format(outdir))


    def make_full_length_statistics_table(self, loci, total_counter, productive_counter, full_length_counter, full_length_prod_counter):
        """ Make full length statistics table"""

        header = "##Proportion of full-length sequences of all recovered sequences##\n\n\t"
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

            all_string = "{}/{} ({}%)\t".format(full_length_seqs, all_seqs, all_percent)
            prod_string = "{}/{} ({}%)\t".format(prod_full_length_seqs, prod_seqs, prod_percent)
            all_outstring += all_string
            prod_outstring += prod_string

        outstring = "{header}\n{all_outstring}\n{prod_outstring}\n".format(header=header, all_outstring=all_outstring, prod_outstring=prod_outstring)
        full_length_statistics_table = outstring
        return (full_length_statistics_table, all_dict, prod_dict)


    def make_full_length_plots(self):
        # Plot proportions of sequences that are full-length

        """D_prod = prod_cdr3_counter
        highest = None
        lowest = None

        H_values = int(D_prod["H"])
        K_values = int(D_prod["K"])
        L_values = int(D_prod["L"])
        all_values = H_values + K_values + L_values
        highest = max(all_values)
        lowest = min(all_values)

        n_groups = len(range(lowest, highest+1))

        fig, ax = plt.subplots()
        x_ticks = tuple(range(lowest, highest+1, 1))
        index = np.arange(n_groups)
        bar_width = 0.25
        opacity = 0.2

        rects1 = plt.bar(index, tuple(H_values), bar_width, color='#000000', label='H')
        rects2 = plt.bar(index + bar_width, tuple(K_values), bar_width, color='#cccccc', label='K')
        rects3 = plt.bar(index + bar_width + bar_width, tuple(L_values), bar_width, color='#cccccc', label='L')

        plt.xlabel('Locus')
        plt.ylabel('Frequency')
        plt.title('CDR3 length distribution (aa)')
        plt.xticks(index + bar_width, x_ticks)
        plt.legend()

        #plt.tight_layout()
        plt.savefig("{}/full_length_seqs.pdf".format(outdir))


        #D = prod_cdr3_counter
        #for l in self.loci:
            #D = dictionary[l]
            #lengths = []
            #counts = []


        for l in self.loci:
            lns = lengths[l]
            if len(lns) > 1:
                plt.figure()
                sns.distplot(lns)
                sns.despine()
                plt.xlabel("{receptor}_{locus} CDR3 length (aa)".format(receptor=self.receptor_name,
                                                                                 locus=l))
                plt.ylabel("Density")
                plt.savefig("{}_{}.pdf".format(cdr3_filename_root, l))
            if len(lns) > 0:
                with open("{}_{}.txt".format(cdr3_filename_root,l), 'w') as f:
                        for l in sorted(lns):
                            f.write("{}\n".format(l))"""

    def count_isotype_usage(self, cells):
        """Counts isotype usage of cells"""

        prod_counters = defaultdict(Counter)
        cell_isotypes = []
        isotype_counter = dict()
        for cell in cells.values():

            #for cell_name, cell in six.iteritems(cells):
            """l = "H"
            productive = cell.count_productive_recombinants(self.receptor_name, l)
            prod_counters[l].update({productive: 1})"""
            isotype = cell.isotype
            if isotype == None:
                isotype = "None"
            cell_isotypes.append(isotype)
            #if productive > 0:
                #cell_isotypes.append(isotype)
        for isotype in cell_isotypes:
            if not isotype in isotype_counter:
                isotype_counter[isotype] = 1
            else:
                isotype_counter[isotype] += 1
        return (isotype_counter)

    def make_isotype_table(self, cell_recovery_count, isotype_counter):
        prod_H = cell_recovery_count["H"]
        header = "##Isotype of cells with productive heavy chain##\n\nIsotype\tcells\t% of cells\n"
        outstring = ""
        for isotype, number in six.iteritems(isotype_counter):
            number = str(number)
            if isotype == "None":
                isotype = "Unknown"
            percent = float(number)/int(prod_H)*100
            percent = format(percent, '.2f')
            string = "{isotype}\t{number}\t{percent}\n".format(isotype=isotype, number=number, percent=percent)
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
            plt.bar(range(len(D)), D.values(), width=w, color='black', align='center')
            plt.xticks(range(len(D)), list(D.keys()))
            plt.xlabel("Isotype")
            plt.ylabel("Cell count")
            plt.savefig("{}/isotype_distribution.pdf".format(outdir))


    def count_cdr3_length_distributions(self):
        pass
        """#count cdr3 length distributions

        prod_counters = defaultdict(Counter)
        #all_cdr3_counter = dict()
        prod_cdr3_counter = dict()
        for l in self.loci:

            #all_cdr3_counter[l] = dict()
            prod_cdr3_counter[l] = dict()

        for cell_name, cell in six.iteritems(cells):
            for l in self.loci:
                prod_lengths = cell.get_prod_cdr3_lengths(self.receptor_name, l)

                for length in prod_lengths:
                    if not length in prod_cdr3_counter[l]:
                        prod_cdr3_counter[l][length] = 1
                    else:
                        prod_cdr3_counter[l][length] += 1

        #print (all_cdr3_counter)
        print (prod_cdr3_counter)"""

    def plot_cdr3_length_distributions(self):
        pass

        """# plot cdr3 distributions

        dictionary = prod_cdr3_counter
        for l in self.loci:
            D = dictionary[l]
            lengths = []
            counts = []


            for length, count in six.iteritems(D):
                lengths.append(length)
                counts.append(count)
                shortest = min(lengths)
                longest = max(lengths)

            if len(lengths) > 1:
                plt.figure()
                w = 0.85
                plt.bar(range(shortest, longest), D.values(), width=w, color='black', align='center')
                plt.xticks(range(len(D)), 1)
                plt.xlabel("CDR3 length (aa)")
                plt.ylabel("Frequency")
                plt.savefig("{}/{}cdr3_distribution.pdf".format(outdir, locus))
        # plot cdr3 distributions

        lengths = defaultdict(list)
        for cell in cells.values():
            for l in self.loci:
                print(cell)
                print(cell.get_prod_cdr3_lengths(self.receptor_name, l))
                lengths[l].extend(cell.get_prod_cdr3_lengths(self.receptor_name, l))"""


class Tester(TracerTask):

    def __init__(self, **kwargs):
        if not kwargs:
            parser = argparse.ArgumentParser(description="Test TraCeR installation with small dataset",
                                             parents=[self.base_parser])
            parser.add_argument('--graph_format', '-f', metavar="<GRAPH_FORMAT>", help='graphviz output format [pdf]',
                                default='pdf')
            parser.add_argument('--no_networks', help='skip attempts to draw network graphs', action="store_true")
            parser.add_argument('--resume_with_existing_files', '-r',
                                help='look for existing intermediate files and use those instead of starting from scratch',
                                action="store_true")
            args = parser.parse_args(sys.argv[2:])

            self.ncores = args.ncores
            self.config_file = args.config_file
            self.graph_format = args.graph_format
            self.no_networks = args.no_networks
            self.resume = args.resume_with_existing_files
        else:
            self.ncores = kwargs.get('ncores')
            self.config_file = kwargs.get('config_file')
            self.graph_format = kwargs.get('graph_format', 'pdf')
            self.no_networks = kwargs.get('no_networks')
            self.resume = kwargs.get('resume_with_existing_files')

    def run(self):

        # test_dir = self.resolve_relative_path("test_data")
        test_dir = os.path.join(base_dir, 'test_data')
        test_names = ['cell1']
        out_dir = "{}/results".format(test_dir)
        for name in test_names:
            f1 = "{}/{}_1.fastq".format(test_dir, name)
            f2 = "{}/{}_2.fastq".format(test_dir, name)

            Assembler(ncores=str(self.ncores), config_file=self.config_file, resume_with_existing_files=self.resume,
                      species='Mmus', seq_method='imgt', fastq1=f1, fastq2=f2, cell_name=name, output_dir=out_dir,
                      single_end=False, fragment_length=False, fragment_sd=False, receptor_name='TCR',
                      loci=['A', 'B'], max_junc_len=50).run()

        Summariser(config_file=self.config_file, use_unfiltered=False, keep_invariant=False,
                   graph_format=self.graph_format, no_networks=self.no_networks, root_dir=out_dir, receptor_name='TCR',
                   loci=['A', 'B'], species='Mmus').run()


class Builder(TracerTask):

    """ Build Combinatorial Recombinomes for a given species """

    def __init__(self, **kwargs):
        self.leader_padding = 20
        if not kwargs:
            parser = argparse.ArgumentParser(description="Build resources from sequences", parents=[self.base_parser])
            parser.add_argument('--force_overwrite', '-f', help='force overwrite of existing resources',
                                action='store_true')
            parser.add_argument('species', metavar="<SPECIES>", help='species (eg Mmus)')
            parser.add_argument('receptor_name', metavar="<RECEPTOR_NAME>", help='name of receptor (eg TCR)')
            parser.add_argument('locus_name', metavar="<LOCUS_NAME>", help='name of locus (eg A)')
            parser.add_argument('N_padding', metavar="<N_PADDING>", 
                                 help='number of ambiguous N nucleotides between V and J', type=int)
            parser.add_argument('colour', metavar="<COLOUR>", default = 'random', 
                                help='colour for productive recombinants. Specify as HTML (eg E41A1C)\
                                or use "random"', type = self.check_colour)
            parser.add_argument('V_seqs', metavar="<V_SEQS>", help='fasta file containing V gene sequences')
            parser.add_argument('J_seqs', metavar="<J_SEQS>", help='fasta file containing J gene sequences')
            parser.add_argument('C_seqs', metavar="<C_SEQS>", 
                                help='fasta file containing C gene sequence(s) for creation of recombinomes')

            parser.add_argument('D_seqs', metavar="<D_SEQS>", nargs='?', default=False,
                                help='fasta file containing D gene sequences (optional)')
            parser.add_argument('--C_db', metavar="<ALT_C_SEQS>", nargs='?',
                                help='specify alternative fasta file (if other than the one used to make recombinomes) containing C gene sequences for creation of BLAST database (optional)')
            
            args = parser.parse_args(sys.argv[2:])
            
            self.ncores = args.ncores
            self.force_overwrite = args.force_overwrite
            self.species = args.species
            self.receptor_name = args.receptor_name
            self.locus_name = args.locus_name
            self.N_padding = args.N_padding
            
            self.raw_seq_files = {}
            self.raw_seq_files['V'] = args.V_seqs
            self.raw_seq_files['J'] = args.J_seqs
            self.raw_seq_files['C'] = args.C_seqs
            self.prod_colour = args.colour
            if args.D_seqs:
                self.raw_seq_files['D'] = args.D_seqs
            if args.C_db:
                self.raw_seq_files['c'] = args.C_db
            config_file = args.config_file
            
        else:
            self.ncores = kwargs.get('ncores')
            self.force_overwrite = kwargs.get('force_overwrite')
            self.species = kwargs.get('species')
            self.receptor_name = kwargs.get('receptor_name')
            self.locus_name = kwargs.get('locus_name')
            self.N_padding = kwargs.get('N_padding')
            self.raw_seq_files = {}
            self.raw_seq_files['V'] = kwargs.get('V_seqs')
            self.raw_seq_files['J'] = kwargs.get('J_seqs')
            self.raw_seq_files['C'] = kwargs.get('C_seqs')
            self.prod_colour = kwargs.get('colour')
            if kwargs.get('D_seqs'):
                self.raw_seq_files['D'] = kwargs.get('D_seqs')
            if kwargs.get('C_db'):
                self.raw_seq_files['c'] = kwargs.get('C_db')

            config_file = kwargs.get('config_file')

        self.config = self.read_config(config_file)
        self.species_dir = self.get_resources_root(self.species)
        
        

    def run(self):

        # Check that there will not be git conflicts with inbuilt species
        #assert self.species not in ('Mmus', 'Hsap'), \
        #    "Cannot overwrite inbuilt species. Please choose a unique name " \
        #    "e.g. 'Mmus_1'"
        #
        self.init_dirs()
        
        self.calculate_colours(self.prod_colour)
        VDJC_files = self.copy_raw_files()
        recombinome_fasta = self.make_recombinomes(VDJC_files)
        self.make_bowtie2_index(recombinome_fasta)
        missing_dbs = self.make_igblast_db(VDJC_files)
        for s in missing_dbs:
            print("\nIMPORTANT: there is no IgBLAST database for {receptor}_{segment}\n"\
                  "Run build with {segment} segments for {receptor} before using tracer assemble\n"
                  .format(receptor = self.receptor_name, segment=s))
    
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
                                                            receptor_name = self.receptor_name)     
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
                if (self.receptor_name in colour_map and 
                   self.locus_name in colour_map[self.receptor_name]
                   and not colour_map[self.receptor_name][self.locus_name][0] == c):
                        if c in used_colours and c not in allowed_pal:
                            msg = "{c} already in use. Please specify a different colour.".format(c=c)
                            raise argparse.ArgumentTypeError(msg)
            
            prod_colour = c
        
        prod_rgb = hex2color(prod_colour)
        h, s, v = colorsys.rgb_to_hsv(*prod_rgb)
        
        nonprod_rgb = colorsys.hsv_to_rgb(h, s*0.5, self.adj_v(v))
        nonprod_colour = str(rgb2hex(nonprod_rgb))
        
        d1 = {self.locus_name : (prod_colour, nonprod_colour)}
        
        if self.receptor_name in colour_map:
            colour_map[self.receptor_name].update(d1)
        else:
            colour_map[self.receptor_name] = d1
        
        io.write_colour_file(colour_file, colour_map)
        
                
        
    def adj_v(self, v):
        new_v = v * 1.3
        if new_v > 1:
            new_v = 1
        return(new_v)

    def check_duplicate(self, new_path, segment=None, descriptor="Resource"):
        error_string = "{descriptor} already exists for {receptor}_{locus}" \
            .format(descriptor=descriptor, receptor=self.receptor_name,
                    locus=self.locus_name)
        if segment:
            error_string += "_" + segment
        error_string += ". Use --force_overwrite to replace existing file."
        if os.path.isfile(new_path):
            assert self.force_overwrite, error_string

    def init_dirs(self):

        # Set up output directories
        subdirs = ['igblast_dbs', 'combinatorial_recombinomes', 'raw_seqs']

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
            fn = "{receptor}_{locus}_{s}.fa".format(receptor=self.receptor_name,
                                                    locus=self.locus_name, s=s)
            out_file = os.path.join(self.species_dir, 'raw_seqs', fn)
            VDJC_files[s] = out_file
            self.check_duplicate(out_file, segment=s,
                                 descriptor="Sequence File")
            shutil.copy(self.raw_seq_files[s], out_file)

        return VDJC_files
     
    def load_segment_seqs(self, filename):
        seqs = {}
        with open(filename, 'rU') as fn:
            for record in SeqIO.parse(fn, 'fasta'):
                seqs[record.id] = str(record.seq)

        return seqs
        
    def make_recombinomes(self, VDJC_files):
        
        out_fasta = os.path.join(
            self.species_dir, 'combinatorial_recombinomes',
            '{receptor}_{locus}.fa'.format(receptor=self.receptor_name,
                                           locus=self.locus_name))

        self.check_duplicate(out_fasta, descriptor="Combinatorial recombinome")
        
        seqs = {}
        
        for s in 'VJC':
            in_file = VDJC_files[s]
            seqs[s] = self.load_segment_seqs(in_file)

        """# Logical check for C region
        if len(seqs['C']) > 1:
            print("\nMore than one constant region sequence included in {C_file}." \
                  .format(self.raw_seq_files['C']))
            print("Please only provide one constant sequence.\n")
            sys.exit(1)"""

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
                    seq = N_leader_string + V_seq.lower() + N_junction_string + \
                      J_seq.lower() + const_seq
                    seqs_to_write.append("{chr_name}\n{seq}\n"
                                     .format(seq=seq, chr_name=chr_name))
        
        elif len(seqs['C']) > 1:
            for V_name, V_seq in six.iteritems(seqs['V']):
                for J_name, J_seq in six.iteritems(seqs['J']):
                    for C_name, C_seq in six.iteritems(seqs['C']):
                        chr_name = ">chr={V_name}_{J_name}_{C_name}_C_region".format(J_name=J_name, V_name=V_name, C_name=C_name)
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
            '{receptor}_{locus}'.format(receptor=self.receptor_name,
                                        locus=self.locus_name))

        self.check_duplicate(index_base + ".1.bt2", descriptor="Bowtie2 index")
        
        command = [bowtie2_build, '-q', recombinome_fasta, index_base]
        try:
            subprocess.check_call(command)
        except subprocess.CalledProcessError:
            print("bowtie2-build failed")
    
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
            fn = "{receptor}_{segment}.fa".format(receptor=self.receptor_name,
                                                  segment=s)
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
                          '{receptor}_{segment} already found in {file}.'
                          .format(receptor=self.receptor_name, segment=s,
                                  file=fasta_file))
                    print('These sequences were not overwritten. '
                          'Use --force_overwrite to replace with new ones')
                    for seq in non_overwritten_seqs:
                        print(seq)
                
                command = [makeblastdb, '-parse_seqids', '-dbtype', 'nucl',
                           '-in', fasta_file]
                try:
                    subprocess.check_call(command)
                except subprocess.CalledProcessError:
                    print("makeblastdb failed for {receptor}_{segment}"
                          .format(receptor=self.receptor_name, segment=s))
