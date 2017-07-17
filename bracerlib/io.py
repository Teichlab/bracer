from __future__ import print_function

import csv
import distutils
import shutil

import pandas as pd
import os
import re
from collections import defaultdict

import six
import sys
from Bio import SeqIO

from bracerlib.bracer_func import process_chunk, find_possible_alignments, extract_blast_info
import glob
import pdb

import json

def makeOutputDir(output_dir_path):
    if not os.path.exists(output_dir_path):
        os.mkdir(output_dir_path)


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def clean_file_list(file_list):
    return_list = []
    trinity_pattern = re.compile(r"(.+)\.Trinity\.fasta")
    for file_name in file_list:
        clean_name = os.path.split(file_name)[1]
        trinity_match = trinity_pattern.search(clean_name)
        if trinity_match:
            clean_name = trinity_match.group(1)
        return_list.append(clean_name)

    return (sorted(return_list))


def get_filename_and_locus(name):
    name = name.split("_")
    cell = name[0]
    locus = "_".join(name[1:3])
    return ([cell, locus])


def sort_locus_names(dictionary_to_sort):
    for key, value in six.iteritems(dictionary_to_sort):
        sorted_value = sorted(value)
        dictionary_to_sort[key] = sorted_value
    return (dictionary_to_sort)


def load_IMGT_seqs(file):
    seqs = {}
    with open(file, 'rU') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seqs[record.id] = str(record.seq)
    return (seqs)


def parse_assembled_file(output_dir, cell_name, assembled_file):
    """Creates unique sequence names for each sequence if already assembled
    sequences are provided as input to BraCeR Assemble"""

    if os.path.exists(assembled_file):
        outfile = "{output_dir}/Trinity_output/{cell_name}.fasta".format(
                   output_dir=output_dir, cell_name=cell_name)
    write = False
    count = 0    
    with open(assembled_file, "r") as input:
        with open(outfile, "w") as output:
            for line in input:
                if line.startswith(">"):
                    write = True
                    header = ">TRINITY_DN0_c0_g0_i{}\n".format(str(count))
                    output.write(header)
                    count += 1
                     
                elif write == True:
                    if len(line) > 0:
                        output.write(line)


def parse_IgBLAST(loci, output_dir, cell_name, raw_seq_dir, species, 
                  assembled_file, max_junc_len=100):
    
    IMGT_seqs = dict()
    loci_for_segments = defaultdict(list)
    

    for locus in loci:
        seq_files = glob.glob(os.path.join(raw_seq_dir, "BCR_{}_*.fa".format(locus)))
        for f in seq_files:
            segment_name = os.path.splitext(os.path.split(f)[1])[0]
            IMGT_seqs[segment_name] = load_IMGT_seqs(f)
            loci_for_segments[segment_name.split("_")[2]].append(locus)
                    
    
    locus_names = ["_".join(["BCR",x]) for x in loci]
    all_locus_data = defaultdict(dict)

    for locus in locus_names:
        if assembled_file is not None:
            file = "{output_dir}/IgBLAST_output/{cell_name}.IgBLASTOut".format(
                                    output_dir=output_dir, cell_name=cell_name)
        else:    
            file = "{output_dir}/IgBLAST_output/{cell_name}_{locus}.IgBLASTOut".format(
                                output_dir=output_dir, cell_name=cell_name, locus=locus)
        if os.path.isfile(file):
            igblast_result_chunks = split_igblast_file(file)
         
            for chunk in igblast_result_chunks:
                (query_name, chunk_details) = process_chunk(chunk)
                if query_name is not None: 
                    if "cell_id=" in query_name:
                        query_name = query_name.split("=")[1]
                all_locus_data[locus][query_name] = chunk_details
                
        else:
            all_locus_data[locus] = None

    cell = find_possible_alignments(all_locus_data, locus_names, cell_name, IMGT_seqs, 
                                    output_dir, species, loci_for_segments, 
                                    loci, max_junc_len, assembled_file)
    return (cell)


def parse_BLAST(loci, output_dir, cell_name, species, assembled_file):
    """Parses BLAST output from output files and writes formatted output to BLAST 
    output summary files"""

    locus_names = ["_".join(["BCR",x]) for x in loci]    

    for locus in loci:
        
        blast_dir = "BLAST_output"
        output_file = "{outdir}/{blast_dir}/blastsummary_{locus}.txt".format(
                       outdir=output_dir, blast_dir=blast_dir, locus=locus)
        input_file = "{output_dir}/{blast_dir}/{cell_name}_BCR_{locus}.xml".format(
                      output_dir=output_dir, blast_dir=blast_dir, cell_name=cell_name, 
                      locus=locus)

        with open(output_file, 'w') as outfile:
            outfile.write("------------------\n##{}##\n------------------\n\n#BCR_{}#\n\n".format(
                                                                    cell_name, locus))
        
            # Split result file into chunks corresponding to results for each query sequence.
            if os.path.isfile(input_file):
                blast_result_chunks = split_blast_file(input_file)

                for chunk in blast_result_chunks:
                    message = False
                    for line_x in chunk:
                        line_x= line_x.strip()

                        if line_x.startswith("<Iteration_query-def>"):
                            line = line_x.split(">")[1]
                            blast_query_name = line.split("<")[0]
                            if assembled_file == True or " " in blast_query_name:
                                blast_query_name = blast_query_name.split()[0]
                            
                        elif line_x.startswith("<Hsp_evalue>"):
                            evalue = extract_blast_info(line_x)
                            evalue = format(float(evalue), '.0e')

                        elif line_x.startswith("<Hit_accession>"):
                            C_segment = extract_blast_info(line_x)
                            if "C-REGION" or "CH1" in C_segment:
                                C_segment = C_segment.split("_")[0]
                        elif line_x.startswith("<Hsp_bit-score>"):
                            bit_score = extract_blast_info(line_x)
                              
                        elif line_x.startswith("<Hsp_query-from>"):
                            q_start = extract_blast_info(line_x)
                           
                        elif line_x.startswith("<Hsp_query-to>"):
                            q_end = extract_blast_info(line_x)
                         
                        elif line_x.startswith("<Hsp_hit-from>"):
                            s_start = extract_blast_info(line_x)
                            
                        elif line_x.startswith("<Hsp_hit-to>"):
                            s_end = extract_blast_info(line_x)
                            
                        elif line_x.startswith("<Iteration_query-len>"):
                            query_length = extract_blast_info(line_x)
                            
                        elif line_x.startswith("<Hsp_align-len>"):
                            align_length = extract_blast_info(line_x)
                            
                        elif line_x.startswith("<Hsp_gaps>"):
                            gaps = extract_blast_info(line_x)

                        elif line_x.startswith("<Hsp_identity>"):
                            identity = extract_blast_info(line_x)

                        elif line_x.startswith("<Iteration_message>No hits found"):
                            message = True
                            out_string = "##{blast_query_name}##\nNo C segment found\n\n".format(
                                                                blast_query_name=blast_query_name)
                            outfile.write(out_string)
                        
                        # Create output string when reaching end of BLAST 
                        # iteration result (marked by </Iteration>) and write 
                        # to BLAST summary file
                        elif line_x.startswith("</Iteration>") and message is not True:
                            identity_pro = float(identity)/int(align_length)*100
                            identity_pro = format(identity_pro, '.2f')
                            mismatches = int(align_length) - int(identity)

                            #Account for reversed sequences
                            if int(s_start) > int(s_end):
                                blast_query_name = "reversed|" + blast_query_name
                                x, y = int(q_start), int(q_end)
                                q_start = int(query_length) - y + 1
                                q_end = int(query_length) - x + 1
                                s_start, s_end = s_end, s_start
                               
                            intro_string = "##{}##\nC segment:\t{}\n\n".format(
                                            blast_query_name, C_segment)
                            header_string = ("Segment\tquery_id\tsubject_id\t% identity\talignment length\t"
                                            "mismatches\tgap opens\tgaps\tq start\tq end\ts start\ts end\t"
                                            "evalue\tbit score\n")
                            out_string = ("C\t{blast_query_name}\t{C_segment}\t{identity_pro}\t{align_length}\t{mismatches}\tNA\t{gaps}\t{q_start}\t{q_end}\t{s_start}\t{s_end}\t{evalue}\t{bit_score}\n\n").format(
                                            blast_query_name=blast_query_name, 
                                            C_segment=C_segment, identity_pro=identity_pro, align_length=align_length, 
                                            evalue=evalue, mismatches=mismatches, gaps=gaps, q_start=q_start, 
                                            q_end=q_end, s_start=s_start, s_end=s_end, bit_score=bit_score)
                            string_to_write = intro_string + header_string + out_string
                            outfile.write(string_to_write)                           
                                  

def split_igblast_file(filename):
    # code adapted from http://stackoverflow.com/questions/19575702/pythonhow-to-split-file-into-chunks-by-the-
    #occurrence-of-the-header-word
    token = '# IGBLASTN'
    chunks = []
    current_chunk = []

    with open(filename) as fh:
        for line in fh:
            line = line.rstrip()
                
            if line.startswith(token) and current_chunk and not line.startswith("Total queries"):
                # if line starts with token and the current chunk is not empty
                chunks.append(current_chunk[:])  # add not empty chunk to chunks
                current_chunk = []  # make current chunk blank
            # just append a line to the current chunk on each iteration
            if not line.startswith("Total queries"):
                current_chunk.append(line)

        chunks.append(current_chunk)  # append the last chunk outside the loop
    return (chunks)


def split_blast_file(filename):
    # code adapted from http://stackoverflow.com/questions/19575702/pythonhow-to-split-file-into-chunks-by-the-
    #occurrence-of-the-header-word
    token = '<Iteration>'
    chunks = []
    current_chunk = []

    with open(filename) as fh:
        for line in fh:
            line = line.rstrip()

            if line.startswith(token) and current_chunk:
                chunks.append(current_chunk[:])
                current_chunk = []  
            if not line.startswith("Total queries"):
                current_chunk.append(line)

        chunks.append(current_chunk)  
    return (chunks)


def check_binary(name, user_path=None):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    if user_path:
        if is_exe(user_path):
            return user_path
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, name)
            if is_exe(exe_file):
                return exe_file

    raise OSError("Required binary not found: {name}. Please add to PATH or specify location in config file."
                  .format(name=name))


def read_colour_file(filename, return_used_list=False, receptor_name=None):
    
    colour_map = dict()
    used_colours = set()
    with open(filename) as f:
        for line in f:
            line = line.rstrip()
            receptor, locus, prod_colour, nonprod_colour = line.split(",")
            d = {locus : (prod_colour, nonprod_colour)}
            
            if receptor in colour_map:
                colour_map[receptor].update(d)
            else:
                colour_map[receptor] = d
            if receptor_name is not None and receptor == receptor_name:
                used_colours.add(prod_colour)
            elif receptor_name is None:
                used_colours.add(prod_colour)
    if return_used_list:
        t = (colour_map, used_colours)
        return t
    else:
        return colour_map
            
def write_colour_file(filename, colour_map):
    sorted_receptors = sorted(colour_map.keys())
    with open(filename, 'w') as f:
        for receptor in sorted_receptors:
            sorted_loci = sorted(colour_map[receptor].keys())
            for l in sorted_loci:
                colours = colour_map[receptor][l]
                f.write("{},{},{},{}\n".format(receptor, l, colours[0], colours[1]))
