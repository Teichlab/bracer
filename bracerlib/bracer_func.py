##############################################################################

#                         Functions for use with                             #
# BraCeR - a tool to reconstruct BCR sequences from single-cell RNA-seq data #    
#                                                                            #
# Please see README and LICENCE for details of use and licence conditions.   #
# This software was written by Mike Stubbington (ms31@sanger.ac.uk) from the #
# Teichmann Lab, EMBL-EBI and WTSI (www.teichlab.org). Latest versions are   #
# available for download at www.github.com/teichlab/bracer.                  #
#                                                                            #
#      Copyright (c) 2015, 2016 EMBL - European Bioinformatics Institute     #
#      Copyright (c) 2016 Genome Research Ltd.                               #
#      Author: M.J.T. Stubbington, I. Lindeman ida.lindeman@medisin.uio.no   #
##############################################################################

from __future__ import print_function

import glob
import os
import re
import shutil
import subprocess
from collections import defaultdict, Counter
import csv
from time import sleep

import Levenshtein
import networkx as nx
import six
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

from bracerlib.core import Cell, Recombinant
import bracerlib.io

import copy

import pdb
import gzip


def hamming_dist(str1, str2):
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

def process_chunk(chunk):
    store_VDJ_rearrangement_summary = False
    store_junction_details = False
    store_alignment_summary = False
    store_hit_table = False
    alignment_summary = []
    hit_table = []
    looking_for_end = False
    return_dict = defaultdict(list)
    query_name = None
    for line_x in chunk:
        

        if store_VDJ_rearrangement_summary:
            VDJ_rearrangement_summary = line_x.split("\t")
            for i in VDJ_rearrangement_summary:
                return_dict['VDJ_rearrangement_summary'].append(i)
            store_VDJ_rearrangement_summary = False

        elif store_junction_details:
            junction_details = line_x.split("\t")
            for i in junction_details:
                return_dict["junction_details"].append(i)
            store_junction_details = False

        elif store_alignment_summary:
            if not line_x.startswith("#"):
                if line_x.startswith("Total"):
                    store_alignment_summary = False
                else:
                    return_dict['alignment_summary'].append(line_x)

        elif store_hit_table:
            if not looking_for_end:
                if not line_x.startswith("#"):
                    return_dict['hit_table'].append(line_x)
                    looking_for_end = True
            else:
                if line_x.startswith("#") or line_x.startswith("\n"):
                    store_hit_table = False
                else:
                    return_dict['hit_table'].append(line_x)

        elif line_x.startswith('# Query'):
            query_name = line_x.split()[2]
            if "cell_id=" in query_name:
                query_name = query_name.split("=")[1]
            try:
                query_length = query_name.split("Length_")[1]
                return_dict['query_length'] = int(query_length)
            except:
                try:
                    query_length = line_x.split()[3]
                    return_dict['query_length'] = int(query_length.split("=")[1])
                except:
                    return_dict['query_length'] = None

        elif line_x.startswith('# V-(D)-J rearrangement summary'):
            store_VDJ_rearrangement_summary = True

        elif line_x.startswith('# V-(D)-J junction details'):
            store_junction_details = True

        elif line_x.startswith('# Alignment summary'):
            store_alignment_summary = True

        elif line_x.startswith('# Hit table'):
            store_hit_table = True
    return (query_name, return_dict)


      
def extract_blast_info(line):
    line = line.split()[0]
    info = line.split(">")[1]
    info = info.split("<")[0]
    return (info)

def get_C_gene(output_dir, locus, contig_name):
    blast_dir = "BLAST_output"

    locus = locus.split("_")[1]
    blast_summary_file = "{output_dir}/{blast_dir}/blastsummary_{locus}.txt".format(
                            output_dir=output_dir, blast_dir=blast_dir, locus=locus)

    store_details = False
    C_gene = None
    info_line = None
    start_position = None

    with open(blast_summary_file, 'r') as input:
        for line in input:

            if line.startswith("C\t{contig_name}".format(
                contig_name=contig_name)) or line.startswith("C\treversed|{contig_name}".format(
                                                            contig_name=contig_name)):
                C_gene = line.split("\t")[2]

                if "_CH1" or "_C-REGION" in C_gene:
                    C_gene = C_gene.split("_")[0]
                info_line = line
                start_position = int(line.split("\t")[8])

    return (C_gene, info_line, start_position)

 
def parse_alignment_summary(alignment_summary):
    start = None
    cdr3_start = None
    for entry in alignment_summary:
        info = entry.split("\t")
        region = info[0]
        if not region.startswith("FR3-IMGT"):
            if start == None:
                start = int(info[1])
        elif region.startswith("FR3-IMGT"):
             cdr3_start = int(info[2]) + 1
            
    return (start, cdr3_start)    


def get_gapped_db_string(locus, query_name, cell_name, output_dir):
    """Get entry from IgBlast gapped output file"""
    gapped_db_file = "{}/IgBLAST_output/{}_{}_db-pass.tab".format(
                        output_dir, cell_name, locus)
    gapped_db_string = ""
    if os.path.isfile(gapped_db_file):
        with open(gapped_db_file, "r") as input:
            for line in input:
                if line.startswith("TRINITY"):
                    info = line.split("\t")
                    if info[0].split()[0] == query_name:
                        gapped_db_string = line
                        break
    return (gapped_db_string)

def get_gapped_db_header(locus, cell_name, output_dir):
    """Get header from IgBlast gapped output file to get order of information
    which may vary between IgBlast versions"""
    gapped_db_file = "{}/IgBLAST_output/{}_{}_db-pass.tab".format(
                            output_dir, cell_name, locus)
    if os.path.isfile(gapped_db_file):
        with open(gapped_db_file, "r") as input:
            for line in input:
                gapped_db_header = line.split("\t")
                break
    return (gapped_db_header)

def parse_gapped_db_string(gapped_db_string, gapped_db_header):
    """Extracts information from IgBlast gapped output using the headers"""
    n = 0
    info = gapped_db_string.split("\t")
    for header in gapped_db_header:
        header = header.strip()
        if header == "FUNCTIONAL":
            productive = info[n].strip()
        elif header == "IN_FRAME":
            in_frame = info[n].strip()
        elif header == "STOP":
            stop = info[n].strip()
        elif header == "INDELS":
            indels = info[n].strip()
        elif header == "JUNCTION":
            junction = info[n].strip()
        elif header == "CDR3_IMGT":
            cdr3_seq = info[n].strip()
        n +=1

    # Translate "F/T" to True/False
    if productive == "T":
        productive = True
    else:
        productive = False
    if in_frame == "T":
        in_frame = True
    else:
        in_frame = False
    if stop == "T":
        stop = True
    else:
        stop = False
    if indels == "T":
        indels = True
    else:
        indels = False

    return (productive, in_frame, stop, indels, junction, cdr3_seq)



def find_possible_alignments(sample_dict, locus_names, cell_name, IMGT_seqs, 
                output_dir, species, loci_for_segments, loci,
                max_junc_string_length, assembled_file):


    alignment_dict = defaultdict(dict)
    recombinants = {}

    for locus in locus_names:
        recombinants[locus] = []
        data_for_locus = sample_dict[locus]

        if data_for_locus is not None:
            for query_name, query_data in six.iteritems(data_for_locus):
                processed_hit_table = process_hit_table(query_name, query_data, 
                                                                         locus)
                if processed_hit_table is not None:
                    (returned_locus, good_hits, rearrangement_summary) = \
                                                        processed_hit_table
                    junction_list = query_data['junction_details']

                    best_V = remove_allele_stars(
                        rearrangement_summary[0].split(",")[0])

                    junc_string = "".join(junction_list)
                    junc_string = remove_NA(junc_string)
                    junc_string= remove_parentheses(junc_string)
                    
                    locus_letter = returned_locus.split("_")[1]
                    
                    if locus_letter in loci_for_segments['D']:
                        has_D = True
                    else:
                        has_D = False
                    
                    if has_D:
                        best_J = remove_allele_stars(
                                 rearrangement_summary[2].split(",")[0])
                    else:
                        best_J = remove_allele_stars(
                                 rearrangement_summary[1].split(",")[0])

                    # Get junctional string from gapped_db_string
                    gapped_db_string = get_gapped_db_string(locus, query_name, cell_name, output_dir)
                    try:
                        gapped_db_header = get_gapped_db_header(locus, cell_name, output_dir)
                    except:
                        gapped_db_header = None
                    if "\t" in gapped_db_string:
                        (productive, in_frame, stop, indels, junc_string, cdr3_seq) = \
                                parse_gapped_db_string(gapped_db_string, gapped_db_header)
                    

                    identifier = best_V + "_" + junc_string + "_" + best_J

                    ##line attempting to add alignment summary to data for 
                    #use with PCR comparisons

                    alignment_summary = query_data['alignment_summary']
                    (C_gene, C_info_line) = (None, None)

                    align_start, cdr3_start = parse_alignment_summary(alignment_summary)
                    (C_gene, C_info_line, C_start) = get_C_gene(output_dir, 
                                                        locus, query_name)
                    
                    if has_D:
                        all_J_names = [remove_allele_stars(x) for x in 
                        rearrangement_summary[2].split(',')]
                    else:
                        all_J_names = [remove_allele_stars(x) for x in 
                        rearrangement_summary[1].split(',')]


                    # get original sequence from Trinity file
                    # needed for summary of reconstructed lengths.
                    # Only use the VDJ portion found by IgBLAST or trim for BCRs:
                    if assembled_file is not None:
                        trinity_file = "{}/Trinity_output/{}.fasta".format(
                                                    output_dir, cell_name)
                    else:
                        trinity_file = "{output_dir}/Trinity_output/{cell_name}_{locus}.Trinity.fasta".format(
                                locus=locus, output_dir=output_dir, cell_name=cell_name)

                    with open(trinity_file, 'rU') as tf:
                        for record in SeqIO.parse(tf, 'fasta'):
                            if query_name in record.id:
                                trinity_seq = record

                    query_length = query_data["query_length"]
                    if query_length == None:
                        query_length = len(trinity_seq)
                    if 'reversed' in good_hits[0][1]:
                        trinity_seq = str(trinity_seq.reverse_complement().seq)
                    else:
                        trinity_seq = str(trinity_seq.seq)

                    start_coord, end_coord = get_coords(good_hits)

                    if C_gene is not None:
                        if (C_start - 1) > end_coord:
                            end_coord = C_start - 1
                            

                    #Assess if rearrangement is full-length 
                    #(from start of V gene to start of C gene)
                    full_length = is_rearrangement_full_length(trinity_seq, 
                                  query_data["hit_table"], query_name, 
                                  query_length, output_dir, locus)
                    untrimmed_seq = trinity_seq
                    trinity_seq = trinity_seq[start_coord:end_coord]
                    cdr3_seq = None

                    # Get recombinant info from ungapped sequences
                    fasta_line_for_contig = trinity_seq
                    (is_productive, bestVJNames, cdr3, cdr3_seq) = (
                        get_fasta_line_for_contig_assembly(trinity_seq, good_hits,
                        returned_locus, IMGT_seqs, cell_name, query_name, 
                        loci_for_segments, full_length, alignment_summary, 
                        rearrangement_summary))
                    productive = is_productive[0]
                    stop_codon = is_productive[1]
                    in_frame = is_productive[2]

                    # Replace with productivity info from IMGT-gapped sequences
                    # if available
                    if "\t" in gapped_db_string:
                        (productive, in_frame, stop_codon, indels, junc_string, cdr3_seq) = \
                            parse_gapped_db_string(gapped_db_string, gapped_db_header)
                        if "-" in cdr3_seq:
                            new_cdr3_seq = ""
                            for l in cdr3_seq:
                                if not l == "-":
                                    new_cdr3_seq += l
                            cdr3_seq = new_cdr3_seq
                        cdr3 = Seq(str(cdr3_seq), generic_dna).translate()

                    #Identify the most likely V and J genes
                    if locus in ["H", "BCR_H"]:
                        threshold_percent = 0.04
                    else:
                        threshold_percent = 0.01
                    all_V_names = find_V_genes_based_on_bit_score(trinity_seq, 
                        query_data["hit_table"], query_name, threshold_percent)
                    all_J_names = find_J_genes_based_on_bit_score(trinity_seq, 
                        query_data["hit_table"], query_name, threshold_percent)
                                      
                    V_genes = all_V_names
                    J_genes = all_J_names
                    
                    all_poss_identifiers = set()
                    for V in all_V_names:
                        for J in all_J_names:    
                            i = V + "_" + J
                            all_poss_identifiers.add(i)
                    

                    if len(junc_string) < int(max_junc_string_length):
                        rec = Recombinant(contig_name=query_name, 
                                locus=returned_locus, identifier=identifier, 
                                all_poss_identifiers=all_poss_identifiers, 
                                productive=productive, stop_codon=stop_codon, 
                                in_frame=in_frame, TPM=0.0, 
                                dna_seq=fasta_line_for_contig, hit_table=good_hits, 
                                summary=rearrangement_summary, 
                                junction_details=junction_list, 
                                best_VJ_names=bestVJNames, 
                                alignment_summary=alignment_summary, 
                                trinity_seq=trinity_seq, 
                                has_D=has_D, output_dir=output_dir, 
                                full_length=full_length, query_length=query_length, 
                                V_genes=V_genes, J_genes=J_genes, cdr3=cdr3, 
                                C_gene=C_gene, C_info_line=C_info_line, 
                                cdr3_seq=cdr3_seq, junc_string=junc_string,
                                untrimmed_seq=untrimmed_seq, 
                                gapped_db_string=gapped_db_string, gapped_db_header=gapped_db_header)
                        recombinants[locus].append(rec)

    if recombinants:
        for locus, rs in six.iteritems(recombinants):
            # Collapse sequences with identical CDR3 nt sequence
            # or CDR3 nt/junctional Hamming distance <= 0.07 normalised by length
            # and overlapping V- and J segment assignments. 
            recombinants[locus] = collapse_close_sequences(rs, locus)
        cell = Cell(cell_name, recombinants, species=species, loci=loci)
        
    else:
        cell = Cell(cell_name, None, species=species, loci=loci)

    
    return (cell)



def find_V_genes_based_on_bit_score(seq, hit_table, query_name, threshold_percent):
    """Identifies the most likely V genes based on IgBlast bit scores"""
    found_V = False
    V_genes = []
    threshold = None
    for hit in hit_table:
        if hit == "":
            continue
        info = hit.split()
        segment = info[0]
        allele = info[2]
        V_gene = allele.split("*")[0]
        bit_score = float(info[13])
        if segment == "V":
            if found_V == False:
                top_bit_score = bit_score
                found_V = True
                threshold = bit_score - bit_score*threshold_percent 
                V_genes.append(V_gene)
            elif found_V == True:
                if bit_score >= threshold and V_gene not in V_genes:
                    V_genes.append(V_gene)
    return(V_genes)
        

def find_J_genes_based_on_bit_score(seq, hit_table, query_name, threshold_percent):
    """Identifies the most likely J genes based on IgBlast bit scores"""
    found_J = False
    J_genes = []
    threshold = None
    for hit in hit_table:
        if hit == "":
            continue
        info = hit.split()
        segment = info[0]
        allele = info[2]
        J_gene = allele.split("*")[0]
        bit_score = float(info[13])
        if segment == "J":
            if found_J == False:
                top_bit_score = bit_score
                found_J = True
                threshold = bit_score - bit_score*threshold_percent
                J_genes.append(J_gene)
            elif found_J == True:
                if bit_score >= threshold and J_gene not in J_genes:
                    J_genes.append(J_gene)
    return(J_genes)


def parse_rearrangement_summary(rearrangement_summary):
    """Returns a tuple of (stop_codon, in_frame, productive) from IgBlast output.
    Only used if no info can be parsed from IMGT-gapped IgBlast output by Change-O"""
    i = 1
    chain_type = rearrangement_summary[0][0:4]
    if "IGH" in chain_type:
        i = 0

    stop_codon = rearrangement_summary[4-i]
    in_frame = rearrangement_summary[5-i]
    productive = rearrangement_summary[6-i]

    if stop_codon == "No":
        stop_codon = False
    elif stop_codon == "Yes":
        stop_codon = True
    else:
        stop_codon = None

    if "In" in in_frame:
        in_frame = True
    else:
        in_frame = False
    
    if productive == "No":
        productive = False
    elif productive == "Yes":
        productive = True
    else:
        productive = None
   
    return (stop_codon, in_frame, productive)

def get_coords(hit_table):
    """Gets start and end coordinates for V(D)J sequence determined by IgBlast"""
    found_V = False
    found_J = False
    for entry in hit_table:
        if entry == "":
            continue
        if entry[0] == 'V':
            if not found_V:
                start = int(entry[8]) - 1
                found_V = True
        if entry[0] == 'J':
            if not found_J:
                end = int(entry[9])
                found_J = True
                bitscore = float(entry[13])
            # Look for other J gene hits with equal bitscore but higher J gene 
            #end position        
            if found_J:
                next_bitscore = float(entry[13])
                if next_bitscore == bitscore:
                    new_end = int(entry[9])
                    if new_end > end:
                        end = new_end
    return (start, end)


def remove_NA(junc_string):
    new_string = junc_string.replace("N/A", "")
    return (new_string)


def remove_parentheses(junc_string):
    new_string = junc_string.replace("(", "")
    new_string = new_string.replace(")", "")
    return (new_string)


def remove_allele_stars(segment):
    p = re.compile(r"(.+)\*\d+")
    m = p.search(segment)
    return (m.group(1))


def process_hit_table(query_name, query_data, locus):
    hit_table = query_data['hit_table']
    rearrangement_summary = query_data['VDJ_rearrangement_summary']


    found_V = set()
    found_D = set()
    found_J = set()

    good_hits = []

    locus_name = locus.split("_")[1]
    
    for entry in hit_table:
        if not entry == "":
          
            entry = entry.split("\t")
            segment = entry[2]
            segment_locus = segment[2]
            segment_type = segment[3]
            e_value = float(entry[12])
            if locus in ["H", "BCR_H"] and segment_type == "V":
                e_value_cutoff = 5e-4
            else:
                e_value_cutoff = 5e-3
            if locus_name in segment_locus:
                if e_value < e_value_cutoff:
                    if segment_type == "V":
                        found_V.add(locus)
                        good_hits.append(entry)
                    elif segment_type == "J":
                        found_J.add(locus)
                        good_hits.append(entry)
                else:
                    if segment_type == "D":
                        percent_identity = float(entry[3])
                        if percent_identity == 100:
                            found_D.add(locus)
                            good_hits.append(entry)
                            
    if locus in found_V and locus in found_J:
        return (locus, good_hits, rearrangement_summary)
    else:
        return (None)



def is_rearrangement_full_length(seq, hit_table, query_name, query_length, 
                                                        output_dir, locus):
    found_V = False
    found_J = False
    full_5_prime = False
    full_3_prime = False
    ref_V_start = None
    J_end_pos = None
    V_hit = None
    J_hit = None
    (C_gene, C_info_line, C_position) = get_C_gene(output_dir, locus, query_name)
    del (C_info_line)
    del (C_position)
   
    full_length = False

    for hit in hit_table:
        if hit == "":
            continue
        info = hit.split()
        segment = info[0]
        if segment == "V" and found_V == False:
            ref_V_start = int(info[10])
            V_hit = hit
            found_V = True
        elif segment == "J" and found_J == False:
            J_end_pos = int(info[9])
            J_ref_end_pos = int(info[11])
            J_hit = hit
            found_J = True
     
        if ref_V_start == 1:
            full_5_prime = True
        if C_gene is not None:
            full_3_prime = True
        if full_5_prime and full_3_prime:
            full_length = True
        elif full_5_prime == True:
            if J_end_pos is not None:
                if ("HJ" in J_hit and J_ref_end_pos > 40) \
                    or (("KJ" or "LJ") in J_hit and J_ref_end_pos > 30):
                    if int(query_length)>= (J_end_pos - 1):
                        full_length = True
                    
    return (full_length)


def get_segment_name(name, pattern):
    match = pattern.search(name)
    number = match.group(1)
    if match.group(3):
        sub_number = match.group(3)
    else:
        sub_number = ""
    return (number)


def get_fasta_line_for_contig_assembly(trinity_seq, hit_table, locus, 
    IMGT_seqs, sample_name, query_name, loci_for_segments, full_length, 
    alignment_summary, rearrangement_summary):

    found_best_V = False
    found_best_D = False
    found_best_J = False

    V_pattern = re.compile(r".+{potential_loci}V.+".format(
                potential_loci='[' + "".join(loci_for_segments['V']) + ']'))
    D_pattern = re.compile(r".+{potential_loci}D.+".format(
                potential_loci='[' + "".join(loci_for_segments['D']) + ']'))
    J_pattern = re.compile(r".+{potential_loci}J.+".format(
                potential_loci='[' + "".join(loci_for_segments['J']) + ']'))

    for hit in hit_table:
        segment = hit[2]
        if V_pattern.search(segment) and not found_best_V:
            V_locus_key = V_locus_key = "_".join([locus, 'V'])
            best_V_name = segment
            ref_V_seq = IMGT_seqs[V_locus_key][segment]
            found_best_V = True
        elif J_pattern.search(segment) and not found_best_J:
            J_locus_key = "_".join([locus, 'J'])
            best_J_name = segment
            ref_J_seq = IMGT_seqs[J_locus_key][segment]
            found_best_J = True
    
    # work out if sequence that exists is in frame
    found_V = False
    found_J = False
    for entry in hit_table:
        if entry[0] == 'V':
            if not found_V:
                ref_V_start = int(entry[10])
                found_V = True
        if entry[0] == 'J':
            if not found_J:
                ref_J_end = int(entry[11])
                J_end = int(entry[9])
                found_J = True
                bitscore = float(entry[13])
            if found_J:
                next_bitscore = float(entry[13])
                if next_bitscore == bitscore:
                    new_J_end = int(entry[9])
                    if new_J_end > J_end:
                        ref_J_end = int(entry[11])

    start_padding = ref_V_start - 1
    ref_J_length = len(ref_J_seq)
    cdr3_length = None
    cdr3_seq = None
    

    (contains_stop, in_frame, productive) = parse_rearrangement_summary(
                                                    rearrangement_summary)
    (start, stop) = get_coords(hit_table)
    del (stop)
 
    (align_start, cdr3_start) = parse_alignment_summary(alignment_summary)

    if cdr3_start is None:
        cdr3 = "Couldn't find CDR3"
        in_frame = False
        
    else:
        cdr3_start = cdr3_start - start
        in_frame = True
    
    # Trim sequence to look for CDR3 in correct region
    if in_frame:
        seq = trinity_seq[cdr3_start - 7:]
        cdr3 = get_cdr3(seq, locus)
        if "Couldn't" in cdr3:
            # Look for out-of-frame cdr3
            in_frame = False
            cdr3 = None
            cdr3_length = None
            cdr3_seq = None
            dna_seq = trinity_seq[cdr3_start - 1:]
            cdr3_2, motif_2 = get_out_of_frame_cdr3(dna_seq, locus, 2)
            cdr3_3, motif_3 = get_out_of_frame_cdr3(dna_seq, locus, 3)
            if cdr3_2 is None and cdr3_3 is None:
                cdr3 = None
                cdr3_length = None
                cdr3_seq = None
            elif cdr3_2 and cdr3_3:
                if (motif_2 == motif_3) or ((motif_2 and motif_3) in [".G.G", "FSDG", "WSQG"]):
                    if len(cdr3_2) > len(cdr3_3):
                        cdr3 = cdr3_2
                        frame = 2
                    else:
                        cdr3 = cdr3_3
                        frame = 3
                elif motif_2 in ["WG.G", "FG.G"]:
                    cdr3 = cdr3_2
                    frame = 2
                elif motif_3 in ["WG.G", "FG.G"]:
                    cdr3 = cdr3_3
                    frame = 3
                    
            elif cdr3_2:
                cdr3 = cdr3_2
                frame = 2
            else:
                cdr3 = cdr3_3
                frame = 3
            if not cdr3 is None:
                cdr3_length = len(cdr3)
                cdr3_seq = trinity_seq[cdr3_start-1:cdr3_length*3+frame-1+cdr3_start-16]
        else:    
            cdr3_length = len(cdr3)
            # Get CDR3 nt sequence excluding conserved Cys and XGXG motif
            cdr3_seq = trinity_seq[cdr3_start-1:cdr3_length*3+cdr3_start-16]
            in_frame = is_cdr3_in_frame(cdr3, locus)

    if in_frame and not contains_stop:
        productive = True
    else:
        productive = False

    
    print()
    
    

    productive_rearrangement = (productive, contains_stop, in_frame)
    bestVJ = [best_V_name, best_J_name]
    

    return (productive_rearrangement, bestVJ, cdr3, cdr3_seq)

def get_out_of_frame_cdr3(dna_seq, locus, frame):
    dna_seq = dna_seq[frame-1:]
    aaseq = Seq(str(dna_seq), generic_dna).translate()
    if locus in ["BCR_H", "H"]:
        motif_start = "W"
    else:
        motif_start = "F"
    motif = motif_start + "G.G"
    if re.findall(motif, str(aaseq)):
        indices = [i for i, x in enumerate(aaseq)]
        if len(re.findall(motif, str(aaseq))) > 1 and str(aaseq).find(re.findall(motif, str(aaseq))[0]) < 15:
            upper = str(aaseq).find(re.findall(motif, str(aaseq))[1])
        else:
            upper = str(aaseq).find(re.findall(motif, str(aaseq))[0])
        cdr3 = aaseq[0:upper + 4]

    elif re.findall(".G.G", str(aaseq)):
        indices = [i for i, x in enumerate(aaseq)]
        upper = str(aaseq).find(re.findall(".G.G", str(aaseq))[0])
        try:
            upper2 = str(aaseq).find(re.findall(".G.G", str(aaseq))[1])
        except:
            upper2 = 0
        cdr3 = aaseq[0:upper + 4]
        motif = ".G.G"

    elif re.findall("WSQG", str(aaseq)) and locus in ["BCR_H", "H"]:
        indices = [i for i, x in enumerate(aaseq)]
        upper = str(aaseq).find(re.findall("WSQG", str(aaseq))[0])
        cdr3 = aaseq[0:upper + 4]
        motif = "WSQG"

    elif re.findall("FSDG", str(aaseq)) and locus in ["BCR_K", "K"]:
        indices = [i for i, x in enumerate(aaseq)]
        upper = str(aaseq).find(re.findall("FSDG", str(aaseq))[0])
        cdr3 = aaseq[0:upper + 4]
        motif = "FSDG"

    else:
        cdr3 = None
        motif = None

    return (cdr3, motif)

def get_cdr3(dna_seq, locus):
    aaseq = Seq(str(dna_seq), generic_dna).translate()
    # Specify first amino acid in conserved motif according to receptor and locus
    if locus in ["BCR_H", "H"]:
        motif_start = "W"
    else:
        motif_start = "F"
    motif = motif_start + "G.G"
    lower = False
    lower2 = False
    if re.findall(motif, str(aaseq)) and re.findall('C', str(aaseq)):
        indices = [i for i, x in enumerate(aaseq) if x == 'C']
        # Look for "false ends" of CDR3s - end motif in middle of CDR3 making the CDR3 very short
        if len(re.findall(motif, str(aaseq))) > 1 and str(aaseq).find(re.findall(motif, str(aaseq))[0]) < 15:
            upper = str(aaseq).find(re.findall(motif, str(aaseq))[1])
        else:
            upper = str(aaseq).find(re.findall(motif, str(aaseq))[0])
        for i in indices:
            if not locus in ["BCR_H", "H", "BCR_K", "K", "BCR_L", "L"]:
                if i < upper:
                    lower = i
            else:
                 if i < upper and lower == False:
                     lower = i
        if lower:
            cdr3 = aaseq[lower:upper + 4]
        else:
            cdr3 = "Couldn't find conserved cysteine"
    elif re.findall(".G.G", str(aaseq)) and re.findall('C', str(aaseq)):
        indices = [i for i, x in enumerate(aaseq) if x == 'C']
        upper = str(aaseq).find(re.findall(".G.G", str(aaseq))[0])
        try:
            upper2 = str(aaseq).find(re.findall(".G.G", str(aaseq))[1])
        except:
            upper2 = 0
        for i in indices:
            if i < upper:
                lower = i
            elif i < upper2:
                lower2 = i
        if lower:
            cdr3 = aaseq[lower:upper + 4]
        elif lower2:
            cdr3 = aaseq[lower2:upper2 + 4]
        else:
            cdr3 = "Couldn't find conserved cysteine"

    elif re.findall("WSQG", str(aaseq)) and re.findall('C', str(aaseq)):
        indices = [i for i, x in enumerate(aaseq) if x == 'C']
        upper = str(aaseq).find(re.findall("WSQG", str(aaseq))[0])
        for i in indices:
            if i < upper:
                lower = i
        if lower:
            cdr3 = aaseq[lower:upper + 4]
        else:
            cdr3 = "Couldn't find conserved cysteine"
                        
    
    elif re.findall("FSDG", str(aaseq)) and re.findall('C', str(aaseq)) and locus in ["BCR_K", "K"]:
        indices = [i for i, x in enumerate(aaseq) if x == 'C']
        upper = str(aaseq).find(re.findall("FSDG", str(aaseq))[0])
        lower = False
        for i in indices:
            if i < upper:
                lower = i
        if lower:
            cdr3 = aaseq[lower:upper + 4]
        else:
            cdr3 = "Couldn't find conserved cysteine"

    elif re.findall(".G.G", str(aaseq)):
        cdr3 = "Couldn't find conserved cysteine"
    elif re.findall('C', str(aaseq)):
        cdr3 = "Couldn't find GXG".format(motif_start)
    else:
        cdr3 = "Couldn't find either conserved boundary"
    
    return (cdr3)

def is_cdr3_in_frame(cdr3, locus):
    if "Couldn't" not in cdr3:
        cdr3_in_frame = True
    else:
        cdr3_in_frame = False

    return (cdr3_in_frame)


def collapse_close_sequences(recombinants, locus):
    contig_names = [r.contig_name for r in recombinants]
    filtered_contig_names = [r.contig_name for r in recombinants]
    uncollapsible_contigs = []
    C_genes = dict()
    if len(recombinants) > 1:
        for i in range(len(recombinants) - 1):
            base_name = recombinants[i].contig_name
            base_seq = recombinants[i].dna_seq
            base_V_segment = set(recombinants[i].V_genes)
            base_J_segment = set(recombinants[i].J_genes)
            base_cdr3 = recombinants[i].cdr3_seq
            base_C_gene = None
            if recombinants[i].C_gene:
                base_C_gene = recombinants[i].C_gene
            if base_C_gene is not None:
                if "*" in base_C_gene:
                    base_C_gene = base_C_gene.split("*")[0]
            base_full_length = recombinants[i].full_length

            base_id = recombinants[i].identifier
            base_junc = base_id.split("_")[1]
            base_e_value = float(recombinants[i].hit_table[0][-2])

            for j in range(i + 1, len(recombinants)):
                comp_name = recombinants[j].contig_name
                comp_seq = recombinants[j].dna_seq
                comp_V_segment = set(recombinants[j].V_genes)
                comp_J_segment = set(recombinants[j].J_genes)
                comp_id = recombinants[j].identifier
                comp_junc = comp_id.split("_")[1]
                comp_cdr3 = recombinants[j].cdr3_seq
                comp_C_gene = None
                if recombinants[j].C_gene:
                    comp_C_gene = recombinants[j].C_gene
                if comp_C_gene is not None:
                    if "*" in comp_C_gene:
                        comp_C_gene = comp_C_gene.split("*")[0]
                comp_full_length = recombinants[i].full_length
                comp_e_value = float(recombinants[j].hit_table[0][-2])


                attempt_collapse = False

                if (base_seq in comp_seq) or (comp_seq in base_seq) \
                        and (base_name and comp_name) in filtered_contig_names:
                    attempt_collapse = True 
                elif base_cdr3 is not None and comp_cdr3 is not None:                
                    if base_cdr3 == comp_cdr3 and base_name in filtered_contig_names \
                            and comp_name in filtered_contig_names:
                        if len(base_V_segment.intersection(comp_V_segment)) > 0 \
                                and len(base_J_segment.intersection(comp_J_segment)) > 0:
                            attempt_collapse = True
                    elif (hamming_dist(base_cdr3, comp_cdr3)/len(base_cdr3) <= 0.07 and base_name 
                        in filtered_contig_names and comp_name in filtered_contig_names):
                        if (len(base_V_segment.intersection(comp_V_segment)) > 0 
                            and len(base_J_segment.intersection(comp_J_segment)) > 0):
                            attempt_collapse = True
                        
                elif base_cdr3 is None or comp_cdr3 is None:
                    if len(base_junc) > 0:
                        if (hamming_dist(base_junc, comp_junc)/len(base_junc) <= 0.07 and base_name in filtered_contig_names 
                            and comp_name in filtered_contig_names):
                            if (len(base_V_segment.intersection(comp_V_segment)) > 0 
                                and len(base_J_segment.intersection(comp_J_segment)) > 0):
                                attempt_collapse = True

                if attempt_collapse is False:
                    uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))
                C_gene = None   

                if attempt_collapse:
                    # Assert isotype
                    if locus == "BCR_H":
                        if not base_C_gene == comp_C_gene:
                            if base_C_gene in ["IGHD", "IGHM"] and \
                            comp_C_gene in ["IGHD", "IGHM"]:
                                C_gene = "IGHDM"
                    # find alignment with lowest E value for V match
                    if (base_e_value <= comp_e_value or comp_seq in base_seq) \
                        and comp_name in filtered_contig_names \
                        and base_name in filtered_contig_names:
                        filtered_contig_names.remove(comp_name)
                        if C_gene is not None:
                            C_genes[base_name] = C_gene
                        elif base_C_gene is None and comp_C_gene:
                            C_genes[base_name] = recombinants[j].C_gene
                    else:
                        if comp_name in filtered_contig_names and \
                            base_name in filtered_contig_names:
                            filtered_contig_names.remove(base_name)
                            if C_gene is not None:
                                C_genes[comp_name] = C_gene
                            elif comp_C_gene is None and base_C_gene:
                                C_genes[comp_name] = recombinants[i].C_gene
          
                else:
                    uncollapsible_contigs.append("{}_vs_{}".format(base_name, 
                                                                  comp_name))

    recombinants_to_delete = []

    for r in recombinants:
        if not r.contig_name in filtered_contig_names:
            recombinants_to_delete.append(r)
        else:
            if r.contig_name in C_genes.keys():
                r.C_gene = C_genes[r.contig_name]

    [recombinants.remove(r) for r in recombinants_to_delete]


    return (recombinants)



def load_kallisto_counts(tsv_file):
    counts = defaultdict(lambda: defaultdict(dict))
    with open(tsv_file) as tsvh:
        for line in tsvh:
            if "BRACER" in line:
                line = line.rstrip()
                line = line.split("\t")
                tags = line[0].split("|")
                receptor = tags[1]
                locus = tags[2]
                contig_name = tags[3]
                tpm = float(line[4])
                
                counts[receptor][locus][contig_name] = tpm
    return dict(counts)



def make_cell_network_from_dna(cells, keep_unlinked, shape, dot, neato, 
                   loci, network_colours, cell_contig_clones, 
                   cells_with_clonal_H, no_multiplets, IGH_networks):
    G = nx.MultiGraph()


    # initialise all cells as nodes

    if shape == 'circle':
        for cell in cells:
            if no_multiplets is False or (no_multiplets is True 
                and not cell.has_excess_recombinants):
                G.add_node(cell, shape=shape, 
                        label=cell.html_style_label_for_circles(loci, 
                        network_colours), sep=0.4, fontname="helvetica neue")

                if cell.bgcolor is not None:
                    G.node[cell]['style'] = 'filled'
                    G.node[cell]['fillcolor'] = cell.bgcolor

    else:
        for cell in cells:
            if no_multiplets is False or (no_multiplets is True 
                and not cell.has_excess_recombinants):

                G.add_node(cell, shape=shape, 
                           label=cell.html_style_label_dna(loci, 
                           network_colours),fontname="helvetica neue")
                if cell.bgcolor is not None:
                    G.node[cell]['style'] = 'filled'
                    G.node[cell]['fillcolor'] = cell.bgcolor

    # Create list of cells belonging to a heavy chain clone group
    cell_names = cells_with_clonal_H
   
    # make edges:
    loci_names = loci
    loci = []
    for locus in ["H", "K", "L"]: 
        if locus in loci_names:
            loci.append(locus)

        
    if len(cells) > 1:
        for i in range(len(cell_names)):
            current_cell = cell_names[i]
            comparison_cells = cell_names[i + 1:]

            for locus in loci:
                col = network_colours["BCR"][locus][0]

                for comparison_cell in comparison_cells:
                    shared_identifiers = 0

                    for cell in cells:
                        if current_cell == cell.name:
                            current_cell = cell
                        elif comparison_cell == cell.name:
                            comparison_cell = cell
                    current_contigs = []
                    comparison_contigs = []
                    if current_cell.name in cell_contig_clones[locus].keys() \
                        and comparison_cell.name in cell_contig_clones[locus].keys():
              
                        for contig_name, clone_group in six.iteritems(
                                cell_contig_clones[locus][current_cell.name]):
                            current_contigs.append(contig_name)
                        for contig_name, clone_group in six.iteritems(
                                cell_contig_clones[locus][comparison_cell.name]):
                            comparison_contigs.append(contig_name)

                    if locus == "H":
                        clonal_H = False

                        # Check if cells share a clonal H chain
                        for current_contig in current_contigs:
                            for comparison_contig in comparison_contigs:
                                if cell_contig_clones["H"][comparison_cell.name][comparison_contig] == \
                                    cell_contig_clones["H"][current_cell.name][current_contig]:
                                    clonal_H = True
                                    shared_identifiers += 1 


                        # Check if cells share nonfunctional chains
                        if len(current_cell.recombinants["BCR"][locus]) > 1 \
                            and len(comparison_cell.recombinants["BCR"][locus]) > 1:
                            for current_recombinant in current_cell.recombinants["BCR"][locus]:
                                if not current_recombinant.productive:
                                    current_id_set = current_recombinant.all_poss_identifiers

                                    for comparison_recombinant in \
                                        comparison_cell.recombinants["BCR"][locus]:
                                        if not comparison_recombinant.productive:
                                            comparison_id_set = \
                                            comparison_recombinant.all_poss_identifiers
                                            if len(current_id_set.intersection(comparison_id_set)) > 0:
                                                if shared_identifiers > 0:
                                                    shared_identifiers += 1
                  
                        if shared_identifiers > 0:
                            width = shared_identifiers * 2
                            G.add_edge(current_cell, comparison_cell, locus, penwidth=width, color=col)    

                    else:

                        # Check if cells share clonal H chain
                        if G.number_of_edges(current_cell, comparison_cell) > 0:
 
                            if len(current_contigs) > 0 and len(comparison_contigs) > 0:
                                for current_contig in current_contigs:
                                    for comparison_contig in comparison_contigs:
                                        if cell_contig_clones[locus][comparison_cell.name][comparison_contig] \
                                            == cell_contig_clones[locus][current_cell.name][current_contig]:
                                            shared_identifiers += 1                          
                        
                        if shared_identifiers > 0:
                            width = shared_identifiers * 2
                            G.add_edge(current_cell, comparison_cell, locus, penwidth=width, color=col)

                    
        # Check if cells sharing at least one functional H 
        # chain also share a nonfunctional light chain
        # Non-productive chains are "shared" if their sets of V- and J genes intersect
        for i in range(len(cell_names)):
            current_cell = cell_names[i]
            comparison_cells = cell_names[i + 1:]

            for comparison_cell in comparison_cells:

                for cell in cells:
                    if current_cell == cell.name:
                        current_cell = cell
                    elif comparison_cell == cell.name:
                        comparison_cell = cell

                if not IGH_networks and len(loci) > 0:
                    threshold = 1
                else:
                    threshold = 0
                if G.number_of_edges(current_cell, comparison_cell) > threshold:

                    for locus in loci:
                        if locus is not "H":
                            shared_identifiers = 0
                            col = network_colours["BCR"][locus][0]

                            if len(current_cell.recombinants["BCR"][locus]) > 0 \
                                and len(comparison_cell.recombinants["BCR"][locus]) > 0:
                                for current_recombinant in \
                                    current_cell.recombinants["BCR"][locus]:
                                    if not current_recombinant.productive:
                                        current_id_set = current_recombinant.all_poss_identifiers

                                        for comparison_recombinant in comparison_cell.recombinants["BCR"][locus]:
                                            if not comparison_recombinant.productive:
                                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                                if len(current_id_set.intersection(comparison_id_set)) > 0:
                                                    shared_identifiers += 1
                                if shared_identifiers > 0:
                                    edge_dict = G.get_edge_data(current_cell, comparison_cell)
                                         
                                    if locus in edge_dict.keys():
                                        width = 4
                                        G.add_edge(current_cell, comparison_cell, locus, 
                                                          penwidth=width, color=col)
                                
                                    else:
                                        width = 2 * shared_identifiers
                                        G.add_edge(current_cell, comparison_cell, locus, 
                                                          penwidth=width, color=col, style="dashed")

        # Remove edges between cells that only share clonal heavy chain if --IGH_networks flag is not provided
        #and H is not only locus listed in command line
        if not IGH_networks and len(loci) > 1:
            for i in range(len(cell_names)):
                current_cell = cell_names[i]
                comparison_cells = cell_names[i + 1:]
                for comparison_cell in comparison_cells:
                    for cell in cells:
                        if current_cell == cell.name:
                            current_cell = cell
                        elif comparison_cell == cell.name:
                            comparison_cell = cell
                    if G.number_of_edges(current_cell, comparison_cell) == 1:
                        G.remove_edge(current_cell, comparison_cell)


    deg = G.degree()

    to_remove = [n for n in deg if deg[n] == 0]

    if len(to_remove) < len(G.nodes()):
        if not shape == 'circle':
            G.remove_nodes_from(to_remove)
            drawing_tool = [dot, '-Gsplines=true', '-Goverlap=false', '-Gsep=0.4']

        else:
            drawing_tool = [dot, '-Gsplines=true', '-Goverlap=false']
    else:
        drawing_tool = [neato, '-Gsplines=true', '-Goverlap=false']



    component_counter = 0
    component_groups = list()
    j = 0
    components = nx.connected_components(G)

    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)

        if not members == []:
            component_groups.append(members)

    return (G, drawing_tool, component_groups)




def draw_network_from_cells(cells, output_dir, output_format, dot, neato, 
                           draw_graphs, loci, network_colours, 
                           cell_contig_clones, cells_with_clonal_H, 
                           no_multiplets, IGH_networks):
    cells = list(cells.values())
    
    network, draw_tool, component_groups = make_cell_network_from_dna(cells, 
                            False, "box", dot, neato, loci, network_colours, 
                                  cell_contig_clones, cells_with_clonal_H, 
                                  no_multiplets, IGH_networks)

    network_file = "{}/clonotype_network_with_identifiers.dot".format(output_dir)

    try:
        nx.write_dot(network, network_file)
    except AttributeError:
        import pydotplus
        nx.drawing.nx_pydot.write_dot(network, network_file)
    if draw_graphs:
        command = draw_tool + ['-o', 
                "{output_dir}/clonotype_network_with_identifiers.{output_format}".format(
                output_dir=output_dir, output_format=output_format), "-T", 
                output_format, network_file]
        subprocess.check_call(command)
    
    network, draw_tool, cgx = make_cell_network_from_dna(cells, False, 
                        "circle", dot, neato, loci, network_colours, 
                        cell_contig_clones, cells_with_clonal_H, 
                        no_multiplets, IGH_networks)

    network_file = "{}/clonotype_network_without_identifiers.dot".format(output_dir)
    try:

        nx.write_dot(network, network_file)
    except AttributeError:
        import pydotplus
        nx.drawing.nx_pydot.write_dot(network, network_file)
    if draw_graphs:
        command = draw_tool + ['-o', 
                "{output_dir}/clonotype_network_without_identifiers.{output_format}".format(
                output_dir=output_dir, output_format=output_format), "-T", 
                output_format, network_file]
        subprocess.check_call(command)
    return (component_groups, network)


def get_component_groups_sizes(cells, loci, G):

    # Should provide G from make_cell_network_from_dna?
    cells = list(cells.values())
    components = nx.connected_components(G)
    component_groups = list()

    singlets = []
    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)
            component_groups.append(members)
        else:
            for cell in component:
                singlets.append(cell.name)

    clonotype_size_counter = Counter([len(x) for x in component_groups])
    clonotype_size_counter.update({1: len(singlets)})

    clonotype_sizes = []
    max_size = max(list(clonotype_size_counter.keys()))
    if max_size < 5:
        for x in range(1, max_size + 1):
            clonotype_sizes.append(clonotype_size_counter[x])
        zero_padding = 5 - len(clonotype_sizes)
        clonotype_sizes = clonotype_sizes + [0] * zero_padding
    else:
        for x in range(1, max_size + 1):
            clonotype_sizes.append(clonotype_size_counter[x])

    return (clonotype_sizes)


def check_config_file(filename):
    if not os.path.isfile(filename):
        print()
        print("Couldn't find config file: {}".format(filename))
        print()
        exit(1)


def detect_read_length(fastq1):
    """Detects read length in order to determine if a second round
    of alignment for heavy chain should be performed."""
    sum_length = 0
    count = 0
    read_length = 0

    if fastq1 is not None:
        if os.path.exists(fastq1):
            if ".gz" in fastq1:
                with gzip.open(fastq1, 'r') as f:
                    for line in f:
                        if count > 10000:
                            break
                        if not line.startswith(b'@') and not line.startswith(b'+'):
                            line = line.strip()
                            length = len(line)
                            sum_length += length
                            count += 1
                            
        
            else:
                with open(fastq1, "r") as input:
                    for line in input:
                        if count < 10000:
                            if not line.startswith("@") and not line.startswith("+"):
                                line = line.strip()
                                length = len(line)
                                sum_length += length
                                count += 1
    if count == 0:
        short_reads = True
    else:
        read_length = sum_length / count
        if read_length > 50:
            short_reads = False
        else:
            short_reads = True

    return(short_reads, read_length)


def run_trim_galore(trim_galore, cutadapt, output_dir, cell_name, fastq1,
                    fastq2, should_resume, single_end):
    """Trimming raw reads using Trim Galore/Cutadapt to remove adapter
    sequences and low-quality reads."""

    print("##Trimming raw reads##")
    # Look for existing trimmed reads
    fastq1_base = fastq1.split(".f")[0]
    if not single_end:
        fastq2_base = fastq2.split(".f")[0]
        if "/" in fastq2_base:
            fastq2_base = fastq2_base.split("/")[-1]
    if "/" in fastq1_base:
        fastq1_base = fastq1_base.split("/")[-1]
    
    trimmed_read_path = "{}/trimmed_reads".format(output_dir)
    if ".gz" in (fastq1 or fastq2):
        ending = "fq.gz"
    else:
        ending = "fq"
    if not single_end:
        fastq1_out = "{}/{}_val_1.{}".format(trimmed_read_path, fastq1_base, ending)
        fastq2_out = "{}/{}_val_2.{}".format(trimmed_read_path, fastq2_base, ending)
    else:
        fastq1_out = "{}/R1_trimmed.{}".format(trimmed_read_path, ending)
        fastq2_out = None
    if should_resume:
        if not single_end:
            if os.path.isfile(fastq1_out) and os.path.isfile(fastq2_out):
                print("Resuming with existing trimmed reads")
                return(fastq1_out, fastq2_out)
        else:
            if os.path.isfile(fastq1_out):
                print("Resuming with existing trimmed reads")
                return(fastq1_out, fastq2_out)

    print("Detecting installed version of Cutadapt: ")

    trimmed_read_dir = "{}/trimmed_reads".format(output_dir)

    if not single_end:
        command = [trim_galore, fastq1, fastq2, '--paired', '--suppress_warn',
                  '--path_to_cutadapt', cutadapt, '--no_report_file',
                  '-o', trimmed_read_dir]
    else:
        command = [trim_galore, fastq1, '--suppress_warn', '--path_to_cutadapt',
                  cutadapt, '--no_report_file', '-o', trimmed_read_dir]
    try:
        subprocess.check_call(command)
        print("Trimming completed")
    except subprocess.CalledProcessError:
        print("trim_galore failed")
        
    return(fastq1_out, fastq2_out)



def bowtie2_alignment(bowtie2, ncores, loci, output_dir, cell_name, 
                      synthetic_genome_path, fastq1, fastq2, should_resume, 
                      single_end, bowtie2_build, trimmed_fastq1, trimmed_fastq2):

    print("##Finding recombinant-derived reads##")
    receptor = "BCR"
    initial_locus_names = ["_".join([receptor,x]) for x in loci]
    locus_names = copy.copy(initial_locus_names)
    
    if should_resume:
        for locus in initial_locus_names:
            aligned_read_path = "{}/aligned_reads/{}_{}_".format(output_dir, 
                                                            cell_name, locus)
            fastq1_out = "{}1.fastq".format(aligned_read_path)
            fastq2_out = "{}2.fastq".format(aligned_read_path)
            if os.path.isfile(fastq1_out) and os.path.isfile(fastq2_out):
                print("Resuming with existing {} reads".format(locus))
                locus_names.remove(locus)
    
    if len(locus_names) == 0:
        return
    
    print("Attempting new assembly for {}\n".format(locus_names))
    short_reads, read_length = detect_read_length(fastq1)

    print("Detected average R1 read length: ")
    print(read_length)

    if short_reads:
        print("Short read length detected. BraCeR will run two rounds of alignment for heavy chain\n") 

    for locus in locus_names:
        print("##{}##".format(locus))
        sam_file = "{}/aligned_reads/{}_{}.sam".format(output_dir, cell_name, locus)
        aligned_fasta = "{}/aligned_reads/{}_{}_aligned.fasta".format(
                            output_dir, cell_name, locus)

        if not single_end:
            # Look for trimmed reads
            if (trimmed_fastq1 and trimmed_fastq2) is not None:
                if os.path.isfile(trimmed_fastq1) and os.path.isfile(trimmed_fastq2):
                    fastq1 = trimmed_fastq1
                    fastq2 = trimmed_fastq2
            
            fastq_out_1 = open("{}/aligned_reads/{}_{}_1.fastq".format(
                                    output_dir, cell_name, locus), 'w')
            fastq_out_2 = open("{}/aligned_reads/{}_{}_2.fastq".format(
                                    output_dir, cell_name, locus), 'w')
            
                
            command = [bowtie2, '--no-unal', '-p', ncores, '-k', '1', 
                      '--np', '0', '--rdg', '1,1', '--rfg', '1,1', '-x', 
                      "/".join([synthetic_genome_path, locus]), '-1', 
                      fastq1, '-2', fastq2, '-S', sam_file, '--score-min', 
                      'L,-0.6,-0.6']
            subprocess.check_call(command)

            if locus == "BCR_H" and short_reads == True:
                # Split the sam file for second alignment
                fastq_lines_1, fastq_lines_2, fastq_lines_1_unpaired = split_sam_file_paired(
                                                                        sam_file, fasta=True)
                with open(aligned_fasta, "w") as input:
                    for line in fastq_lines_1:
                        input.write(line)
                    for line in fastq_lines_2:
                        input.write(line)
                    for line in fastq_lines_1_unpaired:
                        input.write(line)
                
                # Create separate sam file for second round of alignment
                sam_file_2 = "{}/aligned_reads/{}_{}_2.sam".format(output_dir,
                                                            cell_name, locus)
                
                # Make Bowtie index of aligned reads
                index_base = os.path.join(
                    output_dir, 'aligned_reads', '{}'.format(locus))
                
                command = [bowtie2_build, '-q', aligned_fasta, index_base]

                if os.path.getsize(aligned_fasta) > 0:
                    try:
                        subprocess.check_call(command)
                    except subprocess.CalledProcessError:
                        print("bowtie2-build failed")
                    try:
                        #Align all reads against aligned reads in local mode
                        command = [bowtie2, '--no-unal', '-p', ncores, '-k', 
                                  '1', '--np', '0', '--rdg', '7,7', '--rfg', 
                                  '7,7', '-x', index_base, '-1', fastq1, '-2', 
                                  fastq2, '-S', sam_file_2, '--local', '--ma', 
                                  '1', '--mp', '20']
                        subprocess.check_call(command)
                        sam_file = sam_file_2
                    except:
                        pass
                

            # Split the sam file for Trinity.
            fastq_lines_1, fastq_lines_2, fastq_lines_1_unpaired = split_sam_file_paired(
                                                                    sam_file, fasta=False)

            for line in fastq_lines_1:
                fastq_out_1.write(line)
            for line in fastq_lines_1_unpaired:
                fastq_out_1.write(line)
            for line in fastq_lines_2:
                fastq_out_2.write(line)
            
            fastq_out_1.close()
            fastq_out_2.close()
        else:
            # Look for trimmed reads
            if trimmed_fastq1 is not None and os.path.isfile(trimmed_fastq1):
                fastq1 = trimmed_fastq1

            fastq_out = open("{}/aligned_reads/{}_{}.fastq".format(output_dir, 
                                                       cell_name, locus), 'w')

            command = [bowtie2, '--no-unal', '-p', ncores, '-k', '1', '--np', 
                      '0', '--rdg', '1,1', '--rfg', '1,1', '-x', 
                      "/".join([synthetic_genome_path, locus]), '-U', fastq1, 
                      '-S', sam_file]

            subprocess.check_call(command)
            
            if locus == "BCR_H":
                # Split the sam file for second alignment
                out_lines = split_sam_file_unpaired(sam_file, fasta=True)

                with open(aligned_fasta, "w") as input:
                    for line in out_lines:
                        input.write(line)

                # Create separate sam file for second round of alignment
                sam_file_2 = "{}/aligned_reads/{}_{}_2.sam".format(output_dir,
                                                        cell_name, locus)

                # Make Bowtie index of aligned reads
                index_base = os.path.join(
                        output_dir, 'aligned_reads', '{}'.format(locus))

                command = [bowtie2_build, '-q', aligned_fasta, index_base]

                if os.path.getsize(aligned_fasta) > 0:
                    try:
                        subprocess.check_call(command)
                    except subprocess.CalledProcessError:
                        print("bowtie2-build failed")
                    try:
                        #Align all reads against aligned reads in local mode
                        command = [bowtie2, '--no-unal', '-p', ncores, '-k',
                                  '1', '--np', '0', '--rdg', '7,7', '--rfg',
                                  '7,7', '-x', index_base, '-U', fastq1, '-S', 
                                  sam_file_2, '--local', '--ma', '1', '--mp', '20']
                        subprocess.check_call(command)
                        sam_file = sam_file_2
                    except:
                        pass

            # Split the sam file for Trinity
            out_lines = split_sam_file_unpaired(sam_file, fasta=False)
            for line in out_lines:
                fastq_out.write(line)
            fastq_out.close()


def split_sam_file_unpaired(sam_file, fasta=False):
    out_lines = []
    with open(sam_file) as sam_in:
        for line in sam_in:
            if not line.startswith("@"):
                line = line.rstrip()
                line = line.split("\t")
                name = line[0]
                seq = line[9]
                qual = line[10]
                flag = int(line[1])
                if not flag == 0:
                    revcomp_flag = "{0:b}".format(flag)[-5]
                else:
                    revcomp_flag = "0"
                if revcomp_flag == "1":
                    seq = str(Seq(seq).reverse_complement())
                    qual = qual[::-1]
                # Add /1 to name to avoid trouble with Trinity
                if not name[len(name)-2:len(name)] == ("/1" or "/2"):
                    name = name + "/1"
                if fasta==False:
                    out_lines.append("@{name}\n{seq}\n+\n{qual}\n".format(
                                            name=name, seq=seq, qual=qual))
                else:
                    out_lines.append(">{}\n{}\n".format(name, seq))

    return(out_lines)


def split_sam_file_paired(sam_file, fasta=False):
    fastq_lines_1 = []
    fastq_lines_2 = []
    fastq_lines_1_unpaired = []
    with open(sam_file) as sam_in:
        for line in sam_in:
            if not line.startswith("@"):
                line = line.rstrip()
                line = line.split("\t")
                name = line[0]
                seq = line[9]
                qual = line[10]
                flag = int(line[1])
                mate_flag = "{0:b}".format(flag)[-7]
                mate_mapped_flag = "{0:b}".format(flag)[-4]
                revcomp_flag = "{0:b}".format(flag)[-5]
                if revcomp_flag == "1":
                    seq = str(Seq(seq).reverse_complement())
                    qual = qual[::-1]
                if mate_mapped_flag == "0":
                    if mate_flag == "1":
                        name_ending = "/1"
                        if fasta==False:
                            fastq_lines_1.append(
                                "@{name}{name_ending}\n{seq}\n+\n{qual}\n".format(
                                name=name, seq=seq, name_ending=name_ending, 
                                qual=qual))
                        else:
                            fastq_lines_1.append(
                                ">{name}{name_ending}\n{seq}\n".format(
                                name=name, seq=seq, name_ending=name_ending))
                    else:
                        name_ending = "/2"
                        if fasta==False:
                            fastq_lines_2.append(
                                "@{name}{name_ending}\n{seq}\n+\n{qual}\n".format(
                                name=name, seq=seq, name_ending=name_ending, 
                                qual=qual))
                        else:
                            fastq_lines_2.append(
                                ">{name}{name_ending}\n{seq}\n".format(
                                name=name, seq=seq, name_ending=name_ending))
                else:
                    name_ending = "/1"
                    if fasta==False:
                        fastq_lines_1_unpaired.append(
                            ">{name}{name_ending}\n{seq}\n+\n{qual}\n".format(
                            name=name, seq=seq, name_ending=name_ending, 
                            qual=qual))
                    else:
                        fastq_lines_1_unpaired.append(
                            "@{name}{name_ending}\n{seq}\n".format(name=name,
                                        seq=seq, name_ending=name_ending))

    return(fastq_lines_1, fastq_lines_2, fastq_lines_1_unpaired)


def assemble_with_trinity(trinity, loci, output_dir, cell_name, ncores, trinity_grid_conf, 
                          JM, version, should_resume, single_end, species, no_normalise):
    print("##Assembling Trinity Contigs##")
    receptor = "BCR"
    if should_resume:
        trinity_report_successful = "{}/Trinity_output/successful_trinity_assemblies.txt".format(
                                                                                    output_dir)
        trinity_report_unsuccessful = "{}/Trinity_output/unsuccessful_trinity_assemblies.txt".format(
                                                                                    output_dir)
        if (os.path.isfile(trinity_report_successful) and os.path.isfile(trinity_report_unsuccessful)) and (
                        os.path.getsize(trinity_report_successful) > 0 or os.path.getsize(
                    trinity_report_unsuccessful) > 0):
            print("Resuming with existing Trinity output")
            successful_files = glob.glob("{}/Trinity_output/*.fasta".format(output_dir))
            return(successful_files)

    base_command = [trinity]
    if trinity_grid_conf:
        base_command = base_command + ['--grid_conf', trinity_grid_conf]

    memory_string = '--max_memory' if (version == '2') else '--JM'
    base_command = base_command + ['--seqType', 'fq', memory_string, JM, 
            '--CPU', ncores, '--full_cleanup']
    if no_normalise:
        base_command = base_command + ['--no_normalize_reads']

    locus_names = ["_".join([receptor,x]) for x in loci]
    
    for locus in locus_names:
        print("##{}##".format(locus))
        trinity_output = "{}/Trinity_output/{}_{}".format(output_dir, cell_name, locus)
        aligned_read_path = "{}/aligned_reads/{}_{}".format(output_dir, cell_name, locus)
        if not single_end:
            file1 = "{}_1.fastq".format(aligned_read_path)
            file2 = "{}_2.fastq".format(aligned_read_path)
            command = base_command + ["--left", file1, "--right", file2, "--output",
                    '{}/Trinity_output/Trinity_{}_{}'.format(output_dir, cell_name, locus)]
        else:
            file = "{}.fastq".format(aligned_read_path)
            command = base_command + ["--single", file, "--output",
                    '{}/Trinity_output/Trinity_{}_{}'.format(output_dir, cell_name, locus)]
        try:
            subprocess.check_call(command)
            shutil.move('{}/Trinity_output/Trinity_{}_{}.Trinity.fasta'.format(
                output_dir, cell_name, locus),
                '{}/Trinity_output/{}_{}.Trinity.fasta'.format(output_dir, cell_name, locus))
        except (subprocess.CalledProcessError, IOError):
            print("Trinity failed for locus")

    # Clean up unsuccessful assemblies
    sleep(10)  # this gives the cluster filesystem time to catch up and stops weird things happening
    successful_files = glob.glob("{}/Trinity_output/*.fasta".format(output_dir))
    unsuccessful_directories = next(os.walk("{}/Trinity_output".format(output_dir)))[1]
    for directory in unsuccessful_directories:
        shutil.rmtree("{}/Trinity_output/{}".format(output_dir, directory))
    successful_file_summary = "{}/Trinity_output/successful_trinity_assemblies.txt".format(output_dir)
    unsuccessful_file_summary = "{}/Trinity_output/unsuccessful_trinity_assemblies.txt".format(output_dir)

    successful_files = bracerlib.io.clean_file_list(successful_files)
    unsuccessful_directories = bracerlib.io.clean_file_list(unsuccessful_directories)

    success_out = open(successful_file_summary, "w")
    fail_out = open(unsuccessful_file_summary, "w")

    successful = defaultdict(list)
    unsuccessful = defaultdict(list)

    successful_ordered_files = set()
    unsuccessful_ordered_files = set()

    for filename in successful_files:
        parsed_name = bracerlib.io.get_filename_and_locus(filename)
        successful[parsed_name[0]].append(parsed_name[1])
        successful_ordered_files.add(parsed_name[0])
    successful_ordered_files = sorted(list(successful_ordered_files))

    for filename in unsuccessful_directories:
        parsed_name = bracerlib.io.get_filename_and_locus(filename)
        unsuccessful[parsed_name[0]].append(parsed_name[1])
        unsuccessful_ordered_files.add(parsed_name[0])
    unsuccessful_ordered_files = sorted(list(unsuccessful_ordered_files))

    successful = bracerlib.io.sort_locus_names(successful)
    unsuccessful = bracerlib.io.sort_locus_names(unsuccessful)

    for file in successful_ordered_files:
        success_out.write("{}\t{}\n".format(file, successful[file]))

    for file in unsuccessful_ordered_files:
        fail_out.write("{}\t{}\n".format(file, unsuccessful[file]))

    success_out.close()
    fail_out.close()

    # Remove pointless .readcount files
    readcount_files = glob.glob("{}/aligned_reads/*.readcount".format(output_dir))
    for f in readcount_files:
        os.remove(f)

    return successful_files



def run_IgBlast(igblast, loci, output_dir, cell_name, index_location, 
                ig_seqtype, species, should_resume, assembled_file):
    """Running IgBlast for reconstructed sequences in a cell using ungapped
    IMGT reference sequences"""

    receptor = "BCR" 

    print("##Running IgBLAST##")

    if assembled_file is None:
        print ("Ig_seqtype:", ig_seqtype)
    
    if 'Hsap' in species:
        igblast_species = 'human'
    elif 'Mmus' in species:
        igblast_species = 'mouse'
    else:
        species_mapper = {
            'Mmus': 'mouse',
            'Hsap': 'human',
            'Rat': 'rat'
        }
        igblast_species = species
        if species in species_mapper.keys():
            igblast_species = species_mapper[species]
        else:
            # Set species to a species recognised by IgBlast
            igblast_species = 'mouse'
    initial_locus_names = ["_".join([receptor,x]) for x in loci]
    locus_names = copy.copy(initial_locus_names)

    if should_resume:
        for locus in initial_locus_names:
            igblast_out = "{}/IgBLAST_output/{}_{}_{}.IgBLASTOut".format(
                                    output_dir, cell_name, receptor, locus)
            if (os.path.isfile(igblast_out) and os.path.getsize(igblast_out) > 0):
                locus_names.remove(locus)
                print("Resuming with existing IgBLAST output for {}".format(locus))
        
        if len(locus_names) == 0:    
            return
    
    print("Performing IgBlast on {}".format(locus_names))

    databases = {}
    for segment in ['V', 'D', 'J']:
        databases[segment] = "{}/{}_{}.fa".format(index_location, receptor, segment)

    # Lines below suppress Igblast warning about not having an auxliary file.
    # Taken from http://stackoverflow.com/questions/11269575/how-to-hide-output-
    #of-subprocess-in-python-2-7
    DEVNULL = open(os.devnull, 'wb')

    num_alignments_V = '20'
    num_alignments_D = '3'
    num_alignments_J = '5'
    
    if assembled_file is not None:
        trinity_fasta = "{}/Trinity_output/{}.fasta".format(output_dir, cell_name)       

    for locus in locus_names:
        if assembled_file is None:
            print("##{}##".format(locus))
            trinity_fasta = "{}/Trinity_output/{}_{}.Trinity.fasta".format(
                                                output_dir, cell_name, locus)
        
        if os.path.isfile(trinity_fasta):
            command = [igblast, '-germline_db_V', databases['V'], 
                      '-germline_db_J', databases['J'], '-germline_db_D', 
                      databases['D'], '-domain_system', 'imgt', '-organism', 
                      igblast_species, '-ig_seqtype', ig_seqtype, 
                      '-show_translation', '-num_alignments_V', 
                      num_alignments_V, '-num_alignments_D', num_alignments_D, 
                      '-num_alignments_J', num_alignments_J, '-outfmt', '7', 
                      '-query', trinity_fasta]

            if assembled_file is None:
                igblast_out = "{}/IgBLAST_output/{}_{}.IgBLASTOut".format(
                                                output_dir, cell_name, locus)
            else:
                igblast_out = "{}/IgBLAST_output/{}.IgBLASTOut".format(
                                                output_dir, cell_name)
                                                                                              
 
            with open(igblast_out, 'w') as out:
                subprocess.check_call(command, stdout=out, stderr=DEVNULL)
            if assembled_file is not None:
                break
    DEVNULL.close()


def run_IgBlast_IMGT_gaps_for_cell(igblast, loci, output_dir, cell_name, 
            ungapped_index_location, gapped_index_location, ig_seqtype, 
            species, assembled_file):
    """Running IgBlast for reconstructed sequences in a cell using IMGT-gapped
    reference sequences in order to determine CDR3 sequences and productivity"""
    receptor = "BCR"
    if 'Hsap' in species:
        igblast_species = 'human'
    elif 'Mmus' in species:
        igblast_species = 'mouse'
    else:
        species_mapper = {
            'Mmus': 'mouse',
            'Hsap': 'human',
            'Rat' : 'rat',
            'rat' : 'rat',
            'Rno' : 'rat'
        }

        igblast_species = species
        if species in species_mapper.keys():
            igblast_species = species_mapper[species]
        else:
            igblast_species = 'mouse'

    initial_locus_names = ["_".join([receptor,x]) for x in loci]
    locus_names = copy.copy(initial_locus_names)


    auxiliary_data = "{}/{}_gl.aux".format(gapped_index_location, species)

    if assembled_file is not None:
        sequence_file = "{}/Trinity_output/{}.fasta".format(output_dir, cell_name)
        output_file = "{}/IgBLAST_output/{}.fmt7".format(output_dir, cell_name)

    for locus in locus_names:
        databases = {}
        databases['V'] = "{}/{}_V.fa".format(gapped_index_location, locus)
        databases['J'] = "{}/{}_J.fa".format(gapped_index_location, locus)
        
        databases['D'] = "{}/BCR_H_D.fa".format(gapped_index_location)

        if assembled_file is None:
            sequence_file = "{}/Trinity_output/{}_{}.Trinity.fasta".format(
                                                output_dir, cell_name, locus)
            output_file = "{}/IgBLAST_output/{}_{}.fmt7".format(
                                                output_dir, cell_name, locus)
            
        if os.path.isfile(sequence_file):
            command = [igblast, '-germline_db_V', databases['V'], '-germline_db_J', 
                         databases['J'], '-germline_db_D', databases['D'],
                         '-domain_system', 'imgt',
                         '-organism', igblast_species, '-ig_seqtype', ig_seqtype,
                         '-outfmt', '7 std qseq sseq btop', '-query', sequence_file]
            
            if os.path.isfile(auxiliary_data):
                command += ['-auxiliary_data', auxiliary_data]
            else:
                print("Warning: No auxiliary_data file found")
            with open(output_file, 'w') as out:
                subprocess.check_call(command, stdout=out)
            if assembled_file is not None:
                break


def run_IgBlast_for_lineage_reconstruction(igblast, locus, output_dir, 
        ungapped_index_location, gapped_index_location, ig_seqtype, species):
    """Runs IgBlast using databases constructed from IMGT-gapped V sequences. 
    Needed for germline reconstruction and lineage reconstruction from 
    clonal sequences"""
    
    if 'Hsap' in species:
        igblast_species = 'human'
    elif 'Mmus' in species:
        igblast_species = 'mouse'
    else:
        species_mapper = {
            'Mmus': 'mouse',
            'Hsap': 'human',
            'Rat' : 'rat'
        }

        igblast_species = species
        if species in species_mapper.keys():
            igblast_species = species_mapper[species]
        else:
            igblast_species = 'mouse'

    databases = {}
    databases['V'] = "{}/BCR_{}_V.fa".format(gapped_index_location, locus)
    databases['J'] = "{}/BCR_{}_J.fa".format(gapped_index_location, locus)
    databases['D'] = "{}/BCR_H_D.fa".format(gapped_index_location)


    auxiliary_data = "{}/{}_gl.aux".format(gapped_index_location, species)
    

    sequence_file = "{}/igblast_input_{}.fa".format(output_dir, locus)
    output_file = "{}/igblast_{}.fmt7".format(output_dir, locus)
    if os.path.isfile(sequence_file):
        
        command = [igblast, '-germline_db_V', databases['V'], '-germline_db_J',
                    databases['J'], '-germline_db_D', databases['D'], 
                    '-domain_system', 'imgt', 
                    '-organism', igblast_species, '-ig_seqtype', ig_seqtype,
                    '-outfmt', '7 std qseq sseq btop', '-query', sequence_file]
        if os.path.isfile(auxiliary_data):
            command += ['-auxiliary_data', auxiliary_data]
        else:
            print("Warning: No auxiliary_data file found")

        with open(output_file, 'w') as out:    
            subprocess.check_call(command, stdout=out)


def run_Blast(blast, loci, output_dir, cell_name, index_location, species,
                should_resume, assembled_file):
    receptor = "BCR"

    print("##Running BLAST##") 
    
    initial_locus_names = ["_".join([receptor,x]) for x in loci]
    locus_names = copy.copy(initial_locus_names)

    print("Performing Blast on {locus_names}".format(locus_names = locus_names))

    databases = {}
    
    
    for segment in ['c', 'C']:
        if segment == 'c':
            seg = 'C_all'
        else:
            seg = segment
        databases[segment] = "{}/{}_{}.fa".format(index_location, receptor, seg)
    
    if (os.path.isfile("{}/{}_C_all.fa".format(index_location, receptor)) and \
        os.path.getsize("{}/{}_C_all.fa".format(index_location, receptor)) > 0):
        database = databases['c']
    else:
        database = databases['C']

  
    if assembled_file is not None:
        bracerlib.io.parse_assembled_file(output_dir, cell_name, assembled_file)
        trinity_fasta = "{}/Trinity_output/{}.fasta".format(output_dir, cell_name)

    for locus in locus_names:
        if assembled_file is None:
            print("##{}##".format(locus))
            trinity_fasta = "{}/Trinity_output/{}_{}.Trinity.fasta".format(
                                            output_dir, cell_name, locus)


        if os.path.isfile(trinity_fasta):
            command = [blast, '-db', database, '-evalue', '0.001', 
                      '-num_alignments', '1', '-outfmt', '5', '-query', 
                      trinity_fasta]
                            
            blast_out = "{}/BLAST_output/{}_{}.xml".format(output_dir, 
                                                    cell_name, locus)

            with open(blast_out, 'w') as out:
                subprocess.check_call(command, stdout=out)



def quantify_with_kallisto(kallisto, cell, output_dir, cell_name, 
                           kallisto_base_transcriptome, fastq1, fastq2, ncores, 
                           should_resume, single_end, fragment_length, 
                           fragment_sd, trimmed_fastq1, trimmed_fastq2,
                           keep_trimmed_reads):

    print("##Running Kallisto##")
    if should_resume:
        if os.path.isfile("{}/expression_quantification/abundance.tsv".format(output_dir)):
            print("Resuming with existing Kallisto output")
            return

    # Look for trimmed reads
    if not single_end:
        if (trimmed_fastq1 and trimmed_fastq2) is not None:
            if os.path.isfile(trimmed_fastq1) and os.path.isfile(trimmed_fastq2):
                fastq1 = trimmed_fastq1
                fastq2 = trimmed_fastq2
    else:
        if trimmed_fastq1 is not None:
            if os.path.isfile(trimmed_fastq1):
                fastq1 = trimmed_fastq1


    print("##Making Kallisto indices##")
    kallisto_dirs = ['kallisto_index']
    for d in kallisto_dirs:
        bracerlib.io.makeOutputDir("{}/expression_quantification/{}".format(
                                                            output_dir, d))
    fasta_filename = "{}/unfiltered_BCR_seqs/{}_BCRseqs.fa".format(output_dir,
                                                                    cell_name)
    fasta_file = open(fasta_filename, 'w')
    fasta_file.write(cell.get_fasta_string())
    fasta_file.close()

    output_transcriptome = "{}/expression_quantification/kallisto_index/{}_transcriptome.fa".format(
                                                            output_dir, cell_name)

    with open(output_transcriptome, 'w') as outfile:
        for fname in [kallisto_base_transcriptome, fasta_filename]:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    idx_file = "{}/expression_quantification/kallisto_index/{}_transcriptome.idx".format(
                                                                output_dir, cell_name)

    index_command = [kallisto, 'index', '--make-unique', '-i', idx_file, output_transcriptome]
    subprocess.check_call(index_command)
    print("##Quantifying with Kallisto##")

    if not single_end:
        if not fragment_length:
            kallisto_command = [kallisto, 'quant', '-i', idx_file, '-t', 
                                ncores, '-o', "{}/expression_quantification".format(
                                output_dir), fastq1, fastq2]
        else:
            kallisto_command = [kallisto, 'quant', '-i', idx_file, '-t', ncores, 
                                '-l', fragment_length, '-o',
                                "{}/expression_quantification".format(
                                output_dir), fastq1, fastq2]

    else: 
        kallisto_command = [kallisto, 'quant', '-i', idx_file, '-t', ncores, 
                            '--single', '-l', fragment_length, '-s', fragment_sd, 
                            '-o', "{}/expression_quantification".format(output_dir), 
                            fastq1] 

    subprocess.check_call(kallisto_command) 

 
    # Delete index file because it's huge and unecessary. Delete transcriptome file 
    shutil.rmtree("{}/expression_quantification/kallisto_index/".format(output_dir))

    # Delete trimmed reads
    if not keep_trimmed_reads:
        for f in [trimmed_fastq1, trimmed_fastq2]: 
            if not f is None:
                if os.path.isfile(f):
                    os.remove(f)


def run_DefineClones(DefineClones, locus, outdir, species, distance): 
    """Runs DefineClones of Change-O"""

    # Set model to Hamming distance if species is not Mmus or Hsap 
    if ("Mmus" or "mouse") in species: 
        model = "mk_rs5nf" 
        dist = "0.2" 
    elif ("Hsap" or "human") in species: 
        model = "hh_s5f" 
        dist = "0.2" 
    else: 
        model = "ham" 
        dist = "0.2" 
 
    if distance:
        dist = distance
      
    changeo_input = "{}/changeo_input_{}.tab".format(outdir, locus) 
    if os.path.isfile(changeo_input) and os.path.getsize(changeo_input) > 0:
        # Check changeo version
        command = [DefineClones, "--version"]
        try:
            changeo_version = subprocess.check_output(command)
        except subprocess.CalledProcessError as err:
            changeo_version = err.output.decode('utf-8')
        print("Running Change-O DefineClones for locus {}".format(locus))
        try:
            changeo_version = changeo_version.decode('utf-8')
            changeo_version = changeo_version.split(".py: ")[1]
            changeo_version = changeo_version.split("-")[0]
            changeo_versions = changeo_version.split(".")
            # Check if changeo-version >= 0.4
            if int(changeo_versions[0])>0 or int(changeo_versions[1]) >= 4:
                command = [DefineClones, '-d', changeo_input, '--mode', 'gene', '--act', 'set',
                '--model', model, '--dist', dist, '--sf', "JUNCTION", '--norm', 'len']
            else:
                command = [DefineClones, "bygroup", '-d', changeo_input, '--mode', 'gene', '--act', 'set',
                '--model', model, '--dist', dist, '--sf', "JUNCTION", '--norm', 'len']
        
        except:
            command = [DefineClones, "bygroup", '-d', changeo_input, '--mode', 'gene', '--act', 'set',  
                    '--model', model, '--dist', dist, '--sf', "JUNCTION", '--norm', 'len'] 
 
        subprocess.check_call(command) 


def run_MakeDb_for_cell(MakeDb, locus, outdir, species, gapped_seq_location, 
                                        ungapped_seq_location, cell_name):
    """Runs MakeDb of Change-O for each cell to identify CDR3 sequences
    during the assembly step"""
    gapped_seqs = {}
    gapped_seqs['V'] = "{}/BCR_{}_V.fa".format(gapped_seq_location, locus)
    gapped_seqs['D'] = "{}/BCR_H_D.fa".format(ungapped_seq_location)
    gapped_seqs['J'] = "{}/BCR_{}_J.fa".format(ungapped_seq_location, locus)

    makedb_input =  "{}/IgBLAST_output/{}_BCR_{}.fmt7".format(outdir, cell_name, locus)
    seq_file = "{}/Trinity_output/{}_BCR_{}.Trinity.fasta".format(outdir, cell_name, locus)    

    if os.path.isfile(makedb_input) and os.path.getsize(makedb_input) > 0:
        if os.path.isfile(seq_file) and os.path.getsize(seq_file) > 0:
            command = [MakeDb, 'igblast', '-i', makedb_input, '-s', seq_file,
                        '-r', gapped_seqs["V"], gapped_seqs["D"],
                        gapped_seqs["J"], '--regions', '--scores']
            subprocess.check_call(command)



def run_MakeDb(MakeDb, locus, outdir, species, gapped_seq_location, 
                                            ungapped_seq_location):
    """Runs MakeDb of Change-O at Summarise level"""

    gapped_seqs = {}
    gapped_seqs['V'] = "{}/BCR_{}_V.fa".format(gapped_seq_location, locus)
    gapped_seqs['D'] = "{}/BCR_H_D.fa".format(ungapped_seq_location)
    gapped_seqs['J'] = "{}/BCR_{}_J.fa".format(ungapped_seq_location, locus)

    makedb_input =  "{}/igblast_{}.fmt7".format(outdir, locus)
    seq_file = "{}/igblast_input_{}.fa".format(outdir, locus)
    

    if os.path.isfile(makedb_input) and os.path.getsize(makedb_input) > 0:
        if os.path.isfile(seq_file) and os.path.getsize(seq_file) > 0:
            command = [MakeDb, 'igblast', '-i', makedb_input, '-s', seq_file, 
                            '-r', gapped_seqs["V"], gapped_seqs["D"],
                            gapped_seqs["J"], '--regions', '--scores']

            subprocess.check_call(command)
                

def run_CreateGermlines(CreateGermlines, locus, outdir, species, 
                        gapped_seq_location, ungapped_seq_location):
    """Runs CreateGermlines of Change-O"""

    gapped_seqs = {}
    gapped_seqs['V'] = "{}/BCR_{}_V.fa".format(gapped_seq_location, locus)
    # Provide D-database even for light chains so IgBlast won't complain
    gapped_seqs['D'] = "{}/BCR_H_D.fa".format(ungapped_seq_location)
    gapped_seqs['J'] = "{}/BCR_{}_J.fa".format(ungapped_seq_location, locus)

    creategermlines_input = "{}/igblast_{}_db-modified.tab".format(outdir, locus)
    if os.path.isfile(creategermlines_input) and os.path.getsize(creategermlines_input) > 0:
        
        command = [CreateGermlines, '-d', creategermlines_input, '-r', gapped_seqs["V"], 
                    gapped_seqs["D"], gapped_seqs["J"], '-g', 'dmask', '--cloned']

        subprocess.check_call(command)
