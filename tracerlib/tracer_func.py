##############################################################################
#                         Functions for use with                             #
# TraCeR - a tool to reconstruct TCR sequences from single-cell RNA-seq data #    
#                                                                            #
# Please see README and LICENCE for details of use and licence conditions.   #
# This software was written by Mike Stubbington (ms31@sanger.ac.uk) from the #
# Teichmann Lab, EMBL-EBI and WTSI (www.teichlab.org). Latest versions are   #
# available for download at www.github.com/teichlab/tracer.                  #
#                                                                            #
#      Copyright (c) 2015, 2016 EMBL - European Bioinformatics Institute     #
#      Copyright (c) 2016 Genome Research Ltd.                               #
#      Author: M.J.T. Stubbington ms31@sanger.ac.uk                          #
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

from tracerlib.core import Cell, Recombinant, Invar_cell
import tracerlib.io

import copy

import pdb

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

def get_C_gene(assembler, output_dir, locus, contig_name):
    if assembler == "basic":
        blast_dir = "Basic_BLAST_output"
    else:
        blast_dir = "BLAST_output"

    locus = locus.split("_")[1]
    blast_summary_file = "{output_dir}/{blast_dir}/blastsummary_{locus}.txt".format(output_dir=output_dir, blast_dir=blast_dir, locus=locus)

    store_details = False
    C_gene = None
    info_line = None
    start_position = None
    with open(blast_summary_file, 'r') as input:
        for line in input:

            if line.startswith("C\t{contig_name}".format(contig_name=contig_name)) or line.startswith("C\treversed|{contig_name}".format(contig_name=contig_name)):
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
        if not region.startswith("CDR3"):
            if start == None:
                start = int(info[1])
        elif region.startswith("CDR3"):
             cdr3_start = int(info[1])
            
    return (start, cdr3_start)    

def find_possible_alignments(sample_dict, locus_names, cell_name, IMGT_seqs, output_dir, species, seq_method,
                             invariant_seqs, loci_for_segments, receptor, loci, max_junc_string_length, assembler):
    alignment_dict = defaultdict(dict)
    recombinants = {}
    for locus in locus_names:
        recombinants[locus] = []
    
    for locus in locus_names:
        data_for_locus = sample_dict[locus]
        if data_for_locus is not None:
            for query_name, query_data in six.iteritems(data_for_locus):
                processed_hit_table = process_hit_table(query_name, query_data, locus)
                if processed_hit_table is not None:
                    (returned_locus, good_hits, rearrangement_summary) = processed_hit_table
                    junction_list = query_data['junction_details']

                    best_V = remove_allele_stars(rearrangement_summary[0].split(",")[0])

                    junc_string = "".join(junction_list)
                    junc_string = remove_NA(junc_string)
                    
                    locus_letter = returned_locus.split("_")[1]
                    
                    if locus_letter in loci_for_segments['D']:
                        has_D = True
                    else:
                        has_D = False
                    
                    if has_D:
                        best_J = remove_allele_stars(rearrangement_summary[2].split(",")[0])
                    else:
                        best_J = remove_allele_stars(rearrangement_summary[1].split(",")[0])


                    identifier = best_V + "_" + junc_string + "_" + best_J

                    ##line attempting to add alignment summary to data for use with PCR comparisons

                    alignment_summary = query_data['alignment_summary']
                    (C_gene, C_info_line) = (None, None)
                    if receptor == "BCR":
                        align_start, cdr3_start = parse_alignment_summary(alignment_summary)
                        (C_gene, C_info_line, C_start) = get_C_gene(assembler, output_dir, locus, query_name)
                    
                    if receptor is not "BCR":
                        all_V_names = [remove_allele_stars(x) for x in rearrangement_summary[0].split(',')]

                    if has_D:
                        all_J_names = [remove_allele_stars(x) for x in rearrangement_summary[2].split(',')]
                    else:
                        all_J_names = [remove_allele_stars(x) for x in rearrangement_summary[1].split(',')]
                    if receptor is not "BCR":
                        all_V_names = [remove_allele_stars(x) for x in rearrangement_summary[0].split(',')]


                    # get original sequence from Trinity/Oases/Basic file - needed for summary of reconstructed lengths.
                    # Only use the VDJ portion found by IgBLAST or trim for BCRs:
                    if assembler == "oases":
                        trinity_file = "{}/Oases_output/{}/MergedAssembly/transcripts.fa".format(output_dir, locus)
                    elif assembler == "basic":
                        if locus == "BCR_H":
                            trinity_file = "{}/Basic_output/basic_H.fa".format(output_dir)
                        else:
                            trinity_file = "{}/Basic_output/basic_L.fa".format(output_dir)
                    else:
                        trinity_file = "{output_dir}/Trinity_output/{cell_name}_{locus}.Trinity.fasta".format(
                        locus=locus, output_dir=output_dir, cell_name=cell_name)
                    with open(trinity_file, 'rU') as tf:
                        for record in SeqIO.parse(tf, 'fasta'):
                            if query_name in record.id:
                    
                                trinity_seq = record
                                if "cell_id=" in query_name:
                                    query_name = query_name.split("=")[1]
                    query_length = query_data["query_length"]
                    if query_length == None:
                        query_length = len(trinity_seq)
                    if 'reversed' in good_hits[0][1]:
                        trinity_seq = trinity_seq.reverse_complement().seq
                    else:
                        trinity_seq = trinity_seq.seq

                    start_coord, end_coord = get_coords(good_hits, receptor)
                    print(end_coord)
                    if receptor == "BCR":
                        if C_gene is not None:
                            #Keep overlapping nt if end of J overlaps with beginning of C 
                            #if (C_start - 1) == end_coord:
                                #end_coord = end_coord + 1 
                            if (C_start - 1) > end_coord:
                                end_coord = C_start - 1
                            
                            print("End coord")
                            print(end_coord)
                    
                    trinity_seq = str(trinity_seq[start_coord:end_coord])
                    (imgt_reconstructed_seq, is_productive, bestVJNames) = get_fasta_line_for_contig_imgt(
                        rearrangement_summary, junction_list, good_hits, returned_locus, IMGT_seqs, cell_name,
                        query_name, species, loci_for_segments)
                    del (is_productive)
                    del (bestVJNames)

                    #Assess if rearrangement is full-length (from start of V gene to start of C gene)
                    full_length = is_rearrangement_full_length(trinity_seq, query_data["hit_table"], query_name, query_length, assembler, output_dir, locus)

                    #query_length = query_data["query_length"]


                    if seq_method == 'imgt':
                        (fasta_line_for_contig, is_productive, bestVJNames) = get_fasta_line_for_contig_imgt(
                            rearrangement_summary, junction_list, good_hits, returned_locus, IMGT_seqs, cell_name,
                            query_name, species, loci_for_segments)
                        cdr3 = get_cdr3(seq, locus)
                        cdr3_seq = None
                        

                    elif seq_method == 'assembly':
                        fasta_line_for_contig = trinity_seq
                        (is_productive, bestVJNames, cdr3, cdr3_seq) = get_fasta_line_for_contig_assembly(trinity_seq, good_hits,
                                                                                          returned_locus, IMGT_seqs,
                                                                                          cell_name, query_name,
                                                                                          loci_for_segments, full_length, receptor, alignment_summary,
                                                                                          rearrangement_summary)

                        print(is_productive, bestVJNames, cdr3)
                    
                    #Identify the most likely V genes if receptor is BCR
                    if receptor == "BCR":
                        if locus in ["H", "BCR_H"]:
                            threshold_percent = 0.05
                        else:
                            threshold_percent = 0.01
                        all_V_names = find_V_genes_based_on_bit_score(trinity_seq, query_data["hit_table"], query_name, threshold_percent)
                          
                    V_genes = all_V_names
                    
                    all_poss_identifiers = set()
                    for V in all_V_names:
                        for J in all_J_names:
                            if receptor != "BCR":
                                i = V + "_" + junc_string + "_" + J
                            else:
                                i = V + "_" + str(len(junc_string)) + "_" + J
                            all_poss_identifiers.add(i)

                    if len(junc_string) < int(max_junc_string_length):
                        rec = Recombinant(contig_name=query_name, locus=returned_locus, identifier=identifier,
                                          all_poss_identifiers=all_poss_identifiers, productive=is_productive[0],
                                          stop_codon=is_productive[1], in_frame=is_productive[2], TPM=0.0,
                                          dna_seq=fasta_line_for_contig, hit_table=good_hits,
                                          summary=rearrangement_summary, junction_details=junction_list,
                                          best_VJ_names=bestVJNames, alignment_summary=alignment_summary,
                                          trinity_seq=trinity_seq, imgt_reconstructed_seq=imgt_reconstructed_seq, 
                                          has_D=has_D, output_dir=output_dir, full_length=full_length, query_length=query_length, 
                                          V_genes=V_genes, cdr3=cdr3, assembler=assembler, C_gene=C_gene, C_info_line=C_info_line, cdr3_seq=cdr3_seq)
                        recombinants[locus].append(rec)
                        print(query_name)

    if recombinants:
        for locus, rs in six.iteritems(recombinants):
            # Adding code to collapse sequences with very low Levenshtein distances caused by confusion between
            # TRAVxD and TRAVx segments with different alignment lengths from IgBlast.
            recombinants[locus] = collapse_close_sequences(rs, locus)

        # cell_name, A_recombinants, B_recombinants, G_recombinants, D_recombinants, is_empty=False, species="Mmus")
        cell = Cell(cell_name, recombinants, species=species, receptor=receptor, loci=loci)
        
    else:
        cell = Cell(cell_name, None, species=species, invariant_seqs=invariant_seqs, receptor=receptor, loci=loci)

    # pdb.set_trace()
    return (cell)



def find_V_genes_based_on_bit_score(seq, hit_table, query_name, threshold_percent):
    """Identifies the most likely V genes if receptor is BCR based on IgBlast bit scores"""
    found_V = False
    V_genes = []
    threshold = None
    for hit in hit_table:
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
        

def parse_rearrangement_summary(rearrangement_summary):
    """Returns a tuple of (stop_codon, in_frame, productive) from IgBlast output for BCRs"""
    chain_type = rearrangement_summary[0][0:4]
    if "IGH" in chain_type:
        i = 0
    else:
        i = 1

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

def get_coords(hit_table, receptor):
    found_V = False
    found_J = False
    for entry in hit_table:
        if entry[0] == 'V':
            if not found_V:
                start = int(entry[8]) - 1
                found_V = True
        if entry[0] == 'J':
            if not found_J:
                end = int(entry[9])
                found_J = True
                bitscore = float(entry[13])
            # Look for other J gene hits with equal bitscore but higher J gene end position
            if receptor == "BCR":
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


def remove_allele_stars(segment):
    p = re.compile(r"(.+)\*\d+")
    m = p.search(segment)
    return (m.group(1))


def process_hit_table(query_name, query_data, locus):
    hit_table = query_data['hit_table']
    rearrangement_summary = query_data['VDJ_rearrangement_summary']

    e_value_cutoff = 5e-3

    found_V = set()
    found_D = set()
    found_J = set()

    good_hits = []

    segment_locus_pattern = re.compile(r"TRAV.+DV.+")
    
    locus_name = locus.split("_")[1]
    
    for entry in hit_table:
        if not entry == "":
          
            entry = entry.split("\t")
            segment = entry[2]
            if segment_locus_pattern.search(segment):
                segment_locus = "AD"
            else:
                segment_locus = segment[2]
            segment_type = segment[3]
            e_value = float(entry[12])
            
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



def get_fasta_line_for_contig_imgt(rearrangement_summary, junction_details, hit_table, locus, IMGT_seqs,
                                   sample_name, query_name, species, loci_for_segments):

    # use first 258 bases of TRBC because they're the same between C1 and C2
    # for TRGC use first 150 bases. Found by aligning the 4 C region transcripts and taking consensus. Ignored start of TCRG-C4-201 because it's only in that one.
    # use first 360 nt of TRBC1 because they're very nearly the same between TRBC1 and TRBCC2

    found_best_V = False
    found_best_D = False
    found_best_J = False
    

    
    V_pattern = re.compile(r".+{potential_loci}V.+".format(potential_loci='[' + "".join(loci_for_segments['V']) + ']'))
    D_pattern = re.compile(r".+{potential_loci}D.+".format(potential_loci='[' + "".join(loci_for_segments['D']) + ']'))
    J_pattern = re.compile(r".+{potential_loci}J.+".format(potential_loci='[' + "".join(loci_for_segments['J']) + ']'))
    

    for hit in hit_table:
        segment = hit[2]
        V_match = V_pattern.search(segment)
        J_match = J_pattern.search(segment)
        if V_match and not found_best_V:
            #V_locus_key = "TR{}V".format(segment[2])
            V_locus_key = "_".join([locus, 'V'])
            best_V_name = segment
            # Remove forward slashes from shared A/D gene names to be the same as in the IMGT files.
            #segment = segment.replace("/", "_")
            best_V_seq = IMGT_seqs[V_locus_key][segment]

            # hit[11] is the end of the V sequence
            best_V_seq = best_V_seq[0:int(hit[11])]
            found_best_V = True
        elif J_match and not found_best_J:
            #J_locus_key = "TR{}J".format(segment[2])
            J_locus_key = "_".join([locus, 'J'])
            best_J_name = segment
            best_J_seq = IMGT_seqs[J_locus_key][segment]
            # hit 10 is the start of the J sequence
            best_J_seq = best_J_seq[int(hit[10]) - 1:]
            found_best_J = True

    junction = []

    parens_pattern = re.compile(r"\([CAGT]+\)")
    
    locus_letter = locus.split("_")[1]
    
    if locus_letter in loci_for_segments['D']:
        # junc_seqs = junction_details[1:3]
        VD_junc = junction_details[1]
        D_region = junction_details[2]
        DJ_junc = junction_details[3]
        if parens_pattern.search(VD_junc):
            VD_junc = re.sub(r'[\(\)]', '', VD_junc)
            length_in_parens = len(VD_junc)
            best_V_seq = best_V_seq[: -length_in_parens]
        if parens_pattern.search(DJ_junc):
            DJ_junc = re.sub(r'[\(\)]', '', DJ_junc)
            length_in_parens = len(DJ_junc)
            best_J_seq = best_J_seq[length_in_parens:]
        junc_seqs = [VD_junc, D_region, DJ_junc]


    else:
        VJ_junc = junction_details[1]
        # junctions in parentheses are represented in the coordinates of the matched segments.
        # Need to trim them then include the NTs in the junction
        if parens_pattern.search(VJ_junc):
            VJ_junc = re.sub(r'[\(\)]', '', VJ_junc)
            length_in_parens = len(VJ_junc)
            best_V_seq = best_V_seq[: -length_in_parens]
            best_J_seq = best_J_seq[length_in_parens:]
        junc_seqs = [VJ_junc]

    for seq in junc_seqs:
        seq = re.sub(r'[\(\)]', '', seq)
        if seq != "N/A":
            junction.append(seq)

    junction = "".join(junction)
    
    constant_seq = list(IMGT_seqs["_".join([locus, 'C'])].values())[0]

    # Editing IMGT V and J sequences to include any alterations from the junction details
    V_end_seq = junction_details[0]
    J_start_seq = junction_details[-1]
    best_V_seq = best_V_seq[:-(len(V_end_seq))]
    
    best_V_seq = best_V_seq + V_end_seq
    best_J_seq = best_J_seq[len(J_start_seq):]
    best_J_seq = J_start_seq + best_J_seq

    full_rearrangement = best_V_seq + junction + best_J_seq + constant_seq
    productive_rearrangement = is_rearrangement_productive(best_V_seq + junction + best_J_seq + constant_seq[0:2])
    # fasta_line = ">chr={}__TCR{}_{}\n{}\n".format(sample_name, locus, query_name, full_rearrangement)

    bestVJ = [best_V_name, best_J_name]

    return (full_rearrangement, productive_rearrangement, bestVJ)


def is_rearrangement_productive(seq):
    # returns a tuple of three true/false values (productive, contains stop, in-frame)
    seq_mod_3 = len(seq) % 3
    if seq_mod_3 == 0:
        in_frame = True
    else:
        in_frame = False

    seq = Seq(seq, IUPAC.unambiguous_dna)
    aa_seq = seq.translate()
    contains_stop = "*" in aa_seq

    if in_frame and not contains_stop:
        productive = True
    else:
        productive = False

    return (productive, contains_stop, in_frame)

 

def is_rearrangement_full_length(seq, hit_table, query_name, query_length, assembler, output_dir, locus):
    found_V = False
    found_J = False
    full_5_prime = False
    full_3_prime = False
    ref_V_start = None
    J_end_pos = None
    V_hit = None
    J_hit = None
    (C_gene, C_info_line, C_position) = get_C_gene(assembler, output_dir, locus, query_name)
    del (C_info_line)
    del (C_position)
   
    for hit in hit_table:
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
        else:
            if full_5_prime == False:
                full_length = False
            else:
                if J_end_pos is not None:
           
                    if ("HJ" in J_hit and J_ref_end_pos > 40) or (("KJ" or "LJ") in J_hit and J_ref_end_pos > 30):
                        if int(query_length)>= (J_end_pos - 1):
                            full_length = True
                    else:
                        full_length = False
                else:
                    full_length = False
    return (full_length)


def get_segment_name(name, pattern):
    match = pattern.search(name)
    number = match.group(1)
    if match.group(3):
        sub_number = match.group(3)
    else:
        sub_number = ""
    return (number)


def get_fasta_line_for_contig_assembly(trinity_seq, hit_table, locus, IMGT_seqs, sample_name, 
                                        query_name, loci_for_segments, full_length, receptor, alignment_summary, rearrangement_summary):
    found_best_V = False
    found_best_D = False
    found_best_J = False

    V_pattern = re.compile(r".+{potential_loci}V.+".format(potential_loci='[' + "".join(loci_for_segments['V']) + ']'))
    D_pattern = re.compile(r".+{potential_loci}D.+".format(potential_loci='[' + "".join(loci_for_segments['D']) + ']'))
    J_pattern = re.compile(r".+{potential_loci}J.+".format(potential_loci='[' + "".join(loci_for_segments['J']) + ']'))

    for hit in hit_table:
        segment = hit[2]
        if V_pattern.search(segment) and not found_best_V:
            V_locus_key = V_locus_key = "_".join([locus, 'V'])
            best_V_name = segment
            # Remove forward slashes from shared A/D gene names to be the same as in the IMGT files.
            segment = segment.replace("/", "_")
            ref_V_seq = IMGT_seqs[V_locus_key][segment]

            # hit[11] is the end of the V sequence
            # best_V_seq = best_V_seq[0:int(hit[11])]
            found_best_V = True
        elif J_pattern.search(segment) and not found_best_J:
            J_locus_key = "_".join([locus, 'J'])
            best_J_name = segment
            ref_J_seq = IMGT_seqs[J_locus_key][segment]
            # hit 10 is the start of the J sequence
            # best_J_seq = best_J_seq[int(hit[10])-1 :]
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
            if receptor == "BCR":
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
    
    if receptor == "TCR":
        cdr3 = get_cdr3(trinity_seq, locus)
        end_padding = (ref_J_length - ref_J_end)
        full_effective_length = start_padding + len(
        trinity_seq) + end_padding + 2  # add two because need first two bases of constant region to put in frame.
        if full_effective_length % 3 == 0:
            in_frame = True
        else:
            in_frame = False
    
        # remove the minimal nucleotides from the trinity sequence to check for stop codons
        start_base_removal_count = (3 - (new_V_start - 1)) % 3
        end_base_removal_count = (1 - end_padding) % 3
    
        seq = trinity_seq[start_base_removal_count:-end_base_removal_count]
        seq = Seq(seq, IUPAC.unambiguous_dna)

        aa_seq = seq.translate()
        contains_stop = "*" in aa_seq

        if in_frame == True and not contains_stop:
            productive = True
        else:
            productive = False

    else:
        (contains_stop, in_frame, productive) = parse_rearrangement_summary(rearrangement_summary)
        (start, stop) = get_coords(hit_table, receptor)
        del (stop)
        if in_frame == False:
            cdr3 = "Couldn't find CDR3"
        else:
            (align_start, cdr3_start) = parse_alignment_summary(alignment_summary)
            cdr3_start = cdr3_start - start
            print("cdr3 start")
            print(cdr3_start)
        
        # remove the minimal nucleotides from the trinity sequence to check for stop codons, working from the position at which the cdr3 region starts (first nt in frame)
        """start_base_removal_count = (cdr3_start - 1) % 3   # There are cdr3_start - 1 nt before beginning of CDR3. These need to be in same reading frame as CDR3.
        print("remove start")
        print(start_base_removal_count)
 
        seq = trinity_seq[start_base_removal_count:]
        print("trimmed start")
        print(seq)
        end_base_removal_count = len(seq) % 3
        print("seq length")
        print(len(seq))
        print("remove end")
        print(end_base_removal_count)
        if end_base_removal_count != 0:
            seq = seq[:-end_base_removal_count]
        print(seq)"""

        # Trim sequence to look for CDR3 in correct region
        if in_frame:
            print("XXXXXXXXXXXX TRIMMING TO DETECT CDR3 XXXXXXXXXXXXXX")
            print(trinity_seq)
            seq = trinity_seq[cdr3_start - 7:]
            print(seq)
            cdr3 = get_cdr3(seq, locus)
            cdr3_length = len(cdr3)
            # Get CDR3 nt sequence excluding conserved Cys and XGXG motif
            cdr3_seq = trinity_seq[cdr3_start-1:cdr3_length*3+cdr3_start-16]
 
        print(cdr3)
        print(cdr3_length)
        print(cdr3_seq)
        #in_frame = is_cdr3_in_frame(cdr3, locus)


    #seq = Seq(seq, IUPAC.unambiguous_dna)
    #aa_seq = seq.translate()

    
    print()

    productive_rearrangement = (productive, contains_stop, in_frame)
    bestVJ = [best_V_name, best_J_name]

    return (productive_rearrangement, bestVJ, cdr3, cdr3_seq)


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
                        

    
    elif re.findall("FSDG", str(aaseq)) and re.findall('C', str(aaseq)):
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
    # pdb.set_trace()
    contig_names = [r.contig_name for r in recombinants]
    filtered_contig_names = [r.contig_name for r in recombinants]
    print("Recs before collapsing close sequences")
    for r in filtered_contig_names:
        print(r)
    
    uncollapsible_contigs = []
    if len(recombinants) > 1:
        for i in range(len(recombinants) - 1):
            base_name = recombinants[i].contig_name
            base_seq = recombinants[i].imgt_reconstructed_seq
            base_V_segment = recombinants[i].best_VJ_names[0]
            base_J_segment = recombinants[i].best_VJ_names[1]

            base_id = recombinants[i].identifier
            base_junc = base_id.split("_")[1]
            base_e_value = float(recombinants[i].hit_table[0][-2])

            for j in range(i + 1, len(recombinants)):
                comp_name = recombinants[j].contig_name
                comp_seq = recombinants[j].imgt_reconstructed_seq
                comp_V_segment = recombinants[j].best_VJ_names[0]
                comp_J_segment = recombinants[j].best_VJ_names[1]
                comp_id = recombinants[j].identifier
                comp_junc = comp_id.split("_")[1]
                comp_e_value = float(recombinants[j].hit_table[0][-2])
                lev_dist = Levenshtein.distance(base_seq, comp_seq)
                # print("{}\t{}\t{}".format(base_id, comp_id, lev_dist))
                if lev_dist < 35 and not base_id == comp_id and base_name in filtered_contig_names \
                        and comp_name in filtered_contig_names:
                    # pdb.set_trace()
                    # define re pattern here to find TRAVx[DN] or TRDVx[DN] depending on locus
                    if locus == "TCRA":
                        duplicate_pattern = re.compile(r"TRAV\d+[DN]")
                        segment_pattern = re.compile(r"TRAV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    elif locus == "TCRD":
                        duplicate_pattern = re.compile(r"DV\d+[DN]")
                        segment_pattern = re.compile(r"DV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    else:
                        uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))
                        attempt_collapse = False
                    if attempt_collapse and (
                        duplicate_pattern.search(base_V_segment) or duplicate_pattern.search(comp_V_segment)):
                        base_segment = get_segment_name(base_V_segment, segment_pattern)
                        comp_segment = get_segment_name(comp_V_segment, segment_pattern)
                        if base_segment == comp_segment:
                            # find alignment with lowest E value for V match
                            if base_e_value <= comp_e_value:
                                filtered_contig_names.remove(comp_name)
                            else:
                                filtered_contig_names.remove(base_name)
                        else:
                            uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))

                    else:
                        uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))

                elif lev_dist < 75 and not base_id == comp_id and base_name in filtered_contig_names \
                        and comp_name in filtered_contig_names:
                    if locus == "TCRA":
                        duplicate_pattern = re.compile(r"TRAV\d+[DN]")
                        segment_pattern = re.compile(r"TRAV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    elif locus == "TCRD":
                        duplicate_pattern = re.compile(r"DV\d+[DN]")
                        segment_pattern = re.compile(r"DV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    else:
                        uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))
                        attempt_collapse = False

                    if attempt_collapse and (
                        duplicate_pattern.search(base_V_segment) or duplicate_pattern.search(comp_V_segment)):
                        base_segment = get_segment_name(base_V_segment, segment_pattern)
                        comp_segment = get_segment_name(comp_V_segment, segment_pattern)
                        if (base_segment == comp_segment) and (base_junc == comp_junc) and (
                            base_J_segment == comp_J_segment):
                            # find alignment with lowest E value for V match
                            if base_e_value <= comp_e_value:
                                filtered_contig_names.remove(comp_name)
                            else:
                                filtered_contig_names.remove(base_name)
                        else:
                            uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))

                elif base_id == comp_id and base_name in filtered_contig_names and comp_name in filtered_contig_names:
                    if base_e_value <= comp_e_value:
                        filtered_contig_names.remove(comp_name)
                    else:
                        filtered_contig_names.remove(base_name)

    recombinants_to_delete = []

    for r in recombinants:
        if not r.contig_name in filtered_contig_names:
            recombinants_to_delete.append(r)

    [recombinants.remove(r) for r in recombinants_to_delete]
    print("Recs after collapsing close sequences")
    for r in recombinants:
        print(r.contig_name)
       

    return (recombinants)


def load_kallisto_counts(tsv_file):
    counts = defaultdict(lambda: defaultdict(dict))
    with open(tsv_file) as tsvh:
        for line in tsvh:
            if "TRACER" in line:
                line = line.rstrip()
                line = line.split("\t")
                tags = line[0].split("|")
                receptor = tags[1]
                locus = tags[2]
                contig_name = tags[3]
                tpm = float(line[4])
                
                counts[receptor][locus][contig_name] = tpm
    return dict(counts)



"""def make_cell_network_from_dna_B_cells(cells, keep_unlinked, shape, dot, neato, receptor, loci,
                               network_colours):
    G = nx.MultiGraph()

    #H_clonal_groups = define_potential_H_clonal_groups(cells, receptor)
    # initialise all cells as nodes
    
    if shape == 'circle':
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_for_circles(receptor, loci, network_colours),
                        sep=0.4, fontname="helvetica neue")
            #print(cell.bgcolor)
            if cell.bgcolor is not None:
                G.node[cell]['style'] = 'filled'
            
                G.node[cell]['fillcolor'] = cell.bgcolor
            #print(cell.bgcolor)
            #print(cell.name, cell.isotype, cell.bgcolor)

    else:
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_dna(receptor, loci, network_colours),
                        fontname="helvetica neue")
            if cell.bgcolor is not None:
                G.node[cell]['style'] = 'filled'

                G.node[cell]['fillcolor'] = cell.bgcolor
    # make edges:
    for i in range(len(cells)):
        current_cell = cells[i]
        
        comparison_cells = cells[i + 1:]
        for locus in loci:
            col = network_colours[receptor][locus][0]

            # current_identifiers = current_cell.getMainRecombinantIdentifiersForLocus(locus)
            for comparison_cell in comparison_cells:
                shared_identifiers = 0
                if current_cell.recombinants[receptor][locus] is not None:
                    for current_recombinant in current_cell.recombinants[receptor][locus]:
                        current_id_set = current_recombinant.all_poss_identifiers
                        #print("Current ID set")
                        #print(current_id_set)
                        if comparison_cell.recombinants[receptor][locus] is not None:
                            for comparison_recombinant in comparison_cell.recombinants[receptor][locus]:
                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                #print("Comp ID set")
                                #print(comparison_id_set)
                                if len(current_id_set.intersection(comparison_id_set)) > 0:
                                    shared_identifiers += 1
                                    
                                    
                
                if shared_identifiers > 0:
                    width = shared_identifiers * 2
                    
                    if locus == "H" or G.has_edge(current_cell, comparison_cell):
                    #print("Shared identifiers")
                    #print(shared_identifiers)
                      
                    
                        G.add_edge(current_cell, comparison_cell, locus, penwidth=width, color=col,
                               weight=shared_identifiers)
                
           

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
                

        component_groups.append(members)

    return (G, drawing_tool, component_groups)"""


def make_cell_network_from_dna_B_cells(cells, keep_unlinked, shape, dot, neato, receptor, loci,
                               network_colours, n_edit_distance):
    G = nx.MultiGraph()

    #H_clonal_groups = define_potential_H_clonal_groups(cells, receptor)
    # initialise all cells as nodes

    if shape == 'circle':
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_for_circles(receptor, loci, network_colours),
                        sep=0.4, fontname="helvetica neue")
            
            if cell.bgcolor is not None:
                G.node[cell]['style'] = 'filled'

                G.node[cell]['fillcolor'] = cell.bgcolor
            print(cell.name, cell.isotype, cell.bgcolor)

    else:
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_dna(receptor, loci, network_colours),
                        fontname="helvetica neue")
            if cell.bgcolor is not None:
                G.node[cell]['style'] = 'filled'

                G.node[cell]['fillcolor'] = cell.bgcolor
    # make edges:

    for clone, clone_data in six.iteritems(n_edit_distance):
        for locus, locus_data in six.iteritems(clone_data):
            col = network_colours[receptor][locus][0]
            shared_identifiers = 0
            for (cell1, cell2), distance in six.iteritems(locus_data):
            
                for cell in cells:
                    if cell1 == cell.name:
                        cell1 = cell
                    elif cell2 == cell.name:
                        cell2 = cell
                if (cell1, cell2, locus) in G.edges():
                    width = 4
                    G.add_edge(cell1, cell2, locus, penwidth=width, color=col)
                else:
                    shared_identifiers = 1

                if len(cell1.recombinants[receptor][locus]) > 1 and len(cell2.recombinants[receptor][locus]) > 1:
                    for current_recombinant in cell1.recombinants[receptor][locus]:
                        if not current_recombinant.productive:
                            current_id_set = current_recombinant.all_poss_identifiers
                         
                            for comparison_recombinant in cell2.recombinants[receptor][locus]:
                                if not comparison_recombinant.productive:
                                    comparison_id_set = comparison_recombinant.all_poss_identifiers
                                    if len(current_id_set.intersection(comparison_id_set)) > 0:
                                        shared_identifiers += 1
                if not (cell1, cell2, locus) in G.edges():
                    width = shared_identifiers * 2
                    distance = float(distance)
                    G.add_edge(cell1, cell2, locus, penwidth=width, color=col,
                               dist=distance)
            

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


        component_groups.append(members)

    return (G, drawing_tool, component_groups)

def make_dict_cells_distance_zero(cells, s_normalised_edit_distances):

    zero_distance = dict()
    zero_distance_names = dict()
    for clone, data in six.iteritems(s_normalised_edit_distances):
        zero_distance[clone] = dict()
        zero_distance_names[clone] = dict()
        cell_list = []
        cell_list_names = []
        counter = 0
        
        for (cell1, cell2), distance in  six.iteritems(s_normalised_edit_distances[clone]):
            
            if distance == 0:
                for cell in cells:
                    if cell1 == cell.name:
                        cella = cell 
                    elif cell2 == cell.name:
                        cellb = cell
                if counter == 0:
                    zero_distance[clone][counter] = [cell1, cell2]
                else:
                    for i in range(counter):
                        if (cell1 or cell2) in zero_distance[clone][i]:
                            if not cell1 in zero_distance[clone][i]:
                                zero_distance[clone][i].append(cell1)
                            elif not cell2 in zero_distance[clone][i]:
                                zero_distance[clone][i].append(cell2)
                cells_append = (cella, cellb) 
                cell_list.append(cells_append)
                names = (cell1, cell2)
                cell_list_names.append(names)
        zero_distance[clone] = cell_list
        zero_distance_names[clone] = cell_list_names
                    
    print(zero_distance_names)
    return(zero_distance)


def create_nodes_of_zero_distance_clones(cells, receptor, loci, network_colours, s_normalised_edit_distances):
    H = nx.Graph()
    # Initialise all zero distance cells as nodes
    zero_distance = make_dict_cells_distance_zero(cells, s_normalised_edit_distances)
    if shape == 'circle':
        for cell in cells:
            
            G.add_node(cell, shape=shape, label=cell.html_style_label_for_circles(receptor, loci, network_colours),
                        sep=0.4, fontname="helvetica neue")
            
            if cell.bgcolor is not None:
                G.node[cell]['style'] = 'filled'

                G.node[cell]['fillcolor'] = cell.bgcolor


    else:
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_dna(receptor, loci, network_colours),
                        fontname="helvetica neue")
            if cell.bgcolor is not None:
                G.node[cell]['style'] = 'filled'

                G.node[cell]['fillcolor'] = cell.bgcolor
    zero_distance = make_dict_cells_distance_zero(cells, s_normalised_edit_distances)

    for clone, data in six.iteritems(zero_distance):
        edges = []

        for (cell1, cell2) in data:

            edge = (cell1, cell2)
            edges.append(edge)
        G.add_edges_from(edges)
        #G.add_edge(current_cell, comparison_cell, weight="100")
        #print(edges)

"""def make_cell_network_from_dna_sum_normalised(cells, keep_unlinked, shape, dot, neato, receptor, loci,
                               network_colours, s_normalised_edit_distances):
    G = nx.Graph()


    # initialise all cells as nodes

    if shape == 'circle':
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_for_circles(receptor, loci, network_colours),
                        sep=0.4, fontname="helvetica neue")
            #print(cell.bgcolor)
            if cell.bgcolor is not None:
                G.node[cell]['style'] = 'filled'

                G.node[cell]['fillcolor'] = cell.bgcolor
            #print(cell.bgcolor)
            #print(cell.name, cell.isotype, cell.bgcolor)

    else:
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_dna(receptor, loci, network_colours),
                        fontname="helvetica neue")
            if cell.bgcolor is not None:
                G.node[cell]['style'] = 'filled'

                G.node[cell]['fillcolor'] = cell.bgcolor

    zero_distance = make_dict_cells_distance_zero(cells, s_normalised_edit_distances)
    print("Zero distance: ")
    print(zero_distance)

    # make edges between cells with distance 0:

    #G.add_edge(current_cell, comparison_cell, distance, color="#000000",
                        #weight=distance)

    H = nx.Graph()
    for clone, data in six.iteritems(zero_distance):
        edges = []
        #print("data")
        #print(data)
        for groups, cells_in_group in six.iteritems(zero_distance[clone]):
             H = nx.Graph(G.subgraph(cells_in_group))
             H.nodes()
             G.add_node(H)   
            #edge = (cell1, cell2)
            #edges.append(edge)
        #G.add_edges_from(edges)
        #G.add_edge(current_cell, comparison_cell, weight="100")
        #print(edges)



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
    #H = nx.Graph()
    #for group in component_groups:
        #subgraph = G.subgraph(group)
        #H.add_node(subgraph)

    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)


        component_groups.append(members)



    return (G, drawing_tool, component_groups)"""




"""def make_cell_network_from_dna_sum_normalised(cells, keep_unlinked, shape, dot, neato, receptor, loci,
                               network_colours, s_normalised_edit_distances):
    G = nx.Graph()
    

    # initialise all cells as nodes

    if shape == 'circle':
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_for_circles(receptor, loci, network_colours),
                        sep=0.4, fontname="helvetica neue")
            #print(cell.bgcolor)
            if cell.bgcolor is not None:
                G.node[cell]['style'] = 'filled'

                G.node[cell]['fillcolor'] = cell.bgcolor
            #print(cell.bgcolor)
            #print(cell.name, cell.isotype, cell.bgcolor)

    else:
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_dna(receptor, loci, network_colours),
                        fontname="helvetica neue")
            if cell.bgcolor is not None:
                G.node[cell]['style'] = 'filled'

                G.node[cell]['fillcolor'] = cell.bgcolor

    zero_distance = make_dict_cells_distance_zero(cells, s_normalised_edit_distances)
    #print(zero_distance)

    # make edges between cells with distance 0:

    #G.add_edge(current_cell, comparison_cell, distance, color="#000000",
                        #weight=distance)
     
    for clone, data in six.iteritems(zero_distance): 
        edges = []
        for (cell1, cell2) in data:
             
            edge = (cell1, cell2, 3)
            edges.append(edge)
        G.add_weighted_edges_from(edges)            
        #G.add_edge(current_cell, comparison_cell, weight="100")
        #print(edges)
        


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
    #H = nx.Graph()
    #for group in component_groups:
        #subgraph = G.subgraph(group)
        #H.add_node(subgraph)

    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)


        component_groups.append(members)


  
    return (G, drawing_tool, component_groups)"""





def make_cell_network_from_dna(cells, keep_unlinked, shape, dot, neato, receptor, loci, 
                               network_colours):
    G = nx.MultiGraph()
    # initialise all cells as nodes

    if shape == 'circle':
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_for_circles(receptor, loci, network_colours), 
                        sep=0.4, fontname="helvetica neue")
    else:
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_dna(receptor, loci, network_colours), 
                        fontname="helvetica neue")
    # make edges:
    for i in range(len(cells)):
        current_cell = cells[i]
        comparison_cells = cells[i + 1:]

        for locus in loci:
            col = network_colours[receptor][locus][0]

            # current_identifiers = current_cell.getMainRecombinantIdentifiersForLocus(locus)
            for comparison_cell in comparison_cells:
                shared_identifiers = 0
                if current_cell.recombinants[receptor][locus] is not None:
                    for current_recombinant in current_cell.recombinants[receptor][locus]:
                        current_id_set = current_recombinant.all_poss_identifiers
                        if comparison_cell.recombinants[receptor][locus] is not None:
                            for comparison_recombinant in comparison_cell.recombinants[receptor][locus]:
                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                if len(current_id_set.intersection(comparison_id_set)) > 0:
                                    shared_identifiers += 1

                if shared_identifiers > 0:
                    width = shared_identifiers * 2
                    G.add_edge(current_cell, comparison_cell, locus, penwidth=width, color=col,
                               weight=shared_identifiers)

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

    bgcolors = ['#8dd3c720', '#ffffb320', '#bebada20', '#fb807220', '#80b1d320', '#fdb46220', '#b3de6920', '#fccde520',
                '#d9d9d920', '#bc80bd20', '#ccebc520', '#ffed6f20']
    component_counter = 0
    component_groups = list()
    j = 0
    components = nx.connected_components(G)

    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)
                G.node[cell]['style'] = 'filled'
                G.node[cell]['fillcolor'] = bgcolors[j]
                cell.bgcolor = bgcolors[j]

            if j < 11:
                j += 1
            else:
                component_counter += 1
                j = 0

        component_groups.append(members)

    return (G, drawing_tool, component_groups)


def draw_network_from_cells(cells, output_dir, output_format, dot, neato, draw_graphs, receptor, loci, network_colours, n_edit_distance):
    cells = list(cells.values())
    if not receptor == "BCR":
        network, draw_tool, component_groups = make_cell_network_from_dna(cells, False, "box", dot,
                                                                      neato, receptor, loci, network_colours)
    else:
        network, draw_tool, component_groups = make_cell_network_from_dna_B_cells(cells, False, "box", dot,
                                                                      neato, receptor, loci, network_colours, n_edit_distance)
    network_file = "{}/clonotype_network_with_identifiers.dot".format(output_dir)
    try:
        nx.write_dot(network, network_file)
    except AttributeError:
        import pydotplus
        nx.drawing.nx_pydot.write_dot(network, network_file)
    if draw_graphs:
        command = draw_tool + ['-o', "{output_dir}/clonotype_network_with_identifiers.{output_format}".format(
            output_dir=output_dir, output_format=output_format), "-T", output_format, network_file]
        subprocess.check_call(command)
    if not receptor == "BCR":
        network, draw_tool, cgx = make_cell_network_from_dna(cells, False, "circle", dot, 
                                                         neato, receptor, loci, network_colours)
    else:
        network, draw_tool, cgx = make_cell_network_from_dna_B_cells(cells, False, "circle", dot,
                                                         neato, receptor, loci, network_colours, n_edit_distance)

    network_file = "{}/clonotype_network_without_identifiers.dot".format(output_dir)
    try:

        nx.write_dot(network, network_file)
    except AttributeError:
        import pydotplus
        nx.drawing.nx_pydot.write_dot(network, network_file)
    if draw_graphs:
        command = draw_tool + ['-o', "{output_dir}/clonotype_network_without_identifiers.{output_format}".format(
            output_dir=output_dir, output_format=output_format), "-T", output_format, network_file]
        subprocess.check_call(command)
    return (component_groups)


def get_component_groups_sizes(cells, receptor, loci):
    cells = list(cells.values())
    G = nx.MultiGraph()
    # initialise all cells as nodes
    for cell in cells:
        G.add_node(cell)
    # make edges:
    for i in range(len(cells)):
        current_cell = cells[i]
        comparison_cells = cells[i + 1:]
        clonality = False
        for locus in loci:
            
            # current_identifiers = current_cell.getMainRecombinantIdentifiersForLocus(locus)
            for comparison_cell in comparison_cells:
                shared_identifiers = 0
                if current_cell.recombinants[receptor][locus] is not None:
                    for current_recombinant in current_cell.recombinants[receptor][locus]:
                        current_id_set = current_recombinant.all_poss_identifiers
                        if comparison_cell.recombinants[receptor][locus] is not None:
                            for comparison_recombinant in comparison_cell.recombinants[receptor][locus]:
                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                if len(current_id_set.intersection(comparison_id_set)) > 0:
                                    shared_identifiers += 1

                # comparison_identifiers = comparison_cell.getAllRecombinantIdentifiersForLocus(locus)
                # common_identifiers = current_identifiers.intersection(comparison_identifiers)
                
                if shared_identifiers > 0:
                    width = shared_identifiers * 2
                    if receptor != "BCR":
                        G.add_edge(current_cell, comparison_cell, locus, penwidth=width, weight=shared_identifiers)
                    if receptor == "BCR":
                        if G.has_edge(current_cell, comparison_cell):
                            clonality = True
                        if locus == "H" or G.has_edge(current_cell, comparison_cell):
                    #print("Shared identifiers")
                    #print(shared_identifiers)


                            G.add_edge(current_cell, comparison_cell, locus, penwidth=width, weight=shared_identifiers)
                        


    deg = G.degree()

    to_remove = [n for n in deg if deg[n] == 0]

    # if len(to_remove) < len(G.nodes()):
    #    G.remove_nodes_from(to_remove)

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


def bowtie2_alignment(bowtie2, ncores, receptor, loci, output_dir, cell_name, synthetic_genome_path, fastq1,
                      fastq2, should_resume, single_end):
    print("##Finding recombinant-derived reads##")
    
    initial_locus_names = ["_".join([receptor,x]) for x in loci]
    locus_names = copy.copy(initial_locus_names)
    
    if should_resume:
        for locus in initial_locus_names:
            aligned_read_path = "{}/aligned_reads/{}_{}_".format(output_dir, cell_name, locus)
            fastq1_out = "{}1.fastq".format(aligned_read_path)
            fastq2_out = "{}2.fastq".format(aligned_read_path)
            if os.path.isfile(fastq1_out) and os.path.isfile(fastq2_out):
                print("Resuming with existing {locus} reads".format(locus=locus))
                locus_names.remove(locus)
    
    
    
    if len(locus_names) == 0:
        return
    
    print("Attempting new assembly for {locus_names}\n".format(locus_names=locus_names))
    
    for locus in locus_names:
        print("##{}##".format(locus))
        sam_file = "{}/aligned_reads/{}_{}.sam".format(output_dir, cell_name, locus)
        if not single_end:
            fastq_out_1 = open("{}/aligned_reads/{}_{}_1.fastq".format(output_dir, cell_name, locus), 'w')
            fastq_lines_1 = []
            fastq_out_2 = open("{}/aligned_reads/{}_{}_2.fastq".format(output_dir, cell_name, locus), 'w')
            fastq_lines_2 = []

            command = [bowtie2, '--no-unal', '-p', ncores, '-k', '1', '--np', '0', '--rdg', '1,1', '--rfg', '1,1',
                       '-x', "/".join([synthetic_genome_path, locus]), '-1', fastq1, '-2', fastq2, '-S', sam_file]

            subprocess.check_call(command)

            # now to split the sam file for Trinity.

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
                                fastq_lines_1.append(
                                    "@{name}{name_ending}\n{seq}\n+\n{qual}\n".format(name=name, seq=seq,
                                                                                      name_ending=name_ending,
                                                                                      qual=qual))
                            else:
                                name_ending = "/2"
                                fastq_lines_2.append(
                                    "@{name}{name_ending}\n{seq}\n+\n{qual}\n".format(name=name, seq=seq,
                                                                                      name_ending=name_ending,
                                                                                      qual=qual))

            for line in fastq_lines_1:
                fastq_out_1.write(line)
            for line in fastq_lines_2:
                fastq_out_2.write(line)

            fastq_out_1.close()
            fastq_out_2.close()
        else:
            fastq_out = open("{}/aligned_reads/{}_{}.fastq".format(output_dir, cell_name, locus), 'w')
            command = [bowtie2, '--no-unal', '-p', ncores, '-k', '1', '--np', '0', '--rdg', '1,1', '--rfg', '1,1',
                       '-x', "/".join([synthetic_genome_path, locus]), '-U', fastq1, '-S', sam_file]

            subprocess.check_call(command)

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
                        fastq_out.write("@{name}\n{seq}\n+\n{qual}\n".format(name=name, seq=seq, qual=qual))
                fastq_out.close()


def assemble_with_trinity(trinity, receptor, loci, output_dir, cell_name, ncores, trinity_grid_conf, JM,
                          version, should_resume, single_end, species):
    print("##Assembling Trinity Contigs##")

    if should_resume:
        trinity_report_successful = "{}/Trinity_output/successful_trinity_assemblies.txt".format(output_dir)
        trinity_report_unsuccessful = "{}/Trinity_output/unsuccessful_trinity_assemblies.txt".format(output_dir)
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
    base_command = base_command + ['--seqType', 'fq', memory_string, JM, '--CPU', ncores, '--full_cleanup']
    
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
            shutil.move('{}/Trinity_output/Trinity_{}_{}.Trinity.fasta'.format(output_dir, cell_name, locus),
                        '{}/Trinity_output/{}_{}.Trinity.fasta'.format(output_dir, cell_name, locus))
        except (subprocess.CalledProcessError, IOError):
            print("Trinity failed for locus")

    # clean up unsuccessful assemblies
    sleep(10)  # this gives the cluster filesystem time to catch up and stops weird things happening
    successful_files = glob.glob("{}/Trinity_output/*.fasta".format(output_dir))
    unsuccessful_directories = next(os.walk("{}/Trinity_output".format(output_dir)))[1]
    for directory in unsuccessful_directories:
        shutil.rmtree("{}/Trinity_output/{}".format(output_dir, directory))
    successful_file_summary = "{}/Trinity_output/successful_trinity_assemblies.txt".format(output_dir)
    unsuccessful_file_summary = "{}/Trinity_output/unsuccessful_trinity_assemblies.txt".format(output_dir)

    successful_files = tracerlib.io.clean_file_list(successful_files)
    unsuccessful_directories = tracerlib.io.clean_file_list(unsuccessful_directories)

    success_out = open(successful_file_summary, "w")
    fail_out = open(unsuccessful_file_summary, "w")

    successful = defaultdict(list)
    unsuccessful = defaultdict(list)

    successful_ordered_files = set()
    unsuccessful_ordered_files = set()

    for filename in successful_files:
        # success_out.write("{}\n".format(filename))
        parsed_name = tracerlib.io.get_filename_and_locus(filename)
        successful[parsed_name[0]].append(parsed_name[1])
        successful_ordered_files.add(parsed_name[0])
    successful_ordered_files = sorted(list(successful_ordered_files))

    for filename in unsuccessful_directories:
        # fail_out.write("{}\n".format(filename))
        parsed_name = tracerlib.io.get_filename_and_locus(filename)
        unsuccessful[parsed_name[0]].append(parsed_name[1])
        unsuccessful_ordered_files.add(parsed_name[0])
    unsuccessful_ordered_files = sorted(list(unsuccessful_ordered_files))

    successful = tracerlib.io.sort_locus_names(successful)
    unsuccessful = tracerlib.io.sort_locus_names(unsuccessful)

    for file in successful_ordered_files:
        success_out.write("{}\t{}\n".format(file, successful[file]))

    for file in unsuccessful_ordered_files:
        fail_out.write("{}\t{}\n".format(file, unsuccessful[file]))

    success_out.close()
    fail_out.close()

    # remove pointless .readcount files
    readcount_files = glob.glob("{}/aligned_reads/*.readcount".format(output_dir))
    for f in readcount_files:
        os.remove(f)

    # if len(unsuccessful_directories) == 2:

    return successful_files


def run_velvet_h(velveth, velvetg, oases, receptor, loci, output_dir, cell_name, ncores, should_resume, single_end, species):
    locus_names = ["_".join([receptor,x]) for x in loci]    
    
    # Create single-k assemblies
    for locus in locus_names:
        print("##{}##".format(locus))
        #Set hash lengths and output paths
        oases_output_path = "{}/Oases_output/{}".format(output_dir, locus)
        tracerlib.io.makeOutputDir(oases_output_path)
        hash_lengths = ["15", "17", "19", "21", "23", "25", "27", "29", "31"]
        hash_length_paths = []
        for hash_length in hash_lengths:
            path = "{}/directory{}".format(oases_output_path, hash_length)
            hash_length_paths.append(path)
            #oases_output = "{}/Oases_output/{}_{}".format(output_dir, cell_name, locus)
            aligned_read_path = "{}/aligned_reads/{}_{}".format(output_dir, cell_name, locus)
            velveth_commands = []
            velvetg_commands = []
            oases_commands = []
        if not single_end:
            file1 = "{}_1.fastq".format(aligned_read_path)
            file2 = "{}_2.fastq".format(aligned_read_path)
            # Create velveth, velvetg and oases commands
            print("##### Creating velveth, velvetg and oases commands#####")
            for i in range(len(hash_lengths)):
                path = hash_length_paths[i]  
                hash_length = hash_lengths[i]
                velveth_command = [velveth, path, hash_length, '-fastq', '-shortPaired', '-separate', file1, file2]
                velveth_commands.append(velveth_command) 
                velvetg_command = [velvetg, path, '-read_trkg', 'yes', '-min_contig_lgth', '400']
                #'-exp_aov', 'auto', '-ins_length', '500'
                velvetg_commands.append(velvetg_command)
                oases_command = [oases, path, '-min_trans_lgth', '400'] 
                #'-cov_cutoff', 'auto'
                oases_commands.append(oases_command)
                
                #else:
                    #file = "{}.fastq".format(aligned_read_path)
                    #command = base_command + ["--single", file, "--output",
                                 #'{}/Oases_output/Oases_{}_{}'.format(output_dir, cell_name, locus)]

        for command_list in [velveth_commands, velvetg_commands, oases_commands]:
            for command in command_list:
                if command_list == velveth_commands:
                    program = "velveth"
                elif command_list == velvetg_commands:
                    program = "velvetg"
                elif command_list == oases_commands:
                    program = "oases"
                try:    
                    subprocess.check_call(command)
                    """shutil.move('{}/Oases_output/pairedEnd*'.format(output_dir, cell_name, locus),
                        '{}/Oases_output/{}_{}.Oases.fasta'.format(output_dir, cell_name, locus))"""
                    print("Single k-mer assembly successful for locus when running {}".format(program))
                except (subprocess.CalledProcessError, IOError):
                    print("Single k-mer assembly failed for locus when running {}".format(program))

    # Merge assemblies
    #for locus in locus_names:
    for locus in locus_names:
        print("Merging assemblies for ##{}##".format(locus))
        #oases_output = "{}/Oases_output/{}_{}".format(output_dir, cell_name, locus)
        oases_output_path = "{}/Oases_output/{}/MergedAssembly".format(output_dir, locus)
        #print("\n\n\n\n\n\n")
        #print("oases output path")
        #print(oases_output_path)
        #print("locus transcript path")
        #print(locus_transcript_path)
        locus_transcript_path = []
        for length in hash_lengths:
            #if len(locus_transcript_path) == 0:
                #locus_transcript_path = "{}/Oases_output/{}/directory{}/transcripts.fa".format(output_dir, locus, length)
            to_append = "{}/Oases_output/{}/directory{}/transcripts.fa".format(output_dir, locus, length)
            locus_transcript_path.append(to_append) 
        
        #locus_transcript_path = "{}/Oases_output/{}/directory*/transcripts.fa".format(output_dir, locus)
        #print(locus_transcript_path)
        commands = []
        if not single_end:
            
            velveth_command = [velveth, oases_output_path, '27', '-long']
            for i in range(len(locus_transcript_path)):
                to_append = locus_transcript_path[i]
                velveth_command.append(to_append)
            #print(velveth_command)
                
            velvetg_command = [velvetg, oases_output_path, '-read_trkg', 'yes', '-conserveLong', 'yes', '-exp_cov', 'auto']
            oases_command = [oases, oases_output_path, '-merge', 'yes']
            commands.append(velveth_command)
            commands.append(velvetg_command)
            commands.append(oases_command)

        for command in commands:
            if command == velveth_command:
                program = "velveth"
            elif command == velvetg_command:
                program = "velvetg"
            elif command == oases_command:
                program = "oases"
            print("Running {}".format(program))
            try:
                subprocess.check_call(command)
                """shutil.move('{}/Oases_output/pairedEnd*'.format(output_dir, cell_name, locus),
                        '{}/Oases_output/{}_{}.Oases.fasta'.format(output_dir, cell_name, locus))"""
            except (subprocess.CalledProcessError, IOError):
                print("Merging of assemblies failed for locus when running {}".format(program))   



def run_IgBlast(igblast, receptor, loci, output_dir, cell_name, index_location, ig_seqtype, species,
                should_resume, assembler):
    print("##Running IgBLAST##")
    print ("Ig_seqtype:", ig_seqtype)
    species_mapper = {
        'Mmus': 'mouse',
        'Hsap': 'human'
    }

    igblast_species = species_mapper[species]
    initial_locus_names = ["_".join([receptor,x]) for x in loci]
    locus_names = copy.copy(initial_locus_names)
    if should_resume:
        for locus in initial_locus_names:
            igblast_out = "{output_dir}/IgBLAST_output/{cell_name}_{receptor}_{locus}.IgBLASTOut".format(
                                                        output_dir=output_dir,cell_name=cell_name, 
                                                        receptor=receptor, locus=locus)
            if (os.path.isfile(igblast_out) and os.path.getsize(igblast_out) > 0):
                locus_names.remove(locus)
                print("Resuming with existing IgBLAST output for {locus}".format(locus=locus))
        
        if len(locus_names) == 0:    
            return
    
    print("Performing IgBlast on {locus_names}".format(locus_names = locus_names))

    databases = {}
    for segment in ['V', 'D', 'J']:
        databases[segment] = "{}/{}_{}.fa".format(index_location, receptor, segment)

    # Lines below suppress Igblast warning about not having an auxliary file.
    # Taken from http://stackoverflow.com/questions/11269575/how-to-hide-output-of-subprocess-in-python-2-7
    DEVNULL = open(os.devnull, 'wb')

    if receptor == "BCR":
        num_alignments_V = '20'
        num_alignments_D = '3'
        num_alignments_J = '5'
    else:
        num_alignments_V = '5'
        num_alignments_D = '5'
        num_alignments_J = '5'
    
    
    for locus in locus_names:
        print("##{}##".format(locus))
        if assembler == "trinity":
            trinity_fasta = "{}/Trinity_output/{}_{}.Trinity.fasta".format(output_dir, cell_name, locus)
        elif assembler == "basic":
            if locus in ["BCR_K", "BCR_L"]:
                trinity_fasta = "{}/Basic_output/basic_L.fa".format(output_dir, cell_name)
            else:
                trinity_fasta = "{}/Basic_output/basic_H.fa".format(output_dir, cell_name)
        else:
            trinity_fasta = "{}/Oases_output/{}/MergedAssembly/transcripts.fa".format(output_dir, locus)

        ############### MOVE TO CONFIG FILE IF WORKS!! ####################
        auxiliary_data = "/nfs/team205/il5/software/human_gl.aux"
        if os.path.isfile(trinity_fasta):
            command = [igblast, '-germline_db_V', databases['V'], '-germline_db_J', databases['J'], '-germline_db_D', 
                        databases['D'], '-domain_system', 'imgt', '-organism', igblast_species,
                       '-ig_seqtype', ig_seqtype, '-auxiliary_data', auxiliary_data, '-show_translation', '-num_alignments_V', num_alignments_V,
                       '-num_alignments_D', num_alignments_D, '-num_alignments_J', num_alignments_J, '-outfmt', '7', '-query', trinity_fasta]
            
            if assembler == "basic":
                igblast_out = igblast_out = "{output_dir}/Basic_IgBLAST_output/{cell_name}_{locus}.IgBLASTOut".format(output_dir=output_dir,
                                                                                              cell_name=cell_name,
                                                                                              locus=locus)
            else:
                igblast_out = "{output_dir}/IgBLAST_output/{cell_name}_{locus}.IgBLASTOut".format(output_dir=output_dir,
                                                                                              cell_name=cell_name,
                                                                                              locus=locus)
            with open(igblast_out, 'w') as out:
                # print(" ").join(pipes.quote(s) for s in command)
                subprocess.check_call(command, stdout=out, stderr=DEVNULL)
    DEVNULL.close()



def run_Blast(blast, receptor, loci, output_dir, cell_name, index_location, species,
                should_resume, assembler):
    print("##Running BLAST##") 

    species_mapper = {
        'Mmus': 'mouse',
        'Hsap': 'human'
    }

    blast_species = species_mapper[species]
    initial_locus_names = ["_".join([receptor,x]) for x in loci]
    locus_names = copy.copy(initial_locus_names)

    print("Performing Blast on {locus_names}".format(locus_names = locus_names))

    databases = {}
    
    
    for segment in ['c', 'C']:
        databases[segment] = "{}/{}_{}.fa".format(index_location, receptor, segment)
    
    if (os.path.isfile("{}/{}_c.fa".format(index_location, receptor)) and os.path.getsize("{}/{}_c.fa".format(index_location, receptor)) > 0):
        database = databases['c']
    else:
        database = databases['C']

    # Lines below suppress Igblast warning about not having an auxliary file.
    # Taken from http://stackoverflow.com/questions/11269575/how-to-hide-output-of-subprocess-in-python-2-7
    DEVNULL = open(os.devnull, 'wb')
  
    for locus in locus_names:
        print("##{}##".format(locus))
        if assembler == "trinity":
            trinity_fasta = "{}/Trinity_output/{}_{}.Trinity.fasta".format(output_dir, cell_name, locus)
        elif assembler == "basic":
            if locus in ["BCR_K", "BCR_L"]:
                trinity_fasta = "{}/Basic_output/basic_L.fa".format(output_dir, cell_name)
            else:
                trinity_fasta = "{}/Basic_output/basic_H.fa".format(output_dir, cell_name)
        else: 
            trinity_fasta = "{}/Oases_output/{}/MergedAssembly/transcripts.fa".format(output_dir, locus)
        if os.path.isfile(trinity_fasta):
            command = [blast, '-db', database, '-evalue', '0.001',
                        '-num_alignments', '1', '-outfmt', '5', '-query', trinity_fasta]
            if assembler == "basic":
                blast_out = "{output_dir}/Basic_BLAST_output/{cell_name}_{locus}.xml".format(output_dir=output_dir, cell_name=cell_name,
                                                                                              locus=locus)
            else:
                blast_out = "{output_dir}/BLAST_output/{cell_name}_{locus}.xml".format(output_dir=output_dir,
                                                                                              cell_name=cell_name,
                                                                                              locus=locus)
            with open(blast_out, 'w') as out:
                # print(" ").join(pipes.quote(s) for s in command)
                subprocess.check_call(command, stdout=out, stderr=DEVNULL)

    DEVNULL.close()



def quantify_with_kallisto(kallisto, cell, output_dir, cell_name, kallisto_base_transcriptome, fastq1, fastq2,
                           ncores, should_resume, single_end, fragment_length, fragment_sd, receptor_name):
    print("##Running Kallisto##")
    if should_resume:
        if os.path.isfile("{}/expression_quantification/abundance.tsv".format(output_dir)):
            print("Resuming with existing Kallisto output")
            return

    print("##Making Kallisto indices##")
    kallisto_dirs = ['kallisto_index']
    for d in kallisto_dirs:
        tracerlib.io.makeOutputDir("{}/expression_quantification/{}".format(output_dir, d))
    fasta_filename = "{output_dir}/unfiltered_{receptor}_seqs/{cell_name}_{receptor}seqs.fa".format(output_dir=output_dir,
                                                                                      cell_name=cell_name, receptor = receptor_name)
    fasta_file = open(fasta_filename, 'w')
    fasta_file.write(cell.get_fasta_string())
    fasta_file.close()

    output_transcriptome = "{}/expression_quantification/kallisto_index/{}_transcriptome.fa".format(output_dir,
                                                                                                    cell_name)
    with open(output_transcriptome, 'w') as outfile:
        for fname in [kallisto_base_transcriptome, fasta_filename]:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    idx_file = "{}/expression_quantification/kallisto_index/{}_transcriptome.idx".format(output_dir, cell_name)

    index_command = [kallisto, 'index', '-i', idx_file, output_transcriptome]
    subprocess.check_call(index_command)
    print("##Quantifying with Kallisto##")

    if not single_end:
        if not fragment_length:
            kallisto_command = [kallisto, 'quant', '-i', idx_file, '-t', ncores, '-o',
                                "{}/expression_quantification".format(output_dir), fastq1, fastq2]
        else:
            kallisto_command = [kallisto, 'quant', '-i', idx_file, '-t', ncores, '-l', fragment_length, '-o',
                                "{}/expression_quantification".format(output_dir), fastq1, fastq2]
    else:
        kallisto_command = [kallisto, 'quant', '-i', idx_file, '-t', ncores, '--single', '-l', fragment_length,
                            '-s', fragment_sd, '-o', "{}/expression_quantification".format(output_dir), fastq1]
    subprocess.check_call(kallisto_command)

    # delete index file because it's huge and unecessary. Delete transcriptome file
    # os.remove(idx_file)
    # os.remove(output_transcriptome)
    shutil.rmtree("{}/expression_quantification/kallisto_index/".format(output_dir))


def run_changeo(changeo, locus, outdir, species):
    
    # Set model to Hamming distance if species is not Mmus or Hsap
    if species == "Mmus":
        model = "m1n"
        dist = "0.2"
    elif species == "Hsap":
        model = "hs5f"
        dist = "0.2"
    else:
        model = "ham"
        dist = "0.2"

    
    changeo_input = "{}/changeo_input_{}.tab".format(outdir, locus)
    if os.path.isfile(changeo_input) and os.path.getsize(changeo_input) > 0:
        command = [changeo, "bygroup", '-d', changeo_input, '--mode', 'gene', '--act', 'set', 
                        '--model', model, '--dist', dist, '--sf', "JUNCTION", '--norm', 'len']

            #changeo_out = "{}/changeo_input_{}_clone-pass.tab".format(outdir, locus)
            #with open(changeo_result, 'w') as out:
                # print(" ").join(pipes.quote(s) for s in command)
        subprocess.check_call(command)

def run_changeo_clonal_alignment(changeo, locus, outdir, species):

    Command = "python DefineClones.py bygroup -d {changeo_input_file} --mode gene --act set --model m1n --dist 0.02 --sf JUNCTION"
    # Set model to Hamming distance if species is not Mmus or Hsap
    if species == "Mmus":
        model = "m1n"
        dist = "1.5"
    elif species == "Hsap":
        model = "hs5f"
        dist = "0.02"
    else:
        model = "ham"
        dist = "0.02"


    changeo_input = "{}/changeo_clonal_alignment_input_{}.tab".format(outdir, locus)
    if os.path.isfile(changeo_input):
        command = [changeo, "bygroup", '-d', changeo_input, '--mode', 'gene', '--act', 'set',
                        '--model', model, '--dist', dist, '--sf', "SEQUENCE_VDJ", '--norm', 'len', '-f', 'CLONE_GROUP']

            #changeo_out = "{}/changeo_input_{}_clone-pass.tab".format(outdir, locus)
            #with open(changeo_result, 'w') as out:
                # print(" ").join(pipes.quote(s) for s in command)
        subprocess.check_call(command)


def run_muscle(muscle, locus, outdir, species, group):
    muscle_input = "{}/muscle_input_{}_{}.fa".format(outdir, locus, group)
    muscle_fasta_out = "{}/muscle_out_{}_{}.afa".format(outdir, locus, group)
    muscle_clw_out = "{}/muscle_out_{}_{}.aln".format(outdir, locus, group)
    
    #changeo_input = "{}/changeo_input_{}.tab".format(outdir, locus)
    if os.path.isfile(muscle_input) and os.path.getsize(muscle_input)>0:
        # seqtype must be protein to allow for substitution matrix!
        command = [muscle, "-in", muscle_input, '-fastaout', muscle_fasta_out, '-clwout', muscle_clw_out, '-quiet'] 
        #'-matrix', matrix, '-seqtype', 'protein'

            #changeo_out = "{}/changeo_input_{}_clone-pass.tab".format(outdir, locus)
            #with open(changeo_result, 'w') as out:
                # print(" ").join(pipes.quote(s) for s in command)
        subprocess.check_call(command)

