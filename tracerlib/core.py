import re
from collections import Counter, defaultdict

import six
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import pdb
from tracerlib import tracer_func


class Cell(object):

    """Class to describe T cells containing A and B loci and B cells containing 
    H, K and L loci"""

    def __init__(self, cell_name, recombinants, is_empty=False, species="Mmus", 
                    receptor=None, loci=None):
        
        self.name = cell_name
        self.species = species
        self.receptor = receptor
        self.recombinants = self._process_recombinants(recombinants, receptor, loci)
        self.is_empty = self._check_is_empty()
        self.two_most_common = self.find_two_most_common()
        self.ranked_recs = self.rank_recombinants()
        self.bgcolor = None
        self.isotype = None
        self.changeodict = self.get_changeo_db_for_locus(receptor, loci)
        self.detailed_identifier_dict = self.create_detailed_identifier_dict()
        self.cdr3_dict = self.find_recs_with_identical_cdr3()
        self.has_excess_recombinants = self.has_excess_recombinants()
    
    def _process_recombinants(self, recombinants, receptor, loci):
        recombinant_dict = defaultdict(dict)
        if recombinants is not None:
            for r_name, r in six.iteritems(recombinants):
                r_name = r_name.split("_")
                receptor = r_name[0]
                locus = r_name[1]
                recombinant_dict[receptor][locus] = r
            
        #normalise this to put None in cases where no receptors found
        for l in loci:
            if l not in recombinant_dict[receptor]:
                recombinant_dict[receptor][l] = None
        return dict(recombinant_dict)


    def get_changeo_db_for_locus(self, receptor, loci):
        changeodict = defaultdict(dict)
        for l in loci:
            changeodict[l] = None
            changeo_string = ""
            recombinants = self.recombinants[receptor][l]
            for rec in recombinants:
                if rec.productive:
                    string = rec.create_changeo_db_string()
                    string = self.name + "_" + string
                    changeo_string += string + "\n"
                else:
                    changeo_string = changeo_string
            changeodict[l] = changeo_string
        return(changeodict)
            

    def determine_isotype(self, ranked_recs):
        """Determines isotype of B cell based on isotype of the two highest ranking 
        heavy chain recs if productive"""

        ranked_H = ranked_recs["H"][0:2]
        H_recombinants = self.recombinants[self.receptor]["H"]
        isotype = None
        isotype_list = []
        for ranked_rec in ranked_H:
            for recombinant in H_recombinants:
                if recombinant.contig_name == ranked_rec and recombinant.productive:
                    C_gene = recombinant.C_gene
                    if C_gene is not None:
                        if "*" in C_gene:
                            C_gene = C_gene.split("*")[0]
                        isotype_list.append(C_gene)
        if len(isotype_list) == 0:
            isotype = None
        elif len(isotype_list) == 1:
            isotype = isotype_list[0]
        elif len(isotype_list) == 2:
            iso_1 = isotype_list[0]
            iso_2 = isotype_list[1]
            if iso_1 is "IGHDM" or iso_2 is "IGHDM":
                isotype = "IGHDM"
            elif iso_1 == iso_2:
                isotype = iso_1
            elif iso_1 is not None and iso_2 is not None:
                if iso_1 in ["IGHM", "IGHD"] and iso_2 in ["IGHM", "IGHD"]:
                    isotype = "IGHDM"
                else:
                    isotype = iso_1
            elif iso_1 is None:
                isotype = iso_2
            elif iso_2 is None:
                isotype = iso_1
         
        return(isotype)

    def assign_bgcolor(self, isotype):
        """Assigns bgcolor for cell according to isotype (B cells)"""

        if self.species == "Mmus":
            isotype_bgcolors = {"IGHD":'#e6f7ff', "IGHM":'#e5ffcc', 
                "IGHA":'#ffe6e6', "IGHE":'#ffffcc', "IGHG1":'#f1e6ff', 
                "IGHG2A":'#e2ccff', "IGHG2B":'#d4b3ff', "IGHG2C":'#c599ff',
                "IGHG3":'a866ff', "IGHDM":'#b3ffff'}
        else:
            # Isotype background colours for Hsap
            isotype_bgcolors = {"IGHD":'#e6f7ff', "IGHM":'#e5ffcc', 
                "IGHA1":'#ffe6e6', "IGHA2":'#ffcccc', "IGHE":'#ffffcc', 
                "IGHG1":'#f1e6ff', "IGHG2":'#e2ccff', "IGHG3":'#d4b3ff', 
                "IGHG4":'#c599ff', "IGHDM":'#b3ffff'}

        if isotype is None:
            bgcolor = None
        else:
            bgcolor = isotype_bgcolors[isotype]

        return (bgcolor)
    
                
    def _check_is_empty(self):
        if (self.recombinants is None or len(self.recombinants) == 0):
            return True
        else:
            return False
    
    def missing_loci_of_interest(self, receptor_name, loci):
        recombinants = self.recombinants[receptor_name]
        loci_of_interest = set(loci)
        loci_in_cell = set()
        for l in loci:
            if l in recombinants and (recombinants[l] is not None and len(recombinants[l])>0):
                loci_in_cell.add(l)
        if len(loci_of_interest.intersection(loci_in_cell)) == 0:
            return True
        else:
            return False
        
    def getAllRecombinantIdentifiersForLocus(self, locus):
        recombinants = self.all_recombinants[locus]
        identifier_list = set()
        if recombinants is not None:
            for recombinant in recombinants:
                all_possible_recombinant_identifiers = recombinant.all_poss_identifiers
                for identifier in all_possible_recombinant_identifiers:
                    identifier_list.add(identifier)
        return (identifier_list)

    def getMainRecombinantIdentifiersForLocus(self, receptor_name, locus):
        recombinants = self.recombinants[receptor_name][locus]
        identifier_list = set()
        if recombinants is not None:
            for recombinant in recombinants:
                identifier_list.add(recombinant.identifier)
        return identifier_list


    def html_style_label_dna(self, receptor, loci, colours):

        recombinants = dict()
        final_string = '<<FONT POINT-SIZE="16"><B>' + self.name + "</B></FONT>"
        for locus, recombinant_list in six.iteritems(self.recombinants[receptor]):
            recombinant_set = set()
            if recombinant_list is not None:
                for recombinant in recombinant_list:
                    if recombinant.productive:
                        i = 0
                    else:
                        i = 1
                    recombinant_set.add("<BR/>" + '<FONT COLOR = "{}">'.format(
                        colours[receptor][locus][i]) + recombinant.identifier + '</FONT>')

                recombinants[locus] = recombinant_set
        for locus in loci:
            if locus in recombinants.keys():
                id_string = "".join(recombinants[locus])
                final_string = final_string + id_string
        final_string = final_string + ">"
        return (final_string)

    def html_style_label_for_circles(self, receptor, loci, colours):
        
        recombinants = dict()
        final_string = '<<table cellspacing="6px" border="0" cellborder="0">'
        for locus, recombinant_list in six.iteritems(self.recombinants[receptor]):
            recombinant_set = list()
            if recombinant_list is not None:
                for recombinant in recombinant_list:
                    if recombinant.productive:
                        i = 0
                    else:
                        i = 1
                    recombinant_set.append(
                        '<tr><td height="10" width="40" bgcolor="{}"></td></tr>'.format(colours[receptor][locus][i]))

                recombinants[locus] = recombinant_set
        strings = []
        for locus in loci:
            if locus in recombinants.keys():
                strings.append("".join(recombinants[locus]))

        id_string = "".join(strings)
        final_string = final_string + id_string
        final_string = final_string + "</table>>"
        return (final_string)

    def __str__(self):
        return (self.name)

    def full_description(self):
        # pdb.set_trace()
        return_list = [self.name, '#TCRA#']

        if not self.A_recombinants is None:
            for recombinant in self.A_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No TCRA recombinants")

        return_list.append('\n#TCRB#')
        if not self.B_recombinants is None:
            for recombinant in self.B_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No TCRB recombinants")

        return_list.append('\n#TCRG#')
        if not self.G_recombinants is None:
            for recombinant in self.G_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No TCRG recombinants")

        return_list.append('\n#TCRD#')
        if not self.D_recombinants is None:
            for recombinant in self.D_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No TCRD recombinants")

        return_list.append('\n#BCRH#')
        if not self.H_recombinants is None:
            for recombinant in self.H_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No BCRH recombinants")

        return_list.append('\n#BCRK#')
        if not self.K_recombinants is None:
            for recombinant in self.K_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No BCRK recombinants")

        return_list.append('\n#BCRL#')
        if not self.L_recombinants is None:
            for recombinant in self.L_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No BCRL recombinants")


        return ("\n".join(return_list))

    def get_fasta_string(self):
        seq_string = []
        
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None:
                    for rec in recombinants:
                        name = ">TRACER|{receptor}|{locus}|{contig_name}|{identifier}".format(contig_name=rec.contig_name,
                                                                receptor=receptor, locus=locus, identifier=rec.identifier)
                        seq = rec.dna_seq
                        seq_string.append("\n".join([name, seq]))
        
        return ("\n".join(seq_string + ["\n"]))



    def summarise_productivity(self, receptor, locus):
        if (self.recombinants is None or locus not in self.recombinants[receptor] or 
            self.recombinants[receptor][locus] is None):
            return("0/0")
        else:
            recs = self.recombinants[receptor][locus]
            prod_count = 0
            total_count = len(recs)
            for rec in recs:
                if rec.productive:
                    prod_count += 1
            return ("{}/{}".format(prod_count, total_count))
                


    def find_recs_with_identical_cdr3(self):
        for receptor, locus_dict in six.iteritems(self.recombinants):
            cdr3_dict = dict()
            for locus, recombinants in six.iteritems(locus_dict):
                cdr3_dict[locus] = dict()
                if recombinants is not None:
                    if len(recombinants) > 1:
                        for rec in recombinants:
                            cdr3 = rec.cdr3
                            if not cdr3 in cdr3_dict[locus].keys():
                                cdr3_dict[locus][cdr3] = [rec.contig_name]
                            else:
                                cdr3_dict[locus][cdr3].append(rec.contig_name)
        return (cdr3_dict)

    def rank_recombinants(self):
        """Ranks recombinants from highest TPM to lowest in case of more than two 
        recombinants for a locus"""
        ranked_recs = dict()
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                most_common = []
                if recombinants is not None:
                    if len(recombinants) > 0:
                        TPM_ranks = Counter()
                        for rec in recombinants:
                            TPM_ranks.update({rec.contig_name: rec.TPM})
                        most_common = [x[0] for x in TPM_ranks.most_common()]
                        
                    ranked_recs[locus] = most_common
        return (ranked_recs)



    def find_two_most_common(self):
        ranked_recs = self.rank_recombinants()
        two_most_common_dict = dict()
        for locus, recs in six.iteritems(ranked_recs):
            if len(recs) > 1:
                two_recs = recs[:2]
            elif len(recs) == 1:
                two_recs = recs
            else:
                two_recs = []
            two_most_common_dict[locus] = two_recs
        return (two_most_common_dict)


    def print_recombinants(self):
        print_dict = dict()
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                print_dict[locus] = dict()
                if recombinants is not None:
                    for rec in recombinants:
                        tpm = rec.TPM
                        if rec.C_gene is not None:
                            C_gene = rec.C_gene.split("*")[0]
                        else:
                            C_gene = rec.C_gene
                        print_dict[locus][rec] = [rec.contig_name, tpm, rec.cdr3, 
                                                          C_gene, rec.productive]

        return (print_dict)

                           
    def create_detailed_identifier_dict(self):
        detailed_identifier_dict = dict()
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                detailed_identifier_dict[locus] = dict()
                if recombinants is not None:
                    for rec in recombinants:
                        
                        detailed_identifier_dict[locus][rec.contig_name] = rec.detailed_identifier

        return (detailed_identifier_dict)

    def create_tpm_dict(self, recombinants, loci):
        tpm_dict = dict()
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                tpm_dict[locus] = dict()
                if recombinants is not None:
                    for rec in recombinants:
                        tpm = rec.TPM
                        tpm_dict[locus][rec.contig_name] = tpm

        return (tpm_dict)


    def filter_recombinants(self):
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None:
                    if len(recombinants) > 2:
                        TPM_ranks = Counter()
                        for rec in recombinants:
                            TPM_ranks.update({rec.contig_name: rec.TPM})
                        two_most_common = [x[0] for x in TPM_ranks.most_common(2)]
                        to_remove = []
                        for rec in recombinants:
                             if rec.contig_name not in two_most_common:
                                 to_remove.append(rec)
                        for rec in to_remove:
                             self.recombinants[receptor][locus].remove(rec)


 
    def count_productive_recombinants(self, receptor, locus):
        recs = self.recombinants[receptor][locus]
        count = 0
        if recs is not None:
            for rec in recs:
                if rec.productive:
                    count += 1
        return (count)

    def count_total_recombinants(self, receptor, locus):
        recs = self.recombinants[receptor][locus]
        count = 0
        if recs is not None:
            count = len(recs)
        return (count)

    def count_full_length_recombinants(self, receptor, locus):
        recs = self.recombinants[receptor][locus]
        count = 0
        if recs is not None:
            for rec in recs:
                if rec.full_length:
                    count += 1
        return (count)

    def count_productive_full_length_recombinants(self, receptor, locus):
        recs = self.recombinants[receptor][locus]
        count = 0
        if recs is not None:
            for rec in recs:
                if rec.full_length and rec.productive:
                    count += 1
        return (count)


    def get_trinity_lengths(self, receptor, locus):
        recs = self.recombinants[receptor][locus]
        lengths = []
        if recs is not None:
            for rec in recs:
                lengths.append(len(rec.trinity_seq))
        return (lengths)
 
 
    def get_prod_cdr3_lengths(self, receptor, locus):
        recs = self.recombinants[receptor][locus]
        lengths = []
        if recs is not None:
            for rec in recs:
                if rec.cdr3 is not None and rec.productive:
                    length = len(rec.cdr3) - 5
                    lengths.append(length)
        return (lengths)
       
    def has_excess_recombinants(self, max_r=2):
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recs in six.iteritems(locus_dict):
                if recs is not None:
                    if len(recs) > max_r:
                        return(True)


class Recombinant(object):

    """Class to describe a recombined TCR or BCR locus as determined from the single-cell pipeline"""

    def __init__(self, contig_name, locus, identifier, all_poss_identifiers, productive, stop_codon, in_frame, TPM,
                 dna_seq, hit_table, summary, junction_details, best_VJ_names, alignment_summary, trinity_seq,
                 imgt_reconstructed_seq, has_D, output_dir, full_length, query_length, V_genes, J_genes, cdr3, 
                 assembler, C_gene, C_info_line, cdr3_seq):
        self.contig_name = contig_name
        self.locus = locus
        self.identifier = identifier
        self.all_poss_identifiers = all_poss_identifiers
        self.productive = productive
        self.TPM = TPM
        self.dna_seq = dna_seq
        self.assembler = assembler
        self.cdr3 = cdr3
        self.hit_table = hit_table
        self.summary = summary
        self.junction_details = junction_details
        self.best_VJ_names = best_VJ_names
        self.alignment_summary = alignment_summary
        self.in_frame = in_frame
        self.stop_codon = stop_codon
        self.trinity_seq = trinity_seq
        self.imgt_reconstructed_seq = imgt_reconstructed_seq
        self.has_D_segment = has_D
        self.output_dir = output_dir
        self.C_gene = C_gene
        self.C_gene_info = C_info_line
        self.J_gene = self.get_J_gene()
        self.full_length = full_length
        self.query_length = query_length
        self.V_genes = V_genes
        self.J_genes = J_genes
        self.cdr3_seq = cdr3_seq
        self.detailed_identifier = self.create_detailed_identifier(self.productive, 
                                         self.cdr3, self.C_gene, self.full_length)        

    def __str__(self):
        return ("{} {} {} {}".format(self.identifier, self.productive, self.TPM))

    def create_detailed_identifier(self, productive, cdr3, C_gene, full_length):
        if productive == True:
            productive = "productive"
        else:
            productive = "nonproductive"
        if C_gene == None:
            C_gene = "none"
        else:
            C_gene = C_gene.split("*")[0]
        if "Could" in cdr3:
            cdr3 = "none"
        else:
            cdr3 = str(cdr3)
        if full_length == True:
            full_length = "true"
        else:
            full_length = "false"
        detailed_identifier = "{}_{}*{}_{}".format(cdr3, productive, C_gene, full_length)
        return (detailed_identifier)

    def get_J_gene(self):
        if not self.has_D_segment:
            J_segment = self.summary[1]
        else:
            J_segment = self.summary[2]
        return (J_segment)

    def get_summary(self):
        summary_string = "##{contig_name}##\n".format(contig_name=self.contig_name)
        locus = self.locus.split("_")[1]
        if locus in ["H", "K", "L"]:
            find_C_gene = True
            C_gene = self.C_gene
        else:
            find_C_gene = False

        if not self.has_D_segment and find_C_gene == False:
            V_segment = self.summary[0]
            J_segment = self.summary[1]
            segments_string = "V segment:\t{V_segment}\nJ segment:\t{J_segment}\n".format(V_segment=V_segment,
                                                                                          J_segment=J_segment)
        elif find_C_gene == False:
            V_segment = self.summary[0]
            D_segment = self.summary[1]
            J_segment = self.summary[2]
            segments_string = "V segment:\t{V_segment}\nD segment:\t{D_segment}\nJ segment:\t{J_segment}\n".format(
                V_segment=V_segment, D_segment=D_segment, J_segment=J_segment)

        elif self.has_D_segment and find_C_gene == True:
            V_segment = self.summary[0]
            D_segment = self.summary[1]
            J_segment = self.summary[2]
            C_segment = C_gene
            segments_string = "V segment:\t{V_segment}\nD segment:\t{D_segment}\nJ segment:\t{J_segment}\
                \nC segment:\t{C_segment}\n".format(V_segment=V_segment, D_segment=D_segment, J_segment=J_segment, 
                                                                                              C_segment=C_segment)

        else:
            V_segment = self.summary[0]
            J_segment = self.summary[1]
            C_segment = C_gene
            segments_string = "V segment:\t{V_segment}\nJ segment:\t{J_segment}\nC segment:\t{C_segment}\n".format(
                                                     V_segment=V_segment, J_segment=J_segment, C_segment=C_segment)


        summary_string += segments_string
        summary_string += "ID:\t{}\n".format(self.identifier)
        summary_string += "TPM:\t{TPM}\nProductive:\t{productive}\nStop codon:\t{stop_codon}\nIn frame:\t{in_frame}\
            \nFull length:\t{full_length}\nSequence length: \t{query_length}\nPossible V genes:\t{V_genes}\\n\n".format(
            TPM=self.TPM, productive=self.productive, stop_codon=self.stop_codon, in_frame=self.in_frame, 
            full_length=self.full_length, query_length=self.query_length, V_genes=self.V_genes)

        summary_string += 'Segment\tquery_id\tsubject_id\t% identity\talignment length\tmismatches\tgap opens\tgaps\
            \tq start\tq end\ts start\ts end\te value\tbit score\n'

        for line in self.hit_table:
            summary_string = summary_string + "\t".join(line) + "\n"
        
        if find_C_gene == True and C_gene != None and self.C_gene_info is not None:
            summary_string = summary_string + self.C_gene_info + "\n"
                
        return (summary_string)


    def create_changeo_db_string(self):
        """For assessment of B cell clonality"""
        changeo_db_header = "SEQUENCE_ID\tV_CALL\tD_CALL\tJ_CALL\tSEQUENCE_VDJ\tJUNCTION_LENGTH\tJUNCTION"
        #Add sequence_ID at cell level to include cell name
        changeo_db_string = ""
        if self.productive:
            V_genes = self.V_genes
            V_call = ",".join(str(x) for x in V_genes)
            D_call = "None"
            if not self.has_D_segment:
                J_call = self.summary[1]
            else:
                J_call = self.summary[2]
            sequence_vdj = self.dna_seq
            J_genes = J_call.split(",")
            J_call = ""
            for J_gene in J_genes:
                J_gene = J_gene.split("*")[0]
                if len(J_call) == 0:
                    J_call = J_gene
                else:
                    if not J_gene in J_call:
                        J_call += ",{}".format(J_gene)

            # Replace junction with CDR3 sequence
            junction = self.cdr3_seq
            junction_length = int(len(junction))
            changeo_db_string = "{}_{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.contig_name, self.identifier, 
                                        V_call, D_call, J_call, sequence_vdj, junction_length, junction)

        return(changeo_db_string)

class Invar_cell(object):
    
    """Class to describe invariant cells and their specific sequences"""
    
    def __init__(self, d):
        self.name = d['cell_name']
        self.receptor_type = d['receptor_type']
        self.invariant_recombinants = d['recombinants']
        self.defining_locus = d['defining_locus']
        self.expected_string = self._get_expected_string()
        
    def check_for_match(self, cell, locus):
        found_identifiers = set()
        found_locus = False
        

        #check for expected recombinants for defining locus
        cell_recs = cell.recombinants[self.receptor_type][locus]
        invariant_recs = self.invariant_recombinants[locus]
        if cell_recs is not None:
            for rec in cell_recs:
                if rec.productive:
                    for ident in rec.all_poss_identifiers:
                        ident = ident.split("_")
                        v = ident[0]
                        j = ident[2]
                        for ivr in invariant_recs:
                            if (v in ivr['V'] or ivr['V']=='*') and (j in ivr['J'] or ivr['J']=='*'):
                                found_locus = True
                                found_identifiers.add("_".join(ident))
        
        return found_locus, found_identifiers
    
    def _get_expected_string(self):
        s = ""
        defining_recs = self.invariant_recombinants[self.defining_locus]
        r = defining_recs[0]
        s = s + "-".join([r['V'], r['J']])
        if len(defining_recs) > 1:
            for r in defining_recs[1:]:
                s = s + " | " + "-".join([r['V'], r['J']])
        
        for l in self.invariant_recombinants.keys():
            if not l==self.defining_locus:
                recs = self.invariant_recombinants[l]
                r = recs[0]
                s = s + "," +  "-".join([r['V'], r['J']])
                if len(recs) > 1:
                    for r in recs[1:]:
                        s = s + " | " + "-".join([r['V'], r['J']])
        
        return s
                
                
        
        
