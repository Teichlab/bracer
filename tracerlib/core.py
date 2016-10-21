import re
from collections import Counter, defaultdict

import six
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import pdb
from tracerlib import tracer_func


class Cell(object):

    """Class to describe T cells containing A and B loci and B cells containing H, K and L loci"""

    def __init__(self, cell_name, recombinants, is_empty=False, species="Mmus", 
                    receptor=None, loci=None):
        
        self.name = cell_name
        self.species = species
        self.recombinants = self._process_recombinants(recombinants, receptor, loci)
        
        self.is_empty = self._check_is_empty()
        
        self.isotype = self.determine_isotype(loci, receptor, self.recombinants)
        self.bgcolor = self.assign_bgcolor(species, self.isotype)
        self.changeodict = self.get_changeo_db_for_locus(self.recombinants, receptor, loci)
        for l in loci:
            print(self.changeodict[l])
        #self.cdr3_comparisons = {'A': None, 'B': None, 'mean_both': None}
        #invariant_types = []
        #if invariant_cells is not None:
        #    for ic in invariant_cells:
        #        itype = ic.check_for_match(self)
        #        if itype is not None:
        #            invariant_types.append(itype)
        

        #self.is_inkt = self._check_if_inkt()
    
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


    def get_changeo_db_for_locus(self, recombinants, receptor, loci):
        changeodict = defaultdict(dict)
        for l in loci:
            changeodict[l] = None
            changeo_string = ""
            recombinants = self.recombinants[receptor][l]
            for rec in recombinants:
                string = rec.create_changeo_db_string()
                string = self.name + "_" + string
                changeo_string += string + "\n"
            changeodict[l] = changeo_string
        return(changeodict)
            



    def determine_isotype(self, loci, receptor, recombinants):
        
        H_recombinants = self.recombinants[receptor]["H"]
        
        isotype_list = []
        for recombinant in H_recombinants:
            if recombinant.productive == True:
                C_gene = recombinant.C_gene
                isotype_list.append(C_gene)
        if len(isotype_list) == 0:
            isotype = None
        elif len(isotype_list) == 1:
            if isotype_list[0] is None:
                isotype = None
            else:
                isotype = isotype_list[0].split("*")[0]
        
        else:
            if isotype_list[0] == isotype_list[1]:
                isotype = isotype_list[0].split("*")[0]
            else:
                isotype = None 
        return(isotype)

    def assign_bgcolor(self, species, isotype):
        """Assigns bgcolor for cell according to isotype (B cells)"""
        
        isotype_bgcolors = {"IGHM":'#99e699', "IGHD":'#66a3ff', "IGHA":'#b366ff', "IGHE":'#ffff66', "IGHG1":'#b30000', "IGHG2A":'#e60000', "IGHG2B":'#ff3333', "IGHG2C":'#ff6666', "IGHG3":'#ffb3b3'}

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
        

    #def _check_if_inkt(self):
    #    A_recombs = self.getMainRecombinantIdentifiersForLocus("A")
    #    inkt_ident = False
    #    for recomb in A_recombs:
    #        for invar_seq in self.invariant_seqs:
    #            if invar_seq['V'] in recomb and invar_seq['J'] in recomb:
    #                inkt_ident = recomb
    #    return (inkt_ident)

    #def reset_cdr3_comparisons(self):
    #    self.cdr3_comparisons = {'A': None, 'B': None, 'mean_both': None}

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

    #def getAllRecombinantCDR3ForLocus(self, locus):
    #    recombinants = self.all_recombinants[locus]
    #    identifier_list = set()
    #    if recombinants is not None:
    #        for recombinant in recombinants:
    #            cdr3 = str(recombinant.cdr3)
    #            if "Couldn't" not in cdr3:
    #                identifier_list.add(cdr3)
    #    return (identifier_list)

    def html_style_label_dna(self, receptor, loci, colours):
        #colours = {'A': {'productive': '#E41A1C', 'non-productive': "#ff8c8e"},
        #           'B': {'productive': '#377eb8', 'non-productive': "#95c1e5"},
        #           'G': {'productive': '#4daf4a', 'non-productive': "#aee5ac"},
        #           'D': {'productive': '#984ea3', 'non-productive': "#deace5"}}
        #locus_names = ['A', 'B', 'G', 'D']
        

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
        # return(self.name)

    def html_style_label_for_circles(self, receptor, loci, colours):
        
        #colours = {'A': {'productive': '#E41A1C', 'non-productive': "#ff8c8e"},
        #           'B': {'productive': '#377eb8', 'non-productive': "#95c1e5"},
        #           'G': {'productive': '#4daf4a', 'non-productive': "#aee5ac"},
        #           'D': {'productive': '#984ea3', 'non-productive': "#deace5"}}
        #locus_names = ['A', 'B', 'G', 'D']
        
        
        
        recombinants = dict()
        final_string = '<<table cellspacing="6px" border="0" cellborder="0">'
        # final_string = "<"
        for locus, recombinant_list in six.iteritems(self.recombinants[receptor]):
            recombinant_set = set()
            if recombinant_list is not None:
                for recombinant in recombinant_list:
                    if recombinant.productive:
                        i = 0
                    else:
                        i = 1
                    recombinant_set.add(
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
                                                                            receptor=receptor, locus=locus,
                                                                            identifier=rec.identifier)
                        seq = rec.dna_seq
                        seq_string.append("\n".join([name, seq]))
        
        #for locus, recombinants in six.iteritems(self.all_recombinants):
        #    if recombinants is not None:
        #        for rec in recombinants:
        #            name = ">TCR|{contig_name}|{identifier}".format(contig_name=rec.contig_name,
        #                                                            identifier=rec.identifier)
        #            seq = rec.dna_seq
        #            seq_string.append("\n".join([name, seq]))
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
 
    """def get_all_cdr3_lengths(self, receptor, locus):
        recs = self.recombinants[receptor][locus]
        lengths = []
        if recs is not None:
            for rec in recs:
                lengths.append(len(rec.cdr3))
        return (lengths)"""
 
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
                 imgt_reconstructed_seq, has_D, output_dir, full_length, query_length, V_genes, cdr3):
        self.contig_name = contig_name
        self.locus = locus
        self.identifier = identifier
        self.all_poss_identifiers = all_poss_identifiers
        self.productive = productive
        self.TPM = TPM
        self.dna_seq = dna_seq
        #self.cdr3 = self._get_cdr3(dna_seq)
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
        self.C_gene = self.get_C_gene()
        self.full_length = full_length
        self.query_length = query_length
        self.V_genes = V_genes
        self.cdr3_in_frame = self.is_cdr3_in_frame(self.cdr3)       
        
    def __str__(self):
        return ("{} {} {} {}".format(self.identifier, self.productive, self.TPM))

    def _get_cdr3(self, dna_seq):
        
        aaseq = Seq(str(dna_seq), generic_dna).translate()
        # Specify first amino acid in conserved motif according to receptor and locus
        if self.locus in ["BCR_H", "H"]:
            motif_start = "W"
        else:
            motif_start = "F"
        motif = motif_start + "G.G"

        if re.findall(motif, str(aaseq)) and re.findall('C', str(aaseq)):
            indices = [i for i, x in enumerate(aaseq) if x == 'C']
            upper = str(aaseq).find(re.findall(motif, str(aaseq))[0])
            lower = False
            for i in indices:
                if i < upper:
                    lower = i
            # If motif not found, allow to search for "FSDG" in kappa sequences (present in IGKJ3)
            if lower == False:
                if self.locus in ["BCR_K", "K"]:
                    motif = "FSDG"
                    upper = str(aaseq).find(re.findall(motif, str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i

            if lower:
                cdr3 = aaseq[lower:upper + 4]
            else:
                cdr3 = "Couldn't find conserved cysteine"
        elif re.findall(motif, str(aaseq)):
            cdr3 = "Couldn't find conserved cysteine"
        elif re.findall('C', str(aaseq)):
            cdr3 = "Couldn't find {}GXG".format(motif_start)
        else:
            cdr3 = "Couldn't find either conserved boundary"

        return (cdr3)
 
    def is_cdr3_in_frame(self, cdr3):
        if len(cdr3) % 3 == 0:
            cdr3_in_frame = True
        else:
            cdr3_in_frame = False
        return (cdr3_in_frame)
       
    def get_C_gene(self):
        
        locus = self.locus.split("_")[1]
        blast_summary_file = "{output_dir}/BLAST_output/blastsummary_{locus}.txt".format(output_dir=self.output_dir, locus=locus)

        store_details = False
        C_gene = None
        with open(blast_summary_file, 'r') as input:
            for line in input:
                if line.startswith("##{contig_name}##".format(contig_name=self.contig_name)) or line.startswith("##reversed|{contig_name}##".format(contig_name=self.contig_name)):
                    store_details = True
                elif store_details == True:
                    if line.startswith("C segment"):
                        C_gene = line.split("\t")[1]
                        store_details = False
                    elif line.startswith("No C segment found"):
                        C_gene = None
        return (C_gene)


    def get_summary(self):
        summary_string = "##{contig_name}##\n".format(contig_name=self.contig_name)

        locus = self.locus.split("_")[1]
        blast_summary_file = "{output_dir}/BLAST_output/blastsummary_{locus}.txt".format(output_dir=self.output_dir, locus=locus)
       
 
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
            segments_string = "V segment:\t{V_segment}\nD segment:\t{D_segment}\nJ segment:\t{J_segment}\nC segment:\t{C_segment}\n".format(
                V_segment=V_segment, D_segment=D_segment, J_segment=J_segment, C_segment=C_segment)

        else:
            V_segment = self.summary[0]
            J_segment = self.summary[1]
            C_segment = C_gene
            segments_string = "V segment:\t{V_segment}\nJ segment:\t{J_segment}\nC segment:\t{C_segment}\n".format(
                V_segment=V_segment, J_segment=J_segment, C_segment=C_segment)


        summary_string += segments_string
        summary_string += "ID:\t{}\n".format(self.identifier)
        summary_string += "TPM:\t{TPM}\nProductive:\t{productive}\nStop codon:\t{stop_codon}\nIn frame:\t{in_frame}\nFull length:\t{full_length}\nSequence length: \t{query_length}\nPossible V genes:\t{V_genes}\n\n".format(
            TPM=self.TPM, productive=self.productive, stop_codon=self.stop_codon, in_frame=self.in_frame, full_length=self.full_length, query_length=self.query_length, V_genes=self.V_genes)

        summary_string += 'Segment\tquery_id\tsubject_id\t% identity\talignment length\tmismatches\tgap opens\tgaps\tq start\tq end\ts start\ts end\te value\tbit score\n'
        for line in self.hit_table:
            summary_string = summary_string + "\t".join(line) + "\n"
        
        if find_C_gene == True:
            store_details = False
            with open(blast_summary_file, 'r') as input:
                for line in input:
                    if line.startswith("##{contig_name}##".format(contig_name=self.contig_name)) or line.startswith("##reversed|{contig_name}##".format(contig_name=self.contig_name)):
                        store_details = True
                    elif store_details == True:
                        if line.startswith("C\t"):
                            summary_string = summary_string + line + "\n"
                            store_details = False
                
        return (summary_string)


    def create_changeo_db_string(self):
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
            junc_string = "".join(self.junction_details)
            junc_string = tracer_func.remove_NA(junc_string)
            junc_string = junc_string.split("(")
            junc_string = "".join(junc_string)
            junc_string = junc_string.split(")")
            junction = "".join(junc_string)
            junction_length = int(len(junction))
            changeo_db_string = "{}_{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.contig_name, self.identifier, V_call, D_call, J_call, sequence_vdj, junction_length, junction)

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
                
                
        
        
