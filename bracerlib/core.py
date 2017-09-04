import re
from collections import Counter, defaultdict

import six
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq


class Cell(object):

    """Class to describe B cells containing H, K and L loci"""

    def __init__(self, cell_name, recombinants, is_empty=False, species="Hsap", 
                loci=None):
        
        self.name = cell_name
        self.species = species
        self.recombinants = self._process_recombinants(recombinants, loci)
        self.is_empty = self._check_is_empty()
        self.ranked_recs = self.rank_recombinants()
        self.bgcolor = None
        self.isotype = None
        self.changeodict = self.get_changeo_db_for_locus(loci)
        self.has_excess_recombinants = self.has_excess_recombinants()
        self.databasedict = self.get_database_for_locus(loci)
    
    def _process_recombinants(self, recombinants, loci):
        recombinant_dict = defaultdict(dict)
        if recombinants is not None:
            for r_name, r in six.iteritems(recombinants):
                r_name = r_name.split("_")
                locus = r_name[1]
                recombinant_dict["BCR"][locus] = r
            
        #normalise this to put None in cases where no receptors found
        for l in loci:
            if l not in recombinant_dict["BCR"]:
                recombinant_dict["BCR"][l] = None
        return dict(recombinant_dict)


    def get_changeo_db_for_locus(self, loci):

        changeodict = defaultdict(dict)
        for l in loci:
            changeodict[l] = None
            changeo_string = ""
            recombinants = self.recombinants["BCR"][l]
            if not recombinants is None:
                for rec in recombinants:
                    if rec.productive:
                        string = rec.create_changeo_db_string()
                        string = self.name + "_" + string
                        changeo_string += string + "\n"
                    else:
                        changeo_string = changeo_string
            changeodict[l] = changeo_string

        return(changeodict)
            
    def get_database_for_locus(self, loci):

        databasedict = defaultdict(dict)
        for l in loci:
            databasedict[l] = None
            changeo_string = ""
            recombinants = self.recombinants["BCR"][l]
            if not recombinants is None:
                databasedict[l] = dict()
                for rec in recombinants:
                    string = rec.create_database_string()
                    string = "{}\t{}_".format(self.name, self.name) + string
                    databasedict[l][rec.contig_name] = string

        return(databasedict)
    

    def determine_isotype(self, ranked_recs):
        """Determines isotype of B cell based on isotype of the two highest ranking 
        heavy chain recs if productive"""

        ranked_H = ranked_recs["H"][0:2]
        H_recombinants = self.recombinants["BCR"]["H"]
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
        """Assigns bgcolor for cell according to isotype"""

        if self.species == "Mmus":
            isotype_bgcolors = {"IGHD":'#e6f7ff', "IGHM":'#e5ffcc', 
                "IGHA":'#ffe6e6', "IGHE":'#ffffcc', "IGHG1":'#f1e6ff', 
                "IGHG2A":'#e2ccff', "IGHG2B":'#d4b3ff', "IGHG2C":'#c599ff',
                "IGHG3":'#a866ff', "IGHDM":'#b3ffff'}
        elif self.species == "Hsap":
            isotype_bgcolors = {"IGHD":'#e6f7ff', "IGHM":'#e5ffcc', 
                "IGHA1":'#ffe6e6', "IGHA2":'#ffcccc', "IGHE":'#ffffcc', 
                "IGHG1":'#f1e6ff', "IGHG2":'#e2ccff', "IGHG3":'#d4b3ff', 
                "IGHG4":'#c599ff', "IGHDM":'#b3ffff'}
        else:
            isotype_bgcolors = {"IGHD":'#e6f7ff', "IGHM":'#e5ffcc', 
                "IGHA":'#ffe6e6', "IGHE":'#ffffcc', "IGHG":'#e2ccff', 
                "IGHDM":'#b3ffff'}

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
    
    def missing_loci_of_interest(self, loci):
        recombinants = self.recombinants["BCR"]
        loci_of_interest = set(loci)
        loci_in_cell = set()
        for l in loci:
            if l in recombinants and (recombinants[l] is not None \
                and len(recombinants[l])>0):
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

    def getMainRecombinantIdentifiersForLocus(self, locus):
        recombinants = self.recombinants["BCR"][locus]
        identifier_list = set()
        if recombinants is not None:
            for recombinant in recombinants:
                identifier_list.add(recombinant.identifier)
        return identifier_list


    def html_style_label_dna(self, loci, colours):

        recombinants = dict()
        final_string = '<<FONT POINT-SIZE="16"><B>' + self.name + "</B></FONT>"
        for locus, recombinant_list in six.iteritems(self.recombinants["BCR"]):
            recombinant_set = set()
            if recombinant_list is not None:
                for recombinant in recombinant_list:
                    if recombinant.productive:
                        i = 0
                    else:
                        i = 1
                    recombinant_set.add("<BR/>" + '<FONT COLOR = "{}">'.format(
                        colours["BCR"][locus][i]) + recombinant.identifier + '</FONT>')

                recombinants[locus] = recombinant_set
        for locus in loci:
            if locus in recombinants.keys():
                id_string = "".join(recombinants[locus])
                final_string = final_string + id_string
        final_string = final_string + ">"
        return (final_string)

    def html_style_label_for_circles(self, loci, colours):
        
        recombinants = dict()
        final_string = '<<table cellspacing="6px" border="0" cellborder="0">'
        for locus, recombinant_list in six.iteritems(self.recombinants["BCR"]):
            recombinant_set = list()
            if recombinant_list is not None:
                for recombinant in recombinant_list:
                    if recombinant.productive:
                        i = 0
                    else:
                        i = 1
                    recombinant_set.append(
                        '<tr><td height="10" width="40" bgcolor="{}"></td></tr>'.format(
                                                            colours["BCR"][locus][i]))

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
        return_list = [self.name, '#BCRH#']

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
                        name = ">BRACER|{receptor}|{locus}|{contig_name}|{identifier}".format(
                                contig_name=rec.contig_name, receptor=receptor, 
                                        locus=locus, identifier=rec.identifier)
                        seq = str(rec.dna_seq)
                        seq_string.append("\n".join([name, seq]))
        
        return ("\n".join(seq_string + ["\n"]))


    def summarise_productivity(self, locus):
        if (self.recombinants is None or locus not in self.recombinants["BCR"] or 
            self.recombinants["BCR"][locus] is None):
            return("0/0")
        else:
            recs = self.recombinants["BCR"][locus]
            prod_count = 0
            total_count = len(recs)
            for rec in recs:
                if rec.productive:
                    prod_count += 1
            return ("{}/{}".format(prod_count, total_count))


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

 
    def count_productive_recombinants(self, locus):
        recs = self.recombinants["BCR"][locus]
        count = 0
        if recs is not None:
            for rec in recs:
                if rec.productive:
                    count += 1
        return (count)

    def count_total_recombinants(self, locus):
        recs = self.recombinants["BCR"][locus]
        count = 0
        if recs is not None:
            count = len(recs)
        return (count)

    def count_full_length_recombinants(self, locus):
        recs = self.recombinants["BCR"][locus]
        count = 0
        if recs is not None:
            for rec in recs:
                if rec.full_length:
                    count += 1
        return (count)

    def count_productive_full_length_recombinants(self, locus):
        recs = self.recombinants["BCR"][locus]
        count = 0
        if recs is not None:
            for rec in recs:
                if rec.full_length and rec.productive:
                    count += 1
        return (count)


    def get_trinity_lengths(self, locus):
        recs = self.recombinants["BCR"][locus]
        lengths = []
        if recs is not None:
            for rec in recs:
                lengths.append(len(rec.trinity_seq))
        return (lengths)
 
       
    def has_excess_recombinants(self, max_r=2):
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recs in six.iteritems(locus_dict):
                if recs is not None:
                    if len(recs) > max_r:
                        return(True)


class Recombinant(object):

    """Class to describe a recombined BCR locus as determined from the 
    single-cell pipeline"""

    def __init__(self, contig_name, locus, identifier, all_poss_identifiers, 
                 productive, stop_codon, in_frame, TPM, dna_seq, hit_table, 
                 summary, junction_details, best_VJ_names, alignment_summary, 
                 trinity_seq, has_D, output_dir, full_length, query_length, 
                 V_genes, J_genes, cdr3, C_gene, C_info_line, cdr3_seq, 
                 junc_string, untrimmed_seq, gapped_db_string):

        self.contig_name = contig_name
        self.locus = locus
        self.identifier = identifier
        self.all_poss_identifiers = all_poss_identifiers
        self.productive = productive
        self.in_frame = in_frame
        self.stop_codon = stop_codon
        self.TPM = TPM
        self.hit_table = hit_table
        self.summary = summary
        self.junction_details = junction_details
        self.best_VJ_names = best_VJ_names
        self.alignment_summary = alignment_summary
        self.trinity_seq = trinity_seq
        self.untrimmed_seq = untrimmed_seq
        self.dna_seq = dna_seq
        self.cdr3 = cdr3
        self.cdr3_seq = cdr3_seq
        self.output_dir = output_dir
        self.V_genes = V_genes
        self.V_gene = self.get_V_gene()
        self.has_D_segment = has_D
        self.D_gene = self.get_D_gene()
        self.J_gene = self.get_J_gene()
        self.J_genes = J_genes
        self.C_gene = C_gene
        self.C_gene_info = C_info_line
        self.full_length = full_length
        self.query_length = query_length
        self.junc_string = junc_string
        self.gapped_db_string = gapped_db_string
        

    def __str__(self):
        return ("{} {} {} {}".format(self.identifier, self.productive, self.TPM))

    def get_V_gene(self):
        V_segment = self.summary[0]
        return(V_segment)

    def get_D_gene(self):
        if self.has_D_segment:
            D_segment = self.summary[1]
        else:
            D_segment = "None"
        return(D_segment)

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
            segments_string = "V segment:\t{V_segment}\nJ segment:\t{J_segment}\n".format(
                                V_segment=V_segment, J_segment=J_segment)

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
               V_segment=V_segment, D_segment=D_segment, J_segment=J_segment, 
               C_segment=C_segment)

        else:
            V_segment = self.summary[0]
            J_segment = self.summary[1]
            C_segment = C_gene
            segments_string = "V segment:\t{V_segment}\nJ segment:\t{J_segment}\nC segment:\t{C_segment}\n".format(
                            V_segment=V_segment, J_segment=J_segment, C_segment=C_segment)

        summary_string += segments_string
        summary_string += "ID:\t{}\n".format(self.identifier)
        summary_string += "TPM:\t{TPM}\nProductive:\t{productive}\nStop codon:\t{stop_codon}\n".format(
                                    TPM=self.TPM, productive=self.productive, stop_codon=self.stop_codon)
        summary_string += "In frame:\t{in_frame}\n""Full length:\t{full_length}\nSequence length: \t{query_length}\n".format(
                                        in_frame=self.in_frame, full_length=self.full_length, query_length=self.query_length)
        V_genes = ", ".join(self.V_genes)
        J_genes = ", ".join(self.J_genes)

        summary_string += "All possible V genes:\t{V_genes}\nAll possible J genes:\t{J_genes}\n\n".format(
                        V_genes=V_genes, J_genes=J_genes)

        summary_string += ("Segment\tquery_id\tsubject_id\t% identity"
            "\talignment length\tmismatches\tgap opens\tgaps\tq start"
            "\tq end\ts start\ts end\te value\tbit score\n")

        for line in self.hit_table:
            summary_string = summary_string + "\t".join(line) + "\n"
        
        if find_C_gene == True and C_gene != None and self.C_gene_info is not None:
            summary_string = summary_string + self.C_gene_info + "\n"
                
        return (summary_string)


    def create_changeo_db_string(self):
        """For assessment of clonality with Change-O"""

        changeo_db_header = "SEQUENCE_ID\tV_CALL\tD_CALL\tJ_CALL\tSEQUENCE_VDJ\tJUNCTION_LENGTH\tJUNCTION"

        #Add sequence_ID at cell level to include cell name
        changeo_db_string = ""
        if self.productive:
            V_genes = self.V_genes
            V_call = ",".join(str(x) for x in V_genes)
            D_call = "None"
            J_genes = self.J_genes
            J_call = ",".join(str(x) for x in J_genes)
            sequence_vdj = self.dna_seq

            # Replace JUNCTION with CDR3 sequence
            junction = self.cdr3_seq
            junction_length = int(len(junction))
            changeo_db_string = "{}_{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    self.contig_name, self.identifier, V_call, D_call, 
                    J_call, sequence_vdj, junction_length, junction)

        return(changeo_db_string)

    def create_database_string(self):
        """Creates string for recombinant to include in tab-delimited Change-O database file"""

        changeo_db_header = ("SEQUENCE_ID\tCONTIG_NAME\tLOCUS\tFUNCTIONAL\tIN_FRAME" +
            "\tSTOP\tINDELS\tV_CALL\tTOP_V_ALLELE\tD_CALL\tTOP_D_ALLELE\tJ_CALL\t" + 
            "TOP_J_ALLELE\tSEQUENCE_INPUT\tSEQUENCE_VDJ\tV_SEQ_START\tV_SEQ_LENGTH" + 
            "\tD_SEQ_START\tD_SEQ_LENGTH\tJ_SEQ_START\tJ_SEQ_LENGTH\t" +
            "\tCDR3\tCDR3_LENGTH\tJUNCTION\tJUNCTION_LENGTH\tTPM\tC_CALL")
        
        locus = self.locus.split("_")[1]

        #Create unique SEQUENCE_ID
        seq_id = "{}_{}".format(self.contig_name, locus)
        changeo_db_string = ""

        V_genes = self.V_genes
        V_call = ",".join(str(x) for x in V_genes)
        D_call = "None"
        J_genes = self.J_genes
        J_call = ",".join(str(x) for x in J_genes)
        sequence_vdj = self.dna_seq

        # Check if CDR3 sequence is detected
        cdr3 = self.cdr3_seq
        if cdr3 is not None:
            cdr3_length = int(len(cdr3))
        else:
            cdr3 = "None"
            cdr3_length = "None"
        
        # Extract IgBlast info for top hits
        found_V = False
        found_D = False
        found_J = False
        V_seq_start = "None"
        V_seq_length = "None"
        J_seq_start = "None"
        J_seq_length = "None"
        D_seq_start = "None"
        D_seq_length = "None"

        for hit in self.hit_table:
            info = hit
            segment = info[0]
            gap_opens = int(info[6])
            gaps = int(info[7])
            seq_start = info[8]
            align_length = info[4]
            if gap_opens > 0 or gaps > 0:
                indels = True
            else:
                indels = False
            if segment == "V" and found_V is False:
                found_V = True
                V_indels = indels
                V_seq_start = seq_start
                V_seq_length = align_length
            elif segment == "D" and found_D is False:
                found_D = True
                D_seq_start = seq_start
                D_seq_length = align_length
            elif segment == "J" and found_J is False:
                found_J = True
                J_indels = indels
                J_seq_start = seq_start
                J_seq_length = align_length

        if V_indels is True or J_indels is True:
            indels = True
        else:
            indels = False


        junction_length = int(len(self.junc_string))
        changeo_db_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            seq_id, self.contig_name, locus, self.productive, self.in_frame, 
            self.stop_codon, indels, V_call, self.V_gene, D_call, self.D_gene,
            J_call, self.J_gene, self.untrimmed_seq, sequence_vdj, V_seq_start, 
            V_seq_length, D_seq_start, D_seq_length, J_seq_start, J_seq_length,
            cdr3, cdr3_length, self.junc_string, junction_length, self.TPM, self.C_gene)
                
                
        return(changeo_db_string)
        
