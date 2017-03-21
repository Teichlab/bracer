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
        self.receptor = receptor
        self.recombinants = self._process_recombinants(recombinants, receptor, loci)
               
        self.is_empty = self._check_is_empty()
        self.two_most_common = self.find_two_most_common()
        self.ranked_recs = self.rank_recombinants()
        self.bgcolor = None
        self.isotype = None
        #self.isotype = self.determine_isotype(loci, receptor, self.recombinants)
        #self.bgcolor = self.assign_bgcolor(species, self.isotype)
        #self.changeodict = self.get_changeo_db_for_locus(self.recombinants, receptor, loci)
        #self.print_dict = self.print_recombinants()
        self.detailed_identifier_dict = self.create_detailed_identifier_dict()
        self.cdr3_dict = self.find_recs_with_identical_cdr3()
        
        #self.rank_recs = self.rank_recombinants()
        #self.two_most_common = self.find_n_most_common(2)
        #self.identical = self.assert_two_most_common()
        
        #self.three_most_common = self.find_n_most_common(3)
        #self.four_most_common = self.find_n_most_common(4)
        #self.replacement_dict = self.assert_third_most_common()
        #self.detailed_identifier_dict = self.create_detailed_identifier_dict()
        #self.tpm_dict = self.create_tpm_dict(self.recombinants, loci)
        #print(self.cdr3_dict)
        #print(self.rank_recs)

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
        #ranked_recs = self.ranked_recs
        ranked_H = ranked_recs["H"][0:2]
        H_recombinants = self.recombinants[self.receptor]["H"]
        isotype = None
        isotype_list = []
        for ranked_rec in ranked_H:
            for recombinant in H_recombinants:
                if recombinant.contig_name == ranked_rec and recombinant.productive:
                    C_gene = recombinant.C_gene
                    if C_gene is not None:
                        C_gene = C_gene.split("*")[0]
                        isotype_list.append(C_gene)
        if len(isotype_list) == 0:
            isotype = None
        elif len(isotype_list) == 1:
            isotype = isotype_list[0]
        elif len(isotype_list) == 2:
            iso_1 = isotype_list[0]
            iso_2 = isotype_list[1]
            if iso_1 == iso_2:
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
            isotype_bgcolors = {"IGHD":'#e6f7ff', "IGHM":'#e5ffcc', "IGHA":'#ffe6e6', 
                "IGHE":'#ffffcc', "IGHG1":'#f1e6ff', "IGHG2A":'#e2ccff', "IGHG2B":'#d4b3ff', "IGHG2C":'#c599ff',
                "IGHG3":'a866ff', "IGHDM":'#b3ffff'}
        else:
            # Isotypoe background colours for Hsap
            isotype_bgcolors = {"IGHD":'#e6f7ff', "IGHM":'#e5ffcc', "IGHA1":'#ffe6e6', "IGHA2":'#ffcccc', 
                "IGHE":'#ffffcc', "IGHG1":'#f1e6ff', "IGHG2":'#e2ccff', "IGHG3":'#d4b3ff', "IGHG4":'#c599ff', "IGHDM":'#b3ffff'}

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
        
        recombinants = dict()
        final_string = '<<table cellspacing="6px" border="0" cellborder="0">'
        # final_string = "<"
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
                                                                            receptor=receptor, locus=locus,
                                                                            identifier=rec.identifier)
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
                

    #Code to detect recombinants with same VDJ rearrangement but different isotypes (naive B cells)

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
        """Ranks recombinants from highest TPM to lowest in case of more than two recombinants for a locus"""
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


    """def find_n_most_common(self, n):
        most_common_dict = dict()
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None:
                    most_common_dict[locus] = []
                    most_common = []
                    if len(recombinants) > 2:
                        TPM_ranks = Counter()
                        for rec in recombinants:
                            TPM_ranks.update({rec.contig_name: rec.TPM})

                        most_com = [x[0] for x in TPM_ranks.most_common(n)]
                        #most_common = []
                        for rec in recombinants:
                            if rec.contig_name in most_com:
                                name = rec.contig_name
                                most_common.append(name)
                        most_common_dict[locus] = most_common
                    
                    elif len(recombinants) >= 1:
                        for rec in recombinants:
                            name = rec.contig_name
                            most_common.append(name)
                            
                        most_common_dict[locus] = most_common

        return (most_common_dict)"""


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


    def assert_two_most_common(self, ranked_recs):
        pass      


    """def assert_two_most_common(self):
        two_identical = dict()
        keep_one_rec = dict()
        detailed_identifier_dict = self.create_detailed_identifier_dict()
        two_most_common_dict = self.find_n_most_common(2)
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None:
                    two_identical[locus] = False
                    keep_one_rec[locus] = False
                    if len(recombinants) > 1:
                        two_most_common = two_most_common_dict[locus]
                        #print("PRINTING TWO MOST COMMON DICT")
                        #print(locus)
                        #print(two_most_common)
                        
                        rec1 = None
                        rec2 = None
                        #print("Should be None")
                        #print(rec1)
                        #print(rec2)
                        
                        for rec in recombinants:
                            if rec.contig_name in two_most_common:
                                if rec1 == None:
                                    rec1 = rec.contig_name
                                elif rec2 == None:
                                    rec2 = rec.contig_name
                            elif rec in two_most_common:
                                if rec1 == None:
                                    rec1 = rec.contig_name
                                elif rec2 == None:
                                    rec2 = rec.contig_name

                        #for i in range(len(two_most_common)):
                            #for rec in recombinants:
                                #if rec.contig_name == two_most_common[i]:
                                    #two_most_common[i] = rec
                        #rec1 = two_most_common[0]
                        #rec2 = two_most_common[1]
                        #print("Should be two most common recs for locus")
                        #print(rec1)
                        #print(rec2)
                        if detailed_identifier_dict[locus][rec1].split("*")[0] == detailed_identifier_dict[locus][rec2].split("*")[0]:
                            two_identical[locus] = True
                            isotype_rec1 = detailed_identifier_dict[locus][rec1].split("*")[1].split("_")[0]
                            isotype_rec2 = detailed_identifier_dict[locus][rec2].split("*")[1].split("_")[0]
                            full_length_rec1 = detailed_identifier_dict[locus][rec1].split("*")[1].split("_")[1]
                            full_length_rec2 = detailed_identifier_dict[locus][rec2].split("*")[1].split("_")[1]
                            if detailed_identifier_dict[locus][rec1].split("*")[1] == detailed_identifier_dict[locus][rec2].split("*")[1]:
                                keep_rec = rec1
                            elif full_length_rec1 == "true":
                                keep_rec = rec1
                            elif full_length_rec2 == "true":
                                keep_rec = rec2
                            if isotype_rec1 != isotype_rec2:
                                if isotype_rec1 == "none":
                                    keep_isotype = isotype_rec2
                                elif isotype_rec2 == "none":
                                    keep_isotype = isotype_rec1
                                elif (isotype_rec1 and isotype_rec2) in ["IGHD", "IGHM"]:
                                    keep_isotype = "IGHDM"
                            keep_one_rec[locus] = (keep_rec, keep_isotype)
        #print("KEEP ONE REC")
        #print(keep_one_rec)
        return (keep_one_rec)""" 
 

    """def assert_two_most_common(self):
        two_identical = dict()
        most_common_dict = self.two_most_common
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None:
                    two_identical [locus] = False
                    if len(recombinants) > 1:
                        most_common = most_common_dict[locus]
                        for i in range(len(most_common)):
                            for rec in recombinants:
                                if rec.contig_name == most_common[i]:
                                    most_common[i] = rec
                        rec1 = most_common[0]
                        rec2 = most_common[1]
                        if rec1.cdr3 == rec2.cdr3:
                            if (rec1.productive and rec2.productive) or not (rec1.productive and rec2.productive):
                                two_identical[locus] = True
                        
        return (two_identical)"""

    def assert_third_most_common(self):
        three_most_common_dict = self.find_n_most_common(3)
        replacement_rec_dict = dict()
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                replacement_rec_dict[locus] = None
                if self.identical[locus] is not False:
                    if len(recombinants) > 2:
                        three_most_common = three_most_common_dict[locus]
                        for i in range(len(three_most_common)):
                            for rec in recombinants:
                                if rec.contig_name == three_most_common[i]:
                                    three_most_common[i] = rec
                                    if rec.contig_name in self.two_most_common[locus]:
                                        comp_rec = rec
                                    else:
                                        third_rec = rec
                        if comp_rec.cdr3 == third_rec.cdr3:
                            replacement_rec = False
                            
                        else:
                            replacement_rec = third_rec
                        replacement_rec_dict[locus] = replacement_rec
        return (replacement_rec_dict)
                            

    """def assert_two_most_common(self):
        two_identical = dict()
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None:
                    two_identical[locus] = False
                    if len(recombinants) > 1:
                        TPM_ranks = Counter()
                        for rec in recombinants:
                            TPM_ranks.update({rec.contig_name: rec.TPM})

                        most_common = [x[0] for x in TPM_ranks.most_common(2)]
                        for i in range(len(most_common)):
                            for rec in recombinants:
                                if rec.contig_name == most_common[i]:
                                    most_common[i] = rec
                        for rec in most_common:
                            print(rec.contig_name)
                            print(rec.TPM)
                            print(rec.cdr3)
                        rec1 = most_common[0]
                        rec2 = most_common[1]
                        if rec1.cdr3 == rec2.cdr3:
                            if (rec1.productive and rec2.productive) or not (rec1.productive and rec2.productive):
                                two_identical[locus] = True

        return (two_identical)"""


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
                        print_dict[locus][rec] = [rec.contig_name, tpm, rec.cdr3, C_gene, rec.productive]

        return (print_dict)
                           
    def create_detailed_identifier_dict(self):
        #detailed identifier dict
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



    def find_two_most_expressed_unique(self):

        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None:
                    if len(recombinants) > 2:
                        TPM_ranks = Counter()
                        for rec in recombinants:
                            TPM_ranks.update({rec.contig_name: rec.TPM})
                        if receptor is "BCR":
                            most_common = [x[0] for x in TPM_ranks.most_common()]
                        for ranked_rec in most_common:
                            for rec in recombinants:
                                if rec.contig_name == ranked_rec:
                                    ranked_rec = rec
                        current_cells = most_common
                        comparison_cells = most_common[1:]
                        for current_cell in current_cells:
                            for comparison_cell in comparison_cells:
                                if current_cell.cdr3 == comparison_cell.cdr3:
                                    pass


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

                        """if receptor == "BCR":
                            two_most_common = [x[0] for x in TPM_ranks.most_common(2)]
                            three_most_common = [x[0] for x in TPM_ranks.most_common(3)]
                            rec1 = two_most_common[0][0]
                            rec2 = two_most_common[1][0]
                            to_remove = []
                        
                                    
                            
                            for rec in recombinants:
                                rec.isotype = rec.isotype.split("*")[0]
                                rec.J_gene = rec.J_gene.split("*")[0]
                                if rec.C_gene:
                                    rec.isotype = rec.C_gene.split("*")[0]
                            keep_third = False
                            naive = False
                            for curr_rec in recombinants:
                                if curr_rec.contig_name == rec1:
                                    rec1 = curr_rec
                                    for comp_rec in recombinants:
                                        if comp_rec.contig_name == rec2:
                                            rec2 = comp_rec
                                            if (curr_rec.productive and comp_rec.productive) or not (curr_rec.productive and comp_rec.productive):
                                                if curr_rec.cdr3 == comp_rec.cdr3:
                                                    keep_third = True
                                                    
                                                    if curr_rec.tpm > comp_rec.tpm:
                                                        to_remove.append(comp_rec)
                                                    else:
                                                        to_remove.append(curr_rec)
                                                    if (curr_rec.isotype != comp_rec.isotype) and (curr_rec.isotype and comp_rec.isotype) in ["IGHM", "IGHD"]:
                                                        if curr_rec.tpm > comp_rec.tpm:
                                                            curr_rec.isotype = "IGHDM"
                                                        else:
                                                            comp_rec.isotype = "IGHDM"
                                                        naive = True 
                    
                                for comparison_rec.contig_name in two_most_common:
                                for current_rec in recombinants:
                                
                                    if current_rec.contig_name == comparison_rec.contig_name:
                                        continue
                                    elif (current_rec.contig_name and comparison_rec.contig_name) in two_most_common and keep_third == False:
                                        # Comparing the two most highly expressed recs with each other
                                        if len(current_rec.V_genes.intersection(comparison_rec.V_genes)) > 0:
                                            if current_rec.J_gene == comparison_rec.J_gene:
                                                if (current_rec.isotype and comparison_rec.isotype) in ["IGHM", "IGHD"]:
                                                    if current_rec.isotype != comparison_rec.isotype:
                                                        if (current_rec.productive and comparison_rec.productive) or not (current_rec.productive and comparison_rec.productive):
                                                            if current_rec.TPM > comparison_rec.TPM:
                                                                to_remove.append(comparison_rec)
                                                            
                                                            else:
                                                                to_remove.append(current_rec)
                                                             
                                                            keep_third = True
                                                            
                                    elif current_rec not in two_most_common:
           
                                        if len(current_rec.V_genes.intersection(comparison_rec.V_genes)) > 0:
                                            if current_rec.J_gene == comparison_rec.J_gene:
                                                if (current_rec.isotype and comparison_rec.isotype) in ["IGHM", "IGHD"]:
                                                    if current_rec.isotype != comparison_rec.isotype:
                                                        if (current_rec.productive and comparison_rec.productive) or not (current_rec.productive and comparison_rec.productive):
                                                           comparison_rec.isotype = "IGHDM
                            
                            for rec in recombinants:
                                if keep_third == True:
                                    if rec.contig_name not in three_most_common:
                                        to_remove.append(rec)
                                else:
                                    if rec.contig_name not in two_most_common:
                                        to_remove.append(rec)
                            for rec in to_remove:
                                self.recombinants[receptor][locus].remove(rec)
                            for rec in recombinants:
                                if keep_third == True:
                                    if rec in two_most_common:
                                        rec.isotype = "IGHDM" """
                                   

 
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
                 imgt_reconstructed_seq, has_D, output_dir, full_length, query_length, V_genes, cdr3, assembler, C_gene, C_info_line, cdr3_seq):
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
        self.cdr3_seq = cdr3_seq
        self.detailed_identifier = self.create_detailed_identifier(self.productive, self.cdr3, self.C_gene, self.full_length)        

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
        summary_string += "TPM:\t{TPM}\nProductive:\t{productive}\nStop codon:\t{stop_codon}\nIn frame:\t{in_frame}\nFull length:\t{full_length}\nSequence length: \t{query_length}\
                           \nPossible V genes:\t{V_genes}\n\n".format(TPM=self.TPM, productive=self.productive, stop_codon=self.stop_codon, 
                           in_frame=self.in_frame, full_length=self.full_length, query_length=self.query_length, V_genes=self.V_genes)

        summary_string += 'Segment\tquery_id\tsubject_id\t% identity\talignment length\tmismatches\tgap opens\tgaps\tq start\tq end\ts start\ts end\te value\tbit score\n'

        for line in self.hit_table:
            summary_string = summary_string + "\t".join(line) + "\n"
        
        if find_C_gene == True and C_gene != None:
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
                
                
        
        
