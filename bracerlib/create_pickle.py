#!/usr/bin/env python

import os
from tasks import TracerTask
import pickle
from bracerlib import io


class CreateCell(TracerTask):

    def __init__(self, **kwargs):
        resource_dir = None
        self.cell_name = kwargs.get('cell_name')
        self.fastq1 = kwargs.get('fastq1')
        self.fastq2 = kwargs.get('fastq2')
        self.ncores = 1
        self.assembled_file = kwargs.get('assembled_file')
        self.species = 'Hsap'
        self.resume_with_existing_files = kwargs.get(
            'resume_with_existing_files')
        self.output_dir = kwargs.get('output_dir')
        self.single_end = False
        self.fragment_length = False
        self.fragment_sd = False
        self.loci = ['H', 'K', 'L']
        self.max_junc_len = 100
        self.no_trimming = True
        self.keep_trimmed_reads = None
        config_file = kwargs.get('config_file')

        self.trimmed_fastq1 = None
        self.trimmed_fastq2 = None
        self.config = self.read_config(config_file)
        self.species_root = self.get_species_root(self.species,
                                                  root=resource_dir)

    def create_cell(self):
        imgt_seq_location = os.path.join(self.species_root, 'raw_seqs')
        cell = io.parse_IgBLAST(self.loci, self.output_dir, self.cell_name,
                                imgt_seq_location, self.species, self.assembled_file,
                                self.max_junc_len)
        if cell.is_empty:
            self.die_with_empty_cell(self.cell_name, self.output_dir,
                                     self.species)

        ranked_recs = cell.rank_recombinants()
        isotype = cell.determine_isotype(ranked_recs)
        bgcolor = cell.assign_bgcolor(isotype)
        cell.bgcolor = bgcolor
        cell.isotype = isotype

        # Create database dictionary
        cell.databasedict = cell.get_database_for_locus(self.loci)
        return cell


def run():
    """
    This function reads creates the cell2 and cell3 pickle files for the test data. This is neccessary to be re-run after
    a change in the API, such as when Bio.Alphabet was removed from BioPython. After running you need to replace the files
    in test_data/results/cell{2,3}/{un}filtered_BCR_Seqs/cell{2,3}.pkl with the cell{2,3}.pkl files produced by this.
    """
    print("Regenerating the cell2 and cell3 pickle files for test_Data")
    cell2_pickle = "cell2.pkl"
    cell3_pickle = "cell3.pkl"
    cell_creator = CreateCell(output_dir='test_data/results/cell2', cell_name='cell2')
    cell2 = cell_creator.create_cell()
    print(cell2.__dict__["recombinants"]["BCR"]["H"][0].__dict__)
    with open(cell2_pickle, 'wb') as pf:
        pickle.dump(cell2, pf, protocol=5)
    cell_creator = CreateCell(output_dir='test_data/results/cell3', cell_name='cell3')
    cell3 = cell_creator.create_cell()
    with open(cell3_pickle, 'wb') as pf:
        pickle.dump(cell3, pf, protocol=5)


if __name__ == '__main__':
    run()
