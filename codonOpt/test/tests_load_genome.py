import unittest
import codonOpt.SeqMake, codonOpt.GenomeAnalysis.Analysis
from codonOpt.GenomeAnalysis.Analysis import Genome, CDS
from testfixtures import ShouldRaise
from lea import *
import logging




class GenomeAnalysisCreation_tests(unittest.TestCase):

    def setUp(self):
        example_data_dir = 'example_data/'
        self.name = 'Escherichia_coli_K12'
        self.genome_dir = example_data_dir + 'GenomeTmp/'
        self.codon_table_dir = example_data_dir + 'Codon_tables/'
        self.pickle_dir = example_data_dir + "PickledGenomes/"

    def test_load_genome(self):
        genome = Genome(self.name)
        genome.load_all_genbanks_in_folder(self.genome_dir)

    def test_load_pickle(self):
        genome = codonOpt.GenomeAnalysis.Analysis.load_pickle(self.name, self.pickle_dir)
        self.assertEqual(genome.name, self.name)

    def test_save_pickle(self):
        save_name = 'test_save'
        genome = codonOpt.GenomeAnalysis.Analysis.load_pickle(self.name, self.pickle_dir)
        codonOpt.GenomeAnalysis.Analysis.save_pickle(genome, save_name, self.pickle_dir)

        genome_loaded = codonOpt.GenomeAnalysis.Analysis.load_pickle(save_name, self.pickle_dir)
        self.assertEqual(genome_loaded.name, self.name)





if __name__ == '__main__':
    unittest.main()