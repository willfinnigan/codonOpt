import unittest
import codonOpt.SeqMake, codonOpt.GenomeAnalysis.Analysis
from codonOpt.GenomeAnalysis.Analysis import Genome, CDS
from testfixtures import ShouldRaise
from lea import *
import logging



class GenomeAnalysis_tests(unittest.TestCase):

    def setUp(self):
        example_data_dir = 'example_data/'
        self.name = 'Thermus_thermophilus_pTT27'
        self.genome_dir = example_data_dir + 'GenomeTmp/'
        self.codon_table_dir = example_data_dir + 'Codon_tables/'
        self.pickle_dir = example_data_dir + "PickledGenomes/"
        self.genome = codonOpt.GenomeAnalysis.Analysis.load_pickle(self.name, self.pickle_dir)

    def test_save_json(self):
        self.genome.save_relative_codon_table_as_json(self.name, dir=self.codon_table_dir)

    def test_analyse_codon_usage(self):
        self.genome.analyse_codon_usage()

    def test_analyse_initiation_secondary_structure(self):
        self.genome.initiation_ss_regions = [(-10,10)]
        self.genome.analyse_initiation_secondary_structure()

    def test_analyse_secondary_structure_in_windows(self):
        self.genome.ss_window_size_and_move = [(100, 25)]
        self.genome.analyse_secondary_structures_in_windows()

    def test_analyse_gc_windows(self):
        self.genome.gc_window_sizes= [25]
        self.genome.calc_gc_windows()

    def test_calc_codon_usage_at_start_of_genes(self):
        self.genome.codons_at_start_of_genes_sizes = [30]
        self.genome.calc_codon_usage_at_start_of_genes()

    def test_codon_context(self):
        self.genome.analyse_codon_context()

        for triplet in self.genome.codon_context_tables:
            print(triplet)
            print(self.genome.codon_context_tables[triplet])



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    #unittest.main()

    suite = unittest.TestLoader().loadTestsFromTestCase(GenomeAnalysis_tests)
    unittest.TextTestRunner(verbosity=2).run(suite)