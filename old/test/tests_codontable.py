import unittest
import codonOpt.SeqMake
from testfixtures import ShouldRaise
from codonOpt.global_vars import ROOT_DIR
import logging

class SeqMake_tests(unittest.TestCase):

    def setUp(self):
        self.json_path = 'Escherichia_coli_K12.json'

    def test_make_codon_table_class(self):
        pass






if __name__ == '__main__':
    unittest.main()