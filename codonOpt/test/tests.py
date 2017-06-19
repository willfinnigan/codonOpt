import unittest
import codonOpt.Sequence_maker
from testfixtures import ShouldRaise




class Sequence_maker_checker_functions_tests(unittest.TestCase):

    def setUp(self):
        self.correct_protein_seq = "vgfENPqmkTSWardLYHic*"
        self.incorrect_protein_seq = "vgfENPqmkTSWardLYHic*123]["
        self.example_dna_seq = 'AGTTAGGCCTGCCTTATATTACCAAGGGCACAGTGAGGTAACCCCCCGGTAAAGTCGTTCAGACACACATAAGTCCATGAGGCGATTGTTGAACGATTGGATGTGGACTGTACGGCTCCTTTTAGCTCTCTATCACAAGGAGGCATGACCCGTCTCAAACGGAATACTCTGTGGTATTACCGCTCCGGGATC'
        self.example_protein_seq = 'S*ACLILPRAQ*GNPPVKSFRHT*VHEAIVERLDVDCTAPFSSLSQGGMTRLKRNTLWYYRSGI'
        self.example_protein_seq_different = 'PRAQ*GNPPKSFRHT*VHEALDVDCTAPFSSLSQGGMLKRNTLWYYRSGI'

    def test_check_protein_seq(self):
        ''' Test that check protein sequence function passes proteins containing correct amino acids '''

        self.assertTrue(codonOpt.Sequence_maker.check_protein_seq(self.correct_protein_seq))

    def test_check_protein_seq_wrong_aa(self):
        ''' Test that check protein sequence function fails proteins containing incorrect amino acids '''
        with ShouldRaise(NameError('Non amino acid character found in protein sequence')):
            self.assertFalse(codonOpt.Sequence_maker.check_protein_seq(self.incorrect_protein_seq))

    def test_translate_can_translate_correctly(self):
        translated = codonOpt.Sequence_maker.translate(self.example_dna_seq)

        self.assertEqual(translated,self.example_protein_seq)

    def test_that_check_dna_back_translation_picks_up_error(self):

        with ShouldRaise(NameError('Translation of dna seq to protein seq does not give the original protein seq')):
            self.assertFalse(codonOpt.Sequence_maker.check_dna_back_translation(self.example_dna_seq, self.example_protein_seq_different))

    def test_that_check_dna_back_translation_passes_correct_seq(self):

        self.assertTrue(codonOpt.Sequence_maker.check_dna_back_translation(self.example_dna_seq,
                                                                                self.example_protein_seq))


class SequenceMaker_TestCase(unittest.TestCase):

    def setUp(self):
        self.codon_table = codonOpt.Sequence_maker.CodonTable(json_file_path='example_data/Tth codon table.json', low_cuttoff=0.1)
        self.protein_seq = "MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**"
        self.dna_seq = codonOpt.Sequence_maker.generate_seq(self.protein_seq, self.codon_table)

    def test_that_dna_seq_made_is_three_times_protein_seq(self):
        """Test DNA seq is 3* length of protein seq"""

        test=False
        if len(self.dna_seq)/len(self.protein_seq) == 3:
            test = True

        self.assertTrue(test)

    def test_dna_backtranslates_to_protein(self):
        translated = codonOpt.Sequence_maker.translate(self.dna_seq)
        self.assertEqual(translated, self.protein_seq)

    def test_that_codon_table_imports_correctly(self):
        self.assertTrue(self.codon_table.check_codon_table())






if __name__ == '__main__':
    unittest.main()