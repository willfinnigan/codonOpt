import unittest
import codonOpt.Sequence_maker



class SequenceMaker_Tests(unittest.TestCase):

    def setUp(self):
        self.codon_table = codonOpt.Sequence_maker.CodonTable(json_file_path='Tth codon table.json', low_cuttoff=0.1)
        self.protein_seq = "MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**"
        self.dna_seq = codonOpt.Sequence_maker.generate_seq(self.protein_seq, self.codon_table)

    def test_that_dna_seq_made_is_three_times_protein_seq(self):
        """ Test DNA seq is 3* length of protein seq """

        test=False
        if len(self.dna_seq)/len(self.protein_seq) == 3:
            test = True

        self.assertTrue(test)

    def test_dna_backtranslates_to_protein(self):
        translated = codonOpt.Sequence_maker.translate(self.dna_seq)
        self.assertEqual(translated, self.protein_seq)




if __name__ == '__main__':
    unittest.main()