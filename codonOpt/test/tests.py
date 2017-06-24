import unittest
import codonOpt.SeqMake
from testfixtures import ShouldRaise
from lea import *
import logging


class SeqMake_Functions_tests(unittest.TestCase):

    def setUp(self):
        self.correct_protein_seq = "vgfENPqmkTSWardLYHic*"
        self.incorrect_protein_seq = "vgfENPqmkTSWardLYHic*123]["
        self.example_dna_seq = 'agtTAGGCCTGCCTTATATTACCAAGGGCACAGTGAGGTAACCCCCCGGTAAAGTCGTTCAGACACACATAAGTCCATGAGGCGATTGTTGAACGATTGGATGTGGACTGTACGGCTCCTTTTAGCTCTCTATCACAAGGAGGCATGACCCGTCTCAAACGGAATACTCTGTGGTATTACCGCTCCGGGATC'
        self.example_protein_seq = 's*ACLILPRAQ*GNPPVKSFRHT*VHEAIVERLDVDCTAPFSSLSQGGMTRLKRNTLWYYRSGI'
        self.example_protein_seq_different = 'pRAQ*GNPPKSFRHT*VHEALDVDCTAPFSSLSQGGMLKRNTLWYYRSGI'

        self.test_seq = codonOpt.SeqMake.Sequence()

        self.input_dir = 'example_data/'

    def test_check_protein_seq(self):
        ''' Test that check protein sequence function passes proteins containing correct amino acids '''

        self.assertTrue(codonOpt.SeqMake.check_protein_seq(self.correct_protein_seq))

    def test_check_protein_seq_wrong_aa(self):
        ''' Test that check protein sequence function fails proteins containing incorrect amino acids '''
        with ShouldRaise(NameError('Non amino acid character found in protein sequence')):
            self.assertFalse(codonOpt.SeqMake.check_protein_seq(self.incorrect_protein_seq))

    def test_translate_can_translate_correctly(self):
        translated = codonOpt.SeqMake.translate(self.example_dna_seq)

        self.assertEqual(translated,self.example_protein_seq.upper())

    def test_that_checks_dna_back_translation_picks_up_error(self):

        with ShouldRaise(NameError('Translation of dna seq to protein seq does not give the original protein seq')):
            self.assertFalse(codonOpt.SeqMake.check_dna_back_translation(self.example_dna_seq, self.example_protein_seq_different))

    def test_that_check_dna_back_translation_passes_correct_seq(self):

        self.assertTrue(codonOpt.SeqMake.check_dna_back_translation(self.example_dna_seq,
                                                                    self.example_protein_seq))

    def test_import_protein_seqs_returns_list_SequenceObjects(self):

        list_seqs = codonOpt.SeqMake.import_protein_seqs(self.input_dir, 'protein_seqs.fasta', 'fasta')

        passed = True

        for seq in list_seqs:
            if type(seq) != type(self.test_seq):
                passed = False

        self.assertTrue(passed, 'Function import_protein_seqs did not return a list of Sequence objects')

    def test_import_protein_seqs_gets_name_and_sequence(self):

        list_seqs = codonOpt.SeqMake.import_protein_seqs(self.input_dir, 'protein_seqs2.fasta', 'fasta')

        passed = True

        if list_seqs[0].seq_name != 'mpCAR':
            passed = False

        if list_seqs[0].protein_seq != 'MASESRDVRLQRRIAELYDTDPQFAAARP':
            passed = False

        self.assertTrue(passed, 'Function import_protein_seqs import names and/or sequences')

    def test_translate_if_dna_not_divisable_by_three(self):
        dna_not_by_three = 'ATGGGTATGT'
        with ShouldRaise(NameError('DNA Sequence to be translated is not divisable by 3')):
            self.assertFalse(codonOpt.SeqMake.translate(dna_not_by_three))

class SeqMake_test_case(unittest.TestCase):

    def setUp(self):
        self.ct_dir = 'example_data/'
        self.codon_table = codonOpt.SeqMake.CodonTable(codon_tables_dir=self.ct_dir, json_file='Tth codon table.json', low_cuttoff=0.1)
        self.protein_seq = "MravVFENKERVAVKevnaPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**"
        self.dna_seq = codonOpt.SeqMake.generate_seq(self.protein_seq, self.codon_table.lea_codon_dict)

        self.test_dict = {"TCT": 1500, "TCC": 2500, "AGT": 2200, "TCA": 1800, "TCG": 750, "AGC": 1250}
        lea_test_dict = Lea.fromValFreqsDict(self.test_dict)
        self.test_lea_dict = {'S' : lea_test_dict}
        self.test_lea_protein = 'S' * 100000

    def test_that_dna_seq_made_is_three_times_protein_seq(self):
        """Test DNA seq is 3* length of protein seq"""

        test=False
        if len(self.dna_seq)/len(self.protein_seq) == 3:
            test = True

        self.assertTrue(test)

    def test_dna_backtranslates_to_protein(self):
        translated = codonOpt.SeqMake.translate(self.dna_seq)
        self.assertEqual(translated, self.protein_seq.upper())

    def test_that_codon_table_imports_correctly(self):
        self.assertTrue(self.codon_table.check_codon_table())

    def test_lea_selects_weighted_triplets_in_generate_seq(self):
        dna_seq = codonOpt.SeqMake.generate_seq(self.test_lea_protein, self.test_lea_dict)

        count_dict = {"TCT": 0, "TCC": 0, "AGT": 0, "TCA": 0, "TCG": 0, "AGC": 0}

        for i in range(0,len(dna_seq),3):
            triplet = dna_seq[i:i+3]
            count_dict[triplet] += 1

        passed = True
        for triplet in count_dict:
            total = count_dict[triplet] / 10

            test_against = self.test_dict[triplet]

            if total*0.9 <= test_against <= total*1.1:
                logging.debug(str(triplet) + ' passed (within 10%). Counted= ' + str(total) + ' vs test= ' + str(test_against))
                passed = True
            else:
                passed = False

        self.assertTrue(passed)

    def test_optimisation_of_dnaseq(self):
        passed = False
        test_dna = 'ATGCTCGACGACGCGCGGGCCGAGCGGCGGGAGCGGCGCATCGCCGACGCCCTCGCCGACGACCAGGTGCGCGAGGCCGCGGCGGACGCCGCGGTGTCCGAGTCGGTGCGGCGGGTGGAGGTCCGGCTGGCGCGCATCGTGGACGCCGTGATGTCCGGCTACGGGGACAGGGCCGCCCTGGCGTGGCGCCGGAGCGAGCTGGTGGACGGGGCCGTCAGGCTCCTCCCGGAGTACTCCACCATGACGTACCGGGAACTCTGGCGGCAGGCCGGCGCGGTGGCCGCCGAGTGGGGGGCGGACGCCGAGAGCCCGGTGAGGGCCGAGGACTTCGTGTGCACCCTCGGCTTTACCAGCCCGGACTACACCGTGGTGGACCTCGCGCTCATGCGGCTTGCCGCCGTGGCCGTGCCCCTGCAGGCGAGCGCCTCCGTGGCGCAGTGGAGGTCCATCATGGCGGAGACCGAGCCCCGCATGCTGGCCGCCTCGGCCGAGACCCTTCCCGCCGCCGTGGAGGCCGTCCTCGGGGGCTTCGCCCCCCGGCGGGTGCTCGTGTTTGACTACCGGCCCGAGCTCGAAGCCCACCGCAGCGCCGTGGACTCCGCCCGCGAGAGGCTCGCCGAGGTGGGCTGCACGGTGGCCACCGTGGCCGACGCCGTCGACCGGGGGGCCAACCTCCCCGCCCCGCTTCGGATCCCCTCCGACCGCGAGCGGCTCGCCCTCCTCATCTACACGTCCGGGTCCACCGGGGCCCCCAAGGGGGCGATGTACACCGACCGGCTCGTGGCGGGCCTGTGGCTGAGCGCCAACGAGATCCGCGTCCCCGCGCTCACCATGAACTACATGCCGCTGAGCCACATCGCCGGCAGGATGTCCCTTTACGGGACCCTGATGCGGGGGGGCACCGCGTACTTCGCGGCCGCCTCCGACATGAGCACGCTCCTCGACGACTTCGGGCTCGCCCGCCCGACCGAACTCTTCCTTGTCCCCCGGGTGTGCGAGCTGCTGCACCAACGGTACCAAAGCGAGCTCGACCGCCGCGTCGTGGCGGGGGAGGACGCCGAGACGGCCGCCACCAACGTGAAGGCCGAGCTCCGGGAGCGGGTCCTTGGGGGGCGCTACCTTACCGCGCTCTCGGGGAGCGCCCCCCTCGCCGCGGAAATGAAGACCTTCATGGAGTCCCTGCTGGACGACGAGCTCCACGACGGCTACGGGTCCACGGAGGCCGGCGGCAGCGTGCTCCTGGACAACCGCATCAAGCGGCCCCCCGTGCTTGACTACCGCCTGGTCGACGTCCCCGAGCTGGGGTACTTCAGGACCGACAAGCCCCACCCCCGCGGGGAGCTCCTCCTGAAGACCGAGTCCATGTTCCCCGGGTACTACAAGCGCCCCGAGATCACCGCGGAAATGTTTGACGCCGACGGGTTTTACCGGACCGGGGACGTGGTCGCCGAGCTCGGCCCCGAGCAGCTCGTCTACGTCGACAGGAGGAACAACGTGCTCAAGCTGTCCCAGGGGGAGTTTGTCACCGTGGCCGCCCTTGAGGCCGTGTACGCCACCTCGCCGCTGATCCGCCAGATCTTCGTGTACGGGTCCTCCGAGCGGGCGTACCTCCTCGCCGTGGTGGTCCCCACCGACGCCGTGCTTGCGCTCCCCGCCGCGCGCGCCCGCGCGGAGGTGTCCGAGTCCCTTCAGCGGATCGCCAAGGAGTCCGGCCTCCGCCCCTACGAAATCCCCAGGGACCTCATCATCGAGAGCGAGCCGTTTACGATCGACAACGGGCTCCTCTCCGGGATCGGGAAGCTGCTCCGGCCGAAGCTCAAGGAGCACTACGGCGAGCGGCTCGAACAGCTCTACGCCGAACTCGCCGAGCAGCGGGAGGACGAGCTCACCGCGCTGAGGCGCGGGGCCCACGACAGGCCCATCCTCGACACGGTGACGCGCGCCGCCGGGGCCGTGCTGGACCTGACCGCGGGCGAGGTGAGCCCCGACGCCCACTTCGCCGACCTGGGGGGCGACTCCCTCTCCGCCCTCTCCTTCAGCACCCTCCTCCGGGACATCTTCGGGGTGGAGGTCCCCGTGGGGTTCATCGTGGGGCCCGCGACCGACCTCGCCAGGATCGCCGAGTACCTCGTCAGCGAGCGCGACAGCGGGAGCCGGCCGACCGCCGCCACGGTGCACGGGGACGACGGCCTCCTCAGGGCCGACGACCTCGCCCTCGAGGCGTTCCTGGACCCCGCCACCCTGGACGCCGCCGCCCACCTCCCCAGCGCCCTGGAGCCCCCCCGCACCGTGCTGCTCACCGGCGCCAACGGCTACCTCGGCCGGTTCCTCGCCCTGGAGTGGCTGCAGCGCCTCGACGTGTCGGGGGGGACCCTGATCTGCCTGATCCGGGGGAGCGACGCCGACTCCGCCCGCCGGCGGCTCGACGCCGTGTTTGCCACCGGGGACCCCGAACTGGAGGCCCACTACCGCGAGCTTGCCGAGAGGCGGCTGCGCGTGCTGCCCGGCGACATCGGGGAGCCCAACCTGGGGCTGCGCGAGCAGGACTGGAGGGACCTCGCGGAGACGGTGGACCTGATCGTCCACCCCGCCGCCCTCGTCAACCACGTCCTGCCCTACGCCCAGCTCTTCGGGCCCAACGTGGTCGGCACCGCCGAGGTCATCAGGCTGGCCCTCACCAGCCGGCTGAAGCCCGTGACCTACCTCAGCACGGTGGCGGTCAGCGCCGGGATCGACCCCGAGACCTTTACCGAGGACGGGGACATCCGGGAGATCAGCCCCGTCCGCCGCCTCGACGACGGGTACGCCAACGGGTACGGCAACAGCAAGTGGGCGGGGGAGGTGCTTCTGCGCAACGCGCACGACCGGTTCGGCCTGCCCGTGGCCGTCTTCCGCTCCGACATGATCCTCGCCCACAGCCGCTACGCCGGCCAACTCAACGTGCCCGACATGTTCACCCGCCTGCTGCTCTCCGTGCTCGCCACCGGCCTCGCCCCCGGGTCCTTCCACGACGCGCACGGCGAGCGGCACCGGGCCCACTACGACGGGCTCCCCGCCGACTTCACCGCCGCCGCGGTCACCACCCTGGGGTCCCGGGTCACCTCCGGGTACGAGACCTACGACGTGCTGAACCCCCACGACGACGGCATCAGCCTCGACACCTTTGTGGACTGGCTTATCGAGGCGGGGCACCCCATCGACAGGATCGACGACTACGCCGAGTGGTTCGCCAGGTTCGACACCGCCCTCCGCGCCCTCCCCGAACACCAGCGGCAGCACTCCCTGCTCCCCCTGCTCCACGCGTACCGGAGGCCCACCCCCCCCCTTCACGGCGTCGCGCTCCCCGCCAAGCACTTCAGGGCCGCCGTCCAGCAGGCCAAGCTCGGCCCCGACGGCGACATCCCCCACGTGACCAGGGAGCTCATCGAAAAGTACGCCTCGGACCTTAGGCTCCTCGGCCTGATCCAGGGGTGA'
        test_protein = 'MLDDARAERRERRIADALADDQVREAAADAAVSESVRRVEVRLARIVDAVMSGYGDRAALAWRRSELVDGAVRLLPEYSTMTYRELWRQAGAVAAEWGADAESPVRAEDFVCTLGFTSPDYTVVDLALMRLAAVAVPLQASASVAQWRSIMAETEPRMLAASAETLPAAVEAVLGGFAPRRVLVFDYRPELEAHRSAVDSARERLAEVGCTVATVADAVDRGANLPAPLRIPSDRERLALLIYTSGSTGAPKGAMYTDRLVAGLWLSANEIRVPALTMNYMPLSHIAGRMSLYGTLMRGGTAYFAAASDMSTLLDDFGLARPTELFLVPRVCELLHQRYQSELDRRVVAGEDAETAATNVKAELRERVLGGRYLTALSGSAPLAAEMKTFMESLLDDELHDGYGSTEAGGSVLLDNRIKRPPVLDYRLVDVPELGYFRTDKPHPRGELLLKTESMFPGYYKRPEITAEMFDADGFYRTGDVVAELGPEQLVYVDRRNNVLKLSQGEFVTVAALEAVYATSPLIRQIFVYGSSERAYLLAVVVPTDAVLALPAARARAEVSESLQRIAKESGLRPYEIPRDLIIESEPFTIDNGLLSGIGKLLRPKLKEHYGERLEQLYAELAEQREDELTALRRGAHDRPILDTVTRAAGAVLDLTAGEVSPDAHFADLGGDSLSALSFSTLLRDIFGVEVPVGFIVGPATDLARIAEYLVSERDSGSRPTAATVHGDDGLLRADDLALEAFLDPATLDAAAHLPSALEPPRTVLLTGANGYLGRFLALEWLQRLDVSGGTLICLIRGSDADSARRRLDAVFATGDPELEAHYRELAERRLRVLPGDIGEPNLGLREQDWRDLAETVDLIVHPAALVNHVLPYAQLFGPNVVGTAEVIRLALTSRLKPVTYLSTVAVSAGIDPETFTEDGDIREISPVRRLDDGYANGYGNSKWAGEVLLRNAHDRFGLPVAVFRSDMILAHSRYAGQLNVPDMFTRLLLSVLATGLAPGSFHDAHGERHRAHYDGLPADFTAAAVTTLGSRVTSGYETYDVLNPHDDGISLDTFVDWLIEAGHPIDRIDDYAEWFARFDTALRALPEHQRQHSLLPLLHAYRRPTPPLHGVALPAKHFRAAVQQAKLGPDGDIPHVTRELIEKYASDLRLLGLIQG*'

        translated = codonOpt.SeqMake.translate(test_dna)
        optimised_dna = codonOpt.SeqMake.generate_seq(translated, self.codon_table.lea_codon_dict)

        if optimised_dna != test_dna:
            passed = True

        translated_two = codonOpt.SeqMake.translate(optimised_dna)
        if translated_two == test_protein:
            passed = True

        self.assertTrue(passed)

class SeqMake_multipleSeqObjs_testcase(unittest.TestCase):

    def setUp(self):
        self.input_dir = 'example_data/'
        self.ct_dir = 'example_data/'
        self.list_seqs = codonOpt.SeqMake.import_protein_seqs(self.input_dir, 'protein_seqs.fasta', 'fasta')
        self.codon_table = codonOpt.SeqMake.CodonTable(codon_tables_dir=self.ct_dir, json_file='Tth codon table.json', low_cuttoff=0.1)
        self.dna_seqs = codonOpt.SeqMake.codon_optimise_seq_list(self.list_seqs, self.codon_table.lea_codon_dict)

    def test_codon_optimise_seq_list_dnalen(self):
        self.assertEqual(len(self.dna_seqs), 5)

    def test_list_seqs_contains_dna_seqs_that_translate_to_protein(self):

        passed = True
        for seq in self.list_seqs:
            translated = codonOpt.SeqMake.translate(seq.dna_seq)

            if translated != seq.protein_seq.upper():
                passed = False

        self.assertTrue(passed, 'Imported list of sequences which has then been codon optimised cannot be back translated from dna')

    def test_dna_output_as_fasta(self):
        pass
        # TODO finish this test










if __name__ == '__main__':
    unittest.main()