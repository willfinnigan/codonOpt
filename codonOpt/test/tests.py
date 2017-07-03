import unittest
import codonOpt.SeqMake
from testfixtures import ShouldRaise
from lea import *
import logging


class SeqMake_tests(unittest.TestCase):

    def setUp(self):
        self.ct_dir = 'test/data/'

    def test_check_protein_seq(self):
        ''' Test that check protein sequence function passes proteins containing correct amino acids '''

        correct_protein_seq = "vgfENPqmkTSWardLYHic*"
        self.assertTrue(codonOpt.SeqMake.check_protein_seq(correct_protein_seq))

    def test_check_protein_seq_wrong_aa(self):
        ''' Test that check protein sequence function fails proteins containing incorrect amino acids '''

        incorrect_protein_seq = "vgfENPqmkTSWardLYHic*123]["

        with ShouldRaise(NameError('Non amino acid character found in protein sequence')):
            self.assertFalse(codonOpt.SeqMake.check_protein_seq(incorrect_protein_seq))

    def test_translate_can_translate_correctly(self):
        dna_seq = 'agtTAGGCCTGCCTTATATTACCAAGGGCACAGTGAGGTAACCCCCCGGTAAAGTCGTTCAGACACACATAAGTCCATGAGGCGATTGTTGAACGATTGGATGTGGACTGTACGGCTCCTTTTAGCTCTCTATCACAAGGAGGCATGACCCGTCTCAAACGGAATACTCTGTGGTATTACCGCTCCGGGATC'
        protein_seq = 's*ACLILPRAQ*GNPPVKSFRHT*VHEAIVERLDVDCTAPFSSLSQGGMTRLKRNTLWYYRSGI'

        translated = codonOpt.SeqMake.translate(dna_seq)

        self.assertEqual(translated,protein_seq.upper())

    def test_that_checks_dna_back_translation_picks_up_error(self):
        dna_seq = 'agtTAGGCCTGCCTTATATTACCAAGGGCACAGTGAGGTAACCCCCCGGTAAAGTCGTTCAGACACACATAAGTCCATGAGGCGATTGTTGAACGATTGGATGTGGACTGTACGGCTCCTTTTAGCTCTCTATCACAAGGAGGCATGACCCGTCTCAAACGGAATACTCTGTGGTATTACCGCTCCGGGATC'
        protein_seq_different = 'pRAQ*GNPPKSFRHT*VHEALDVDCTAPFSSLSQGGMLKRNTLWYYRSGI'

        with ShouldRaise(NameError('Translation of dna seq to protein seq does not give the original protein seq')):
            self.assertFalse(codonOpt.SeqMake.check_dna_back_translation(dna_seq, protein_seq_different))

    def test_that_check_dna_back_translation_passes_correct_seq(self):
        dna_seq = 'agtTAGGCCTGCCTTATATTACCAAGGGCACAGTGAGGTAACCCCCCGGTAAAGTCGTTCAGACACACATAAGTCCATGAGGCGATTGTTGAACGATTGGATGTGGACTGTACGGCTCCTTTTAGCTCTCTATCACAAGGAGGCATGACCCGTCTCAAACGGAATACTCTGTGGTATTACCGCTCCGGGATC'
        protein_seq = 's*ACLILPRAQ*GNPPVKSFRHT*VHEAIVERLDVDCTAPFSSLSQGGMTRLKRNTLWYYRSGI'

        self.assertTrue(codonOpt.SeqMake.check_dna_back_translation(dna_seq, protein_seq))

    def test_translate_if_dna_not_divisable_by_three(self):
        dna_not_by_three = 'ATGGGTATGT'
        with ShouldRaise(NameError('DNA Sequence to be translated is not divisable by 3')):
            self.assertFalse(codonOpt.SeqMake.translate(dna_not_by_three))

    def test_that_optimised_dna_seq_is_three_times_protein_seq(self):
        """Test DNA seq is 3* length of protein seq"""

        codon_table = codonOpt.SeqMake.CodonTable(codon_tables_dir=self.ct_dir, json_file='Tth_codon_table.json',
                                                       low_cuttoff=0.1)

        protein_seq = "MravVFENKERVAVKevnaPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**"
        dna_seq = codonOpt.SeqMake.generate_seq(protein_seq, codon_table.lea_codon_dict)

        test=False
        if len(dna_seq)/len(protein_seq) == 3:
            test = True

        self.assertTrue(test)

    def test_that_codon_table_imports_correctly(self):

        codon_table = codonOpt.SeqMake.CodonTable(codon_tables_dir=self.ct_dir, json_file='Tth_codon_table.json',
                                                       low_cuttoff=0.1)

        self.assertTrue(codon_table.check_codon_table())

    def test_lea_selects_weighted_triplets_in_generate_seq(self):
        test_dict = {"TCT": 1500, "TCC": 2500, "AGT": 2200, "TCA": 1800, "TCG": 750, "AGC": 1250}
        lea_test_dict = Lea.fromValFreqsDict(test_dict)
        test_lea_dict = {'S': lea_test_dict}
        test_lea_protein = 'S' * 100000

        dna_seq = codonOpt.SeqMake.generate_seq(test_lea_protein, test_lea_dict)

        count_dict = {"TCT": 0, "TCC": 0, "AGT": 0, "TCA": 0, "TCG": 0, "AGC": 0}

        for i in range(0,len(dna_seq),3):
            triplet = dna_seq[i:i+3]
            count_dict[triplet] += 1

        passed = True
        for triplet in count_dict:
            total = count_dict[triplet] / 10

            test_against = test_dict[triplet]

            if total*0.9 <= test_against <= total*1.1:
                logging.debug(str(triplet) + ' passed (within 10%). Counted= ' + str(total) + ' vs test= ' + str(test_against))
                passed = True
            else:
                passed = False

        self.assertTrue(passed)










if __name__ == '__main__':
    unittest.main()