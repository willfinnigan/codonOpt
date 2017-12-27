from lea import *
import logging, random, json

from codonOpt.global_vars import empty_codon_table, ROOT_DIR
from codonOpt.AnalysisRedesignTools import MotifFinder, GCTools, MFE
from codonOpt.AnalysisRedesignTools.GeneralFunctions import return_window_in_frame, translate, check_protein_seq, check_dna_back_translation

import os

"""SeqMake Module.

Summary:


Example:

1.  Set up a codon table, and a sequence generator to use that codon table

    codon_table = CodonTable(codon_tables_dir=codon_tables_dir, json_file='Escherichia_coli_K12.json', low_cuttoff=0.1)
    seq_gen = Sequence_Generator(codon_table)
    
2.  Optimise your protein sequence    
    
    protein_seq = "MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSD*"
    dna_seq = seq_gen.optimise(protein_seq)
    
3.  Remove any restriction sites or other sequence motifs you don't want

    motifs_to_remove = ["ATAA", 'GGAGA', "TACAG"]
    dna_seq = seq_gen.motif_removal(dna_seq, motifs_to_remove)
    
4.  Run other types of redesign routines to minimise folding energy, lower gc content ect..

    dna_seq = seq_gen.remove_high_GC_windows(dna_seq, 100, 67)
    dna_seq = seq_gen.minimise_mfe_windows(dna_seq, -10, 30, 5)
    dna_seq = seq_gen.minimise_mfe_five_prime(dna_seq, -5)

"""


class CodonTable():
    """CodonTable class - codontable.lea_codon_dict for optimisation, loads from json file

    This is the CodonTable Class
    It has codon_dict, lea_codon_dict and low_cuttoff variables.

    When initialising it takes the following arguments:
        codon_tables_dir = "" - This is the directory to find the codon table json file.
        json_file = "" - This is the json file containing the codon usage information.
        low_cuttoff = 0 - This is the lowcuttoff which is used to remove codons with low frequency

    """

    def __init__(self, codon_tables_dir = '', json_file='', low_cuttoff=0.0):

        self.codon_dict = {}
        self.lea_codon_dict = {}
        self.low_cuttoff = low_cuttoff

        if json_file != '':
            self.load_codon_table_from_json(codon_tables_dir, json_file)
            self.load_lea_dict()
            self.check_codon_table()

    def load_codon_table_from_json(self, codon_tables_dir, json_file):

        logging.info('')
        logging.info('---Load Codon Table---')
        # Open to codon table json.
        with open(codon_tables_dir+json_file, 'r') as fp:
            self.codon_dict = json.load(fp)

        # Make nested dicts from json into dicts
        for aa in self.codon_dict:
            self.codon_dict[aa] = dict(self.codon_dict[aa])

        # Logging
        logging.info('Loaded codon dict from json at: ' + str(codon_tables_dir + json_file))
        logging.info(str(self.codon_dict))
        logging.info('')

    def load_lea_dict(self):

        logging.info('')
        logging.info('---Make Lea Codon Table for optimisation---')
        logging.info('Low cuttoff has been set at ' + str(self.low_cuttoff) + '  (codons used less than ' + str(self.low_cuttoff*100) + '% of the time)' )
        x_codon_dict = {}
        count = 0
        for aa in self.codon_dict:
            x_codon_dict[aa] = {}
            for triplet in self.codon_dict[aa]:
                if self.codon_dict[aa][triplet] >= self.low_cuttoff:
                    x_codon_dict[aa][triplet] = int((self.codon_dict[aa][triplet] * 10000))
                else:
                    logging.info('Amino Acid: ' + str(aa) + ', triplet: ' + str(triplet) + ' has freq of ' + str(self.codon_dict[aa][triplet]) + ' (less than ' + str(self.low_cuttoff) + '), so not included' )
                    count += 1

            if len(x_codon_dict[aa]) == 0:
                logging.warning('Amino Acid ' + str(aa) + ' has 0 triplets above low cuttoff')

        # Logging

        logging.info(str(count) + ' triplets lower than cuttoff ' + str(self.low_cuttoff) + ' not included')


        logging.debug('Codon dict multiplied by ' + str(10000) + ' for use with lea module')
        logging.debug(str(x_codon_dict))

        # Then convert to lea dict
        self.lea_codon_dict = {}
        for aa in self.codon_dict:
            self.lea_codon_dict[aa] = Lea.fromValFreqsDict(x_codon_dict[aa])

        # Logging
        logging.debug('')
        logging.debug('\nLea codon dict: ')
        logging.debug('Note: The Lea Module creates fractions which it will use to pick triplets, '
                      'where there is only 1 triplet (and where other triplets have been removed due to the cuttoff), the fraction will be 1')

        logging.debug(str(self.lea_codon_dict))
        logging.info('')

    def check_codon_table(self, test_codon_table=empty_codon_table):

        logging.info('')
        logging.info('---Check codon table uses correct triplet codes---')

        test_pass = True
        logging.info('Checking...')
        for aa in self.codon_dict:
            triplets_in_codon_table = []
            triplets_in_test_table = []

            for triplet in self.codon_dict[aa]:
                triplets_in_codon_table.append(triplet.upper())
            for triplet in test_codon_table[aa]:
                triplets_in_test_table.append(triplet.upper())

            if sorted(triplets_in_codon_table) != sorted(triplets_in_test_table):
                logging.warning('Codon Table uses incorrect triplets: ')
                logging.warning('Should be: ' + str(sorted(triplets_in_test_table)))
                logging.warning('But is: ' + str(sorted(triplets_in_codon_table)))
                test_pass = False
                raise NameError('Codon Table uses incorrect triplets', str(sorted(triplets_in_codon_table)))

        logging.info('Passed')
        logging.info('')

        return test_pass

class Sequence_Generator:

    def __init__(self, codon_table, check_dna_to_protein=True):
        self.codon_table = codon_table
        self.motifs_to_avoid = []
        self.check_dna_to_protein = check_dna_to_protein

    def optimise(self, protein_seq, log=False):
        # Make protein_seq all uppercase
        protein_seq = protein_seq.upper()

        # Set variable to hold new dna sequence
        dna_seq = ""

        # Check protein seq is valid
        check_protein_seq(protein_seq)

        # Iterate over the protein sequence, for each amino acid randomly select a triplet from the lea codon table,
        for i in range(len(protein_seq)):
            aa = protein_seq[i]
            dna_seq += self.codon_table.lea_codon_dict[aa].random()

        # Log new DNA sequence
        if log == False:
            logging.debug('DNA Seq Generated: ' + str(dna_seq))
        else:
            logging.info('')
            logging.info('---Optimise DNA sequence using codon table---')
            logging.info('DNA Seq Generated: ' + str(dna_seq))
            logging.info('')

        # Check that the new DNA sequence can be translated back to the same protein sequence
        if self.check_dna_to_protein == True:
            check_dna_back_translation(dna_seq, protein_seq)

        return dna_seq

    def optimise_from_dna(self, dna_seq):
        protein_seq = translate(dna_seq)
        return self.optimise(protein_seq)

    def redesign_section(self, dna_seq, start, end):

        window = dna_seq[start:end]
        new_window = self.optimise_from_dna(window)
        dna_seq = dna_seq[:start] + new_window + dna_seq[end:]

        return dna_seq

    def motif_removal(self, dna_seq, motifs):

        for motif_seq in motifs:
            list_of_motif_starts = MotifFinder.find_motif_starts(dna_seq, motif_seq)
            dna_seq = MotifFinder.remove_motif(dna_seq, list_of_motif_starts, motif_seq, self)

        return dna_seq

    def remove_high_GC_windows(self, dna_seq, window_size, GC_cuttoff):
        # For IDT windows of 100bp must be less than 76% GC
        # For IDT windows of 600bp must be less than 68% GC

        logging.info('')
        logging.info('---Checking GC Content for ' + str(window_size) + 'bp windows is less than ' + str(GC_cuttoff) + '% ---')

        for i in range(len(dna_seq) - window_size):
            start = i
            end = i + window_size
            gc_window = GCTools.calc_GC_content(dna_seq[start:end])

            if gc_window > GC_cuttoff:
                logging.info(str(window_size) + " bp window starting " + str(start) + " has a GC content of " + str(gc_window) + '.. redesigning')

                start, end = return_window_in_frame(dna_seq, i, i + window_size)
                window = GCTools.reduce_GC(dna_seq[start:end], GC_cuttoff, self)
                dna_seq = dna_seq[0:start] + window + dna_seq[end:len(dna_seq)]

        return dna_seq

    def minimise_mfe_windows(self, dna_seq, energy_limit, window_size, window_move, max_loops=10, max_iterations=100):
        count = 0
        changed = True

        while (changed == True) and (count < max_loops):
            dna_seq_old = dna_seq
            dna_seq = MFE.minimise_mfe_in_windows(dna_seq, energy_limit, window_size,
                                                      window_move, self, iterations=max_iterations)
            count +=1

            if dna_seq == dna_seq_old:
                changed=False

        if count >= max_loops:
            logging.warning("Max Loops reaching when minimising MFE in windows")

        return dna_seq

    def minimise_mfe_five_prime(self, dna_seq, energy_limit, five_prime_length=20, max_iterations=100):

        dna_seq = MFE.minimise_mfe(dna_seq, 0, five_prime_length, self, energy_limit, iterations=max_iterations)
        return dna_seq



if __name__ == "__main__":
    """ For testing this module """

    logging.basicConfig(level=logging.INFO)

    # Test directories
    working_dir = ROOT_DIR + '/Data/'

    # Single sequence test
    codon_table = CodonTable(codon_tables_dir=working_dir, json_file='Escherichia_coli_K12.json', low_cuttoff=0.1)
    seq_gen = Sequence_Generator(codon_table, check_dna_to_protein=False)

    protein_seq = "MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**"

    dna_seq = seq_gen.optimise(protein_seq, log=True)

    motifs_to_remove = ["ATTtt", 'GGaAt', "TaAAt"]
    dna_seq = seq_gen.motif_removal(dna_seq, motifs_to_remove)
    # dna_seq = seq_gen.remove_high_GC_windows(dna_seq, 100, 67)
    # dna_seq = seq_gen.minimise_mfe_windows(dna_seq, -10, 30, 5)
    # dna_seq = seq_gen.minimise_mfe_five_prime(dna_seq, -5)

    print('Finished')


    # Codon context
    # Sequential codons
    # NGG codons
    # SD site






    




