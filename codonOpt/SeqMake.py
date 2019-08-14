import json
import logging
import math
import lea
from codonOpt.global_vars import empty_codon_table, uniform_codon_table, amino_acids_list, DNA_to_aa


class CodonTable():

    def __init__(self, json_file='', low_cuttoff=0.0):

        self.codon_dict = {}
        self.codon_dict_to_make_lea_dict = {}
        self.lea_codon_dict = {}
        self.low_cuttoff = low_cuttoff
        self.rscu_dict = {}

        if json_file != '':
            self.codon_dict = self.load_codon_table_from_json(json_file)
            self.check_codon_table(self.codon_dict)
            self.remove_rare_codons_for_lea_dict()
            self.make_lea_dict()


    """ Main functions to make a default codon table using lea dict """
    def load_codon_table_from_json(self, json_file):

        logging.info('---Load Codon Table---')

        # Open to codon table json.
        with open(json_file, 'r') as fp:
            codon_dict = json.load(fp)

        # Make nested dicts from json into dicts
        for aa in codon_dict:
            codon_dict[aa] = dict(codon_dict[aa])

        # Logging
        logging.info('Loaded codon dict from json at: ' + str(json_file))
        logging.info(str(codon_dict))

        return codon_dict

    def remove_rare_codons_for_lea_dict(self):

        logging.info('---Make Lea Codon Table for optimisation---')
        logging.info('Low cuttoff has been set at ' + str(self.low_cuttoff) + '  (codons used less than ' + str(self.low_cuttoff*100) + '% of the time)' )

        self.codon_dict_to_make_lea_dict = {}
        count = 0
        for aa in self.codon_dict:
            self.codon_dict_to_make_lea_dict[aa] = {}
            for triplet in self.codon_dict[aa]:
                if self.codon_dict[aa][triplet] >= self.low_cuttoff:
                    self.codon_dict_to_make_lea_dict[aa][triplet] = int((self.codon_dict[aa][triplet] * 10000))
                else:
                    logging.info('Amino Acid: ' + str(aa) + ', triplet: ' + str(triplet) + ' has freq of ' + str(self.codon_dict[aa][triplet]) + ' (less than ' + str(self.low_cuttoff) + '), so not included' )
                    count += 1

            self.check_triplet_dict_contains_triplets(self.codon_dict_to_make_lea_dict[aa], aa)

        logging.info(str(count) + ' triplets lower than cuttoff ' + str(self.low_cuttoff) + ' not included')
        logging.debug('Codon dict multiplied by ' + str(10000) + ' for use with lea module')
        logging.debug(str(self.codon_dict_to_make_lea_dict))

    def make_lea_dict(self):
        # Then convert to lea dict
        self.lea_codon_dict = {}
        for aa in self.codon_dict:
            self.lea_codon_dict[aa] = lea.pmf(self.codon_dict_to_make_lea_dict[aa])

        # Logging
        logging.debug('\nLea codon dict: ')
        logging.debug('Note: The Lea Module creates fractions which it will use to pick triplets, '
                      'where there is only 1 triplet (and where other triplets have been removed due to the cuttoff), the fraction will be 1')

        logging.debug(str(self.lea_codon_dict))


    """ Functions to check codon tabel is correct """
    def check_triplet_dict_contains_triplets(self, triplet_dict, aa):
        if len(triplet_dict) == 0:
            logging.warning('Amino Acid ' + str(aa) + ' has 0 triplets above low cuttoff')
            return False

        return True

    def check_codon_table(self, codon_table, test_codon_table=empty_codon_table):

        logging.info('---Check codon table uses correct triplet codes---')

        test_pass = True
        logging.info('Checking...')
        for aa in codon_table:
            triplets_in_codon_table = []
            triplets_in_test_table = []

            for triplet in codon_table[aa]:
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

        return test_pass


    """ Other functions """
    def remove_rare_codons_for_another_organism(self, codon_tables_dir = '', json_file='', low_cuttoff=0.0):
        """
         This function will remove rare codons for a second organism in self.codon_dict_to_make_lea_dict
         The lea_dict will then be remade
        """

        logging.info('')
        logging.info('---Load a second codon table from another organism to remove its rare codons---')
        logging.info('')

        second_codon_table = self.load_codon_table_from_json(codon_tables_dir, json_file)

        self.check_codon_table(second_codon_table)

        logging.info('')
        logging.info('---Remove the rare codons from second codon table---')
        for aa in second_codon_table:
            for triplet in second_codon_table[aa]:

                if second_codon_table[aa][triplet] < low_cuttoff:
                    if triplet in self.codon_dict_to_make_lea_dict[aa]:
                        logging.info('Amino Acid: ' + str(aa) + ', triplet: ' + str(triplet) + ' has freq of ' +
                                     str(second_codon_table[aa][triplet]) + ' (less than ' + str(self.low_cuttoff) +
                                     '), so will be removed.  Frequency in first table was ' + str(self.codon_dict[aa][triplet]))


                        self.codon_dict_to_make_lea_dict[aa].pop(triplet)

            self.check_triplet_dict_contains_triplets(self.codon_dict_to_make_lea_dict[aa], aa)

        self.make_lea_dict()

    def give_all_codons_below_cuttoff_equal_weighting(self):

        self.codon_dict_to_make_lea_dict = {}
        count = 0
        for aa in self.codon_dict:
            self.codon_dict_to_make_lea_dict[aa] = {}
            for triplet in self.codon_dict[aa]:
                if self.codon_dict[aa][triplet] >= self.low_cuttoff:
                    self.codon_dict_to_make_lea_dict[aa][triplet] = 1000
                else:
                    logging.info('Amino Acid: ' + str(aa) + ', triplet: ' + str(triplet) + ' has freq of ' + str(self.codon_dict[aa][triplet]) + ' (less than ' + str(self.low_cuttoff) + '), so not included' )
                    count += 1

            self.check_triplet_dict_contains_triplets(self.codon_dict_to_make_lea_dict[aa], aa)

        logging.info(str(count) + ' triplets lower than cuttoff ' + str(self.low_cuttoff) + ' not included')

        logging.debug('Every codon above cuttoff given a score of 1000')
        logging.debug(str(self.codon_dict_to_make_lea_dict))

    def calculate_rscu(self):
        pass

    def make_uniform(self):
        self.codon_dict = uniform_codon_table
        self.remove_rare_codons_for_lea_dict()
        self.make_lea_dict()

    def return_freq(self, triplet):
        aa = DNA_to_aa.get(triplet)
        freq = self.codon_dict[aa][triplet]

        return freq

class Sequence_Generator():

    def __init__(self, codon_table, check_dna_to_protein=False, logging=False):
        self.codon_table = codon_table
        self.motifs_to_avoid = []
        self.check_dna_to_protein = check_dna_to_protein
        self.optimisation_mode = 'default'
        self.logging = logging

    """ These functions can be called to do codon optimisation"""
    def optimise(self, protein_seq):

        if type(protein_seq) == str:
            protein_seq = Protein(protein_seq)
        protein_seq.make_uppercase()

        # Set variable to hold new dna sequence
        dna_seq = DNA([])
        dna_seq.protein = protein_seq

        # Iterate over the protein sequence, for each amino acid randomly select a triplet from the lea codon table,
        for i in range(len(protein_seq)):
            aa = protein_seq[i]
            dna_seq.append(self.codon_table.lea_codon_dict[aa].random())

        # Log new DNA sequence
        if self.logging == True:
            logging.info('---Optimise DNA sequence using codon table---')
            logging.info('DNA Seq Generated: ' + str(dna_seq))

        # Check that the new DNA sequence can be translated back to the same protein sequence
        if self.check_dna_to_protein == True:
            dna_seq.check_dna_equals_protein()

        return dna_seq

    def only_remove_rare_codons(self, dna_seq):
        logging.info('---Optimise DNA sequence by only removing rare codons using codon table---')
        count = 0

        for i in range(0,len(dna_seq)):
            triplet = dna_seq[i]
            aa = dna_seq.protein[i]
            logging.debug('Triplet: ' + str(triplet) + ',  AA: ' + str(aa))

            if triplet not in self.codon_table.codon_dict_to_make_lea_dict[aa]:
                count += 1
                dna_seq = self.redesign_section(dna_seq,i, i+1)
                logging.debug('Triplet ' + str(triplet) + ' not in codon table, replacing..')

        logging.info('Changed ' + str(count) + ' rare codons')

        return dna_seq

    def match_codon_usage(self, dna_seq):
        new_seq = DNA([])

        for i in range(len(dna_seq)):
            freq = dna_seq.freq[i]
            aa = DNA_to_aa.get(dna_seq[i])

            triplet_dict = self.codon_table.codon_dict[aa]
            best_difference = 1
            best_triplet = ''

            for triplet in triplet_dict:
                dif = triplet_dict[triplet] - freq
                dif = math.sqrt(dif*dif)

                if dif < best_difference:
                    best_difference = dif
                    best_triplet = triplet

            new_seq.append(best_triplet)

        return new_seq

class Protein(list):

    def __init__(self, protein, check_protein_valid=False):
        super(Protein, self).__init__()

        if type(protein) == str:
            self.input_as_string(protein)
        elif type(protein) == list:
            self.input_as_list(protein)

        if check_protein_valid == True:
            self.check_protein_seq()

    def input_as_string(self, protein_string):
        for i in range(0,len(protein_string)):
            aa = protein_string[i]
            self.append(aa)

    def input_as_list(self, protein_list):
        for aa in protein_list:
            self.append(aa)

    def check_protein_seq(self):
        """Checks the protein sequence against a list of amino acid characters.

        Returns:
            bool: The return value. True for success,  will raise a NameError is failed.

        """
        logging.debug('Checking protein sequence is valid...')

        for i in range(len(self)):
            self[i] = self[i].upper()

        for aa in self:
            if aa not in amino_acids_list:
                logging.warning('Non amino acid character found in protein sequence: ' + aa)
                raise NameError('Non amino acid character found in protein sequence')

        logging.debug('Passed')

        return True

    def make_uppercase(self):
        # Make protein_seq all uppercase
        for i in range(len(self)):
            self[i] = str(self[i]).upper()

    def __str__(self):
        s = ""
        for i in range(len(self)):
            s += self[i]

        return s

class DNA(list):

    def __init__(self, dna_seq):
        super(DNA, self).__init__()

        if type(dna_seq) == str:
            self.input_as_string(dna_seq)
        elif type(dna_seq) == list:
            self.input_as_list(dna_seq)

        self.protein = Protein(self.back_translate())
        self.freq = []



    def back_translate(self):
        translated_protein = Protein([])
        for triplet in self:
            translated_protein.append(DNA_to_aa.get(triplet))
        return translated_protein

    def check_dna_equals_protein(self):
        if self.back_translate() == self.protein:
            return True
        else:
            return False

    def get_frequencies(self, codon_table):
        self.freq = []

        for triplet in self:
            self.freq.append(round(codon_table.return_freq(triplet),2 ))

        return self.freq

    def input_as_string(self, dna_string):
        for i in range(0,len(dna_string),3):
            triplet = dna_string[i:i+3]
            self.append(triplet)

    def input_as_list(self, dna_list):
        for triplet in dna_list:
            self.append(triplet)

    def __str__(self):
        s = ""
        for i in range(len(self)):
            s += self[i]

        return s

def make_dna(protein, seq_gen):

    dna_seq = seq_gen.optimise(protein)

    return dna_seq