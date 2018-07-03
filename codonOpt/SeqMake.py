import lea
import logging, json
import math

from codonOpt.global_vars import empty_codon_table, uniform_codon_table
from codonOpt.AnalysisRedesignTools import MotifFinder, GCTools, MFE
from codonOpt.AnalysisRedesignTools.GeneralFunctions import return_window_in_frame, translate, check_protein_seq, check_dna_back_translation


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
    """CodonTable class - uses codontable.lea_codon_dict for optimisation, loads from json file

    This is the CodonTable Class
    It has codon_dict, codon_dict_to_make_lea_dict, lea_codon_dict and low_cuttoff variables.

    When initialising it takes the following arguments:
        codon_tables_dir = "" - This is the directory to find the codon table json file.
        json_file = "" - This is the json file containing the codon usage information.
        low_cuttoff = 0 - This is the lowcuttoff which is used to remove codons with low frequency.

    If a json file is specified this will be loaded and automatically converted into a lea_codon_dict for optimisation.


    """

    def __init__(self, codon_tables_dir = '', json_file='', low_cuttoff=0.0):

        self.codon_dict = {}
        self.codon_dict_to_make_lea_dict = {}
        self.lea_codon_dict = {}
        self.low_cuttoff = low_cuttoff
        self.rscu_dict = {}

        if json_file != '':
            self.codon_dict = self.load_codon_table_from_json(codon_tables_dir, json_file)
            self.check_codon_table(self.codon_dict)
            self.remove_rare_codons_for_lea_dict()
            self.make_lea_dict()


    """ Main functions to make a default codon table using lea dict """
    def load_codon_table_from_json(self, codon_tables_dir, json_file):

        logging.info('')
        logging.info('---Load Codon Table---')
        json_path = codon_tables_dir+json_file

        # Open to codon table json.
        with open(json_path, 'r') as fp:
            codon_dict = json.load(fp)

        # Make nested dicts from json into dicts
        for aa in codon_dict:
            codon_dict[aa] = dict(codon_dict[aa])

        # Logging
        logging.info('Loaded codon dict from json at: ' + str(json_path))
        logging.info(str(codon_dict))
        logging.info('')

        return codon_dict

    def remove_rare_codons_for_lea_dict(self):

        logging.info('')
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
            self.lea_codon_dict[aa] = lea.Lea.fromValFreqsDict(self.codon_dict_to_make_lea_dict[aa])

        # Logging
        logging.debug('')
        logging.debug('\nLea codon dict: ')
        logging.debug('Note: The Lea Module creates fractions which it will use to pick triplets, '
                      'where there is only 1 triplet (and where other triplets have been removed due to the cuttoff), the fraction will be 1')

        logging.debug(str(self.lea_codon_dict))
        logging.info('')


    """ Functions to check codon tabel is correct """
    def check_triplet_dict_contains_triplets(self, triplet_dict, aa):
        if len(triplet_dict) == 0:
            logging.warning('Amino Acid ' + str(aa) + ' has 0 triplets above low cuttoff')
            return False

        return True

    def check_codon_table(self, codon_table, test_codon_table=empty_codon_table):

        logging.info('')
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
        logging.info('')

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
        for aa in self.codon_dict:
           triplet_dict = self.codon_dict[aa]

    def make_uniform(self):
        self.codon_dict = uniform_codon_table
        self.remove_rare_codons_for_lea_dict()
        self.make_lea_dict()

    def return_freq(self, triplet):
        aa = translate(triplet)
        freq = self.codon_dict[aa][triplet]

        return freq

#TODO convert seq_gen to use the DNA and Protein classes and get rid of the extra make_dna function
class Sequence_Generator():

    def __init__(self, codon_table, check_dna_to_protein=False):
        self.codon_table = codon_table
        self.motifs_to_avoid = []
        self.check_dna_to_protein = check_dna_to_protein
        self.optimisation_mode = 'default'

    def optimise_from_protein(self, protein_seq, log=False):

        if self.optimisation_mode == 'default':
            dna_seq = self.default_optimise(protein_seq, log=log)

        elif self.optimisation_mode == 'only_rare':
            print("Please run optimise_dna to only remove rare codons")

        elif self.optimisation_mode == 'match_codon_usage':
            print("Please run optimise_dna to match codon_usage")

        return dna_seq

    def optimise_from_dna(self, dna_seq, codon_table={}):

        if self.optimisation_mode == 'default':
            dna_seq = self.default_optimise_from_dna(dna_seq)

        elif self.optimisation_mode == 'only_rare':
            dna_seq = self.only_remove_rare_codons(dna_seq)

        elif self.optimisation_mode == 'match_codon_usage':
            print("Please run optimise_dna to match codon_usage")

    def default_optimise_from_protein(self, protein_seq, log=False):
        # Make protein_seq all uppercase
        for i in range(len(protein_seq)):
            protein_seq[i] = str(protein_seq[i]).upper()

        # Set variable to hold new dna sequence
        dna_seq = DNA([])

        # Check protein seq is valid
        if self.check_dna_to_protein == True:
            check_protein_seq(protein_seq)

        # Iterate over the protein sequence, for each amino acid randomly select a triplet from the lea codon table,
        for i in range(len(protein_seq)):
            aa = protein_seq[i]
            dna_seq.append(self.codon_table.lea_codon_dict[aa].random())

        # Log new DNA sequence
        if log == True:
            logging.info('')
            logging.info('---Optimise DNA sequence using codon table---')
            logging.info('DNA Seq Generated: ' + str(dna_seq))
            logging.info('')

        # Check that the new DNA sequence can be translated back to the same protein sequence
        if self.check_dna_to_protein == True:
            check_dna_back_translation(dna_seq, protein_seq)

        return dna_seq

    def default_optimise_from_dna(self, dna_seq):
        print(dna_seq)
        protein = translate(dna_seq)

        print(protein)
        # Create a codon optimised dna sequenced from the protein sequence.
        return self.default_optimise_from_protein(protein)

    def redesign_section(self, dna_seq, start_triplet, end_triplet):

        # Take a section of the dna_seq between start and end.
        window = DNA(dna_seq[start_triplet:end_triplet])
        print(window)

        # Codon optimise the window
        new_window = self.default_optimise_from_dna(window)

        # Replace the window section of dna_seq with the newly codon optimised new_window
        dna_seq = DNA(dna_seq[:start_triplet] + new_window + dna_seq[end_triplet:])

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
            aa = translate(dna_seq[i])

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

    """ ------- 
    These functions use the AnalysisRedesignTools
        -------"""
    def motif_removal(self, dna_seq, motifs):

        motifs_present = True

        while motifs_present == True:
            for motif_seq in motifs:
                list_of_motif_starts = MotifFinder.find_motif_starts(dna_seq, motif_seq)
                dna_seq = MotifFinder.remove_motif(dna_seq, list_of_motif_starts, motif_seq, self)

            # Check motifs = 0
            if MotifFinder.count_number_of_motifs(dna_seq, motifs) == 0:
                motifs_present = False


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
                logging.info('-- ' + str(window_size) + " bp window starting " + str(start) + " has a GC content of " + str(gc_window) + '.. redesigning to less than ' + str(GC_cuttoff))

                start, end = return_window_in_frame(dna_seq, i, i + window_size)
                window = GCTools.reduce_GC(dna_seq[start:end], GC_cuttoff, self)
                dna_seq = dna_seq[0:start] + window + dna_seq[end:len(dna_seq)]

        return dna_seq

    def minimise_mfe_windows(self, dna_seq, energy_limit, window_size, window_move, max_loops=10, max_iterations=100):
        count = 0
        changed = True

        # This will keep going until no more changes are made, or the max loops is reached
        while (changed == True) and (count < max_loops):
            dna_seq_old = dna_seq
            dna_seq = MFE.minimise_mfe_in_windows(dna_seq, energy_limit, window_size,
                                                      window_move, self, iterations=max_iterations)
            count +=1

            if dna_seq == dna_seq_old:
                changed=False


        if count >= max_loops:
            logging.warning("Max Loops reached when minimising MFE in windows")

        logging.info('')
        logging.info('--- MFE window minimisation complete, number of loops through dna sequence = ' + str(count) + ' ---')
        logging.info('')



        return dna_seq

    def minimise_mfe_five_prime(self, dna_seq, energy_limit, five_prime_length=20, max_iterations=100):
        logging.info('---Minimise five prime folding energy in first '
                     + str(five_prime_length) + ' to be less than '
                     + str(energy_limit) + ' ---')

        dna_seq = MFE.minimise_mfe(dna_seq, 0, five_prime_length, self, energy_limit, iterations=max_iterations)
        return dna_seq

class Protein(list):

    def __init__(self, protein):
        super(Protein, self).__init__()

        if type(protein) == str:
            self.input_as_string(protein)
        elif type(protein) == list:
            self.input_as_list(protein)

    def input_as_string(self, protein_string):
        for i in range(0,len(protein_string)):
            aa = protein_string[i]
            self.append(aa)

    def input_as_list(self, protein_list):
        for aa in protein_list:
            self.append(aa)

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
        protein = []
        for triplet in self:
            translated = translate(triplet)
            protein.append(translated)
        return protein

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
    dna_seq = []

    for aa in protein:
        dna_seq.append(seq_gen.optimise(aa))

    return dna_seq





