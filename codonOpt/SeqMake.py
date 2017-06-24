from lea import *
import logging
import json
from Bio import SeqIO
from codonOpt.global_vars import DNA_to_aa, amino_acids_list, empty_codon_table

# TODO comment everything, research docstrings
# TODO make it work on the command line
# TODO Make a new module for making json codon table from a genome.
# TODO Optimise from DNA

class CodonTable():
    """
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
        # Open to codon table json.
        with open(codon_tables_dir+json_file, 'r') as fp:
            self.codon_dict = json.load(fp)

        # Make nested dicts from json into dicts
        for aa in self.codon_dict:
            self.codon_dict[aa] = dict(self.codon_dict[aa])

        # Logging
        logging.info('Loaded codon dict from json at: ' + str(codon_tables_dir + json_file))
        logging.info(str(self.codon_dict))

    def load_lea_dict(self):
        x_codon_dict = {}
        count = 0
        for aa in self.codon_dict:
            x_codon_dict[aa] = {}
            for triplet in self.codon_dict[aa]:
                if self.codon_dict[aa][triplet] >= self.low_cuttoff:
                    x_codon_dict[aa][triplet] = int((self.codon_dict[aa][triplet] * 10000))
                else:
                    logging.debug('AA: ' + str(aa) + ', Triplet: ' + str(triplet) + ' has freq of ' + str(self.codon_dict[aa][triplet]) + ' so not included in lea dict, as this is less than the cuttoff of ' + str(self.low_cuttoff))
                    count += 1

            if len(x_codon_dict[aa]) == 0:
                logging.warning('AA ' + str(aa) + ' has 0 triplets above low cuttoff')

        # Logging

        logging.info(str(count) + ' triplets lower than cuttoff ' + str(self.low_cuttoff) + ' removed when generating the lea dict for choosing triplets')


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

    def check_codon_table(self, test_codon_table=empty_codon_table):
        test_pass = True
        logging.debug('Checking codon table uses correct triplets...')
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

        logging.debug('Passed')

        return test_pass

class Sequence():

    def __init__(self):
        self.original_dna_seq = ''
        self.dna_seq = ''
        self.protein_seq = ''
        self.seq_name = ''

# -- Functions for checking sequences are correct --
def check_protein_seq(protein_seq, amino_acids=amino_acids_list):
        logging.debug('Checking protein sequence is valid...')

        protein_seq = protein_seq.upper()
        for aa in protein_seq:
            if aa not in amino_acids:
                logging.warning('Non amino acid character found in protein sequence: ' + aa)
                raise NameError('Non amino acid character found in protein sequence')
                return False


        logging.debug('Passed')

        return True

def translate(dna_seq, DNA_to_aa = DNA_to_aa):
        """This function will translate a DNA sequence to a protein sequence using the dna_to_aa_table"""

        translated = ""
        dna_seq = dna_seq.upper()
        dna_length = len(dna_seq)

        if (dna_length%3) != 0:
            raise NameError('DNA Sequence to be translated is not divisable by 3')

        for i in range(0,dna_length,3):
            triplet = dna_seq[i:i+3]
            aa = DNA_to_aa.get(triplet)
            translated += aa
        return translated

def check_dna_back_translation(dna_seq, protein_seq):
        logging.debug('Checking translation back to protein sequence matches...')

        dna_seq = dna_seq.upper()
        protein_seq = protein_seq.upper()

        translated = translate(dna_seq)
        if translated != protein_seq:
            logging.warning('Translation of dna seq to protein seq does not give the original protein seq')
            logging.warning('Original   : ' + str(protein_seq))
            logging.warning('Translation: ' + str(translated))

            raise NameError('Translation of dna seq to protein seq does not give the original protein seq')
            return False

        logging.debug('Passed')
        return True

# -- Functions to generate new codon optimised DNA from a protein sequence --
def generate_seq(protein_seq, lea_dict):

    # Make protein_seq all uppercase
    protein_seq = protein_seq.upper()

    # Set variable to hold new dna sequence
    dna_seq = ""

    # Check protein seq is valid
    check_protein_seq(protein_seq)

    # Iterate over the protein sequence, for each amino acid randomly select a triplet from the lea codon table,
    for i in range(len(protein_seq)):
        aa = protein_seq[i]
        dna_seq += lea_dict[aa].random()

    # Log new DNA sequence
    logging.debug('DNA Seq Generated: ' + str(dna_seq))

    # Check that the new DNA sequence can be translated back to the same protein sequence
    check_dna_back_translation(dna_seq, protein_seq)

    return dna_seq

def codon_optimise_seq_list(sequence_obj_list, lea_dict):
    dna_seqs = []
    for seq in sequence_obj_list:
        seq.dna_seq = generate_seq(seq.protein_seq, lea_dict)
        dna_seqs.append(dna_seqs)

    return dna_seqs

# -- Functions to import/output sequences
def import_protein_seqs(location, file, file_type):
    list_of_seqs = []

    for record in SeqIO.parse(location+file, file_type):
        seq_obj = Sequence()
        seq_obj.seq_name = record.id
        seq_obj.protein_seq = record.seq
        list_of_seqs.append(seq_obj)

    return list_of_seqs

# This function needs a test
def output_dna_sequence_list_as_fasta(location_output, file_name_output, sequence_list):

    file = open(location_output + file_name_output, 'w')

    for seq in sequence_list:
        file.write('\n' + '>'+ seq.seq_name)
        file.write('\n' + seq.dna_seq + '\n')







if __name__ == "__main__":
    """ For testing this module """

    logging.basicConfig(level=logging.DEBUG)

    # Test individual sequence optimisation
    codon_tables_dir = 'example_data/'

    codon_table = CodonTable(codon_tables_dir=codon_tables_dir, json_file='Tth codon table.json', low_cuttoff=0.1)
    a_protein = "MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**"
    seq = generate_seq(a_protein, codon_table.lea_codon_dict)

    # Test importing set of fasta sequences and exporting optimised dna sequences
    input_dir = 'example_data/'
    output_dir = 'example_data/'
    codon_tables_dir = 'example_data/'

    codon_table = CodonTable(codon_tables_dir=codon_tables_dir, json_file='NC_000913.3.json', low_cuttoff=0.1)
    sequences = import_protein_seqs(input_dir, 'protein_seqs.fasta', 'fasta')
    codon_optimise_seq_list(sequences, codon_table.lea_codon_dict)
    for sequence in sequences:
        print()
        print(sequence.seq_name)
        print(sequence.dna_seq)
    output_dna_sequence_list_as_fasta(output_dir, 'test_dna.fasta', sequences)




