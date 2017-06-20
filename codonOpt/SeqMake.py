from lea import *
import logging
import json
from Bio import SeqIO


class CodonTable():

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

    def check_codon_table(self):
        test_codon_table = {"V": {"GTT": 0.031, "GTG": 0.6156, "GTC": 0.3363, "GTA": 0.017},
                              "G": {"GGC": 0.3996, "GGA": 0.0609, "GGG": 0.5145, "GGT": 0.025},
                              "F": {"TTT": 0.177, "TTC": 0.823},
                              "E": {"GAA": 0.125, "GAG": 0.875},
                              "N": {"AAT": 0.027, "AAC": 0.973},
                              "P": {"CCT": 0.068, "CCA": 0.023, "CCG": 0.176, "CCC": 0.733},
                              "Q": {"CAA": 0.13, "CAG": 0.87},
                              "M": {"ATG": 1.0},
                              "K": {"AAG": 0.915, "AAA": 0.085},
                              "T": {"ACT": 0.012, "ACC": 0.7267, "ACA": 0.011, "ACG": 0.2503},
                              "S": {"TCT": 0.024, "TCC": 0.4464, "AGT": 0.011, "TCA": 0.009, "TCG": 0.1101, "AGC": 0.3994},
                              "W": {"TGG": 1.0},
                              "A": {"GCG": 0.2268, "GCC": 0.7323, "GCA": 0.017, "GCT": 0.024},
                              "R": {"CGG": 0.4264, "CGC": 0.3403, "AGG": 0.1892, "CGA": 0.015, "AGA": 0.01, "CGT": 0.019},
                              "D": {"GAT": 0.048, "GAC": 0.952},
                              "L": {"TTA": 0.008, "CTA": 0.022, "CTT": 0.114, "CTC": 0.505, "CTG": 0.285, "TTG": 0.066},
                              "Y": {"TAT": 0.043, "TAC": 0.957},
                              "H": {"CAC": 0.95, "CAT": 0.05},
                              "I": {"ATC": 0.86, "ATA": 0.046, "ATT": 0.094},
                              "C": {"TGC": 0.95, "TGT": 0.05},
                              "*": {"TAG": 0.3676, "TAA": 0.1948, "TGA": 0.4376}}

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
        self.dna_seq = ''
        self.protein_seq = ''
        self.seq_name = ''


def check_protein_seq(protein_seq):
        logging.debug('Checking protein sequence is valid...')
        amino_acids = ['V', 'G', 'F', 'E', 'N', 'P', 'Q', 'M', 'K', 'T', 'S', 'W', 'A', 'R', 'D', 'L', 'Y', 'H', 'I', 'C', '*']

        protein_seq = protein_seq.upper()
        for aa in protein_seq:
            if aa not in amino_acids:
                logging.warning('Non amino acid character found in protein sequence: ' + aa)
                raise NameError('Non amino acid character found in protein sequence')
                return False


        logging.debug('Passed')

        return True

def translate(dna_seq):
        # This function will translate a DNA sequence to a protein sequence using the dna_to_aa_table
        DNA_to_aa = {
                        "GCT":"A", "GCG":"A", "GCC":"A", "GCA":"A",
                        "AGA":"R", "CGA":"R", "CGT":"R", "AGG":"R", "CGC":"R", "CGG":"R",
                        "AAT":"N", "AAC":"N",
                        "GAT":"D", "GAC":"D",
                        "TGT":"C", "TGC":"C",
                        "TAA":"*", "TAG":"*", "TGA":"*",
                        "CAA":"Q", "CAG":"Q",
                        "GAA":"E", "GAG":"E",
                        "GGT":"G", "GGA":"G", "GGC":"G", "GGG":"G",
                        "CAT":"H", "CAC":"H",
                        "ATA":"I", "ATT":"I", "ATC":"I",
                        "TTA":"L", "TTG":"L", "CTT":"L", "CTG":"L", "CTC":"L", "CTA":"L",
                        "AAA":"K", "AAG":"K",
                        "ATG":"M",
                        "TTT":"F", "TTC":"F",
                        "CCA":"P", "CCT":"P", "CCG":"P", "CCC":"P",
                        "TCA":"S", "AGT":"S", "TCT":"S", "TCG":"S", "AGC":"S", "TCC":"S",
                        "ACA":"T", "ACT":"T", "ACG":"T", "ACC":"T",
                        "TGG":"W",
                        "TAT":"Y", "TAC":"Y",
                        "GTA":"V", "GTT":"V", "GTC":"V", "GTG":"V"}

        translated = ""
        dna_seq = dna_seq.upper()
        dna_length = len(dna_seq)
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


def import_protein_seqs(location, file, file_type):
    list_of_seqs = []

    for record in SeqIO.parse(location+file, file_type):
        seq_obj = Sequence()
        seq_obj.seq_name = record.id
        seq_obj.protein_seq = record.seq
        list_of_seqs.append(seq_obj)

    return list_of_seqs

def output_dna_sequence_list_as_fasta(location_output, file_name_output, sequence_list):

    file = open(location_output + file_name_output, 'w')

    for seq in sequence_list:
        file.write('\n' + '>'+ seq.seq_name)
        file.write('\n' + seq.dna_seq + '\n')






if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Test individual sequence optimisation
    codon_tables_dir = 'example_data/'

    codon_table = CodonTable(codon_tables_dir=codon_tables_dir, json_file='Tth codon table.json', low_cuttoff=0.1)
    a_protein = "MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**"
    seq = generate_seq(a_protein, codon_table.lea_codon_dict)



    # Test importing set of fasta sequences and exporting optimised dna sequences
    input_dir = 'example_data/'
    output_dir = 'example_data/'
    codon_tables_dir = 'example_data/'

    codon_table = CodonTable(codon_tables_dir=codon_tables_dir, json_file='Tth codon table.json', low_cuttoff=0.1)
    sequences = import_protein_seqs(input_dir, 'protein_seqs.fasta', 'fasta')
    codon_optimise_seq_list(sequences, codon_table.lea_codon_dict)
    for sequence in sequences:
        print()
        print(sequence.seq_name)
        print(sequence.dna_seq)
    output_dna_sequence_list_as_fasta(output_dir, 'test_dna.fasta', sequences)


