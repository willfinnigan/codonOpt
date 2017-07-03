from codonOpt.global_vars import DNA_to_aa, empty_codon_table, amino_acids_dict
from Bio import SeqIO, Entrez
import logging
from copy import deepcopy
import json
import os.path
import pickle
Entrez.email = 'wjafinnigan@gmail.com'


def count_codons(codon_dict, cds):
    """
    Take an codon_dict (like in the format of empty_codon_table), and a cds.
    Count the triplet use in the CDS and add it to the codon_dict
    Return the codon_dict
    """

    cds_string = str(cds.dna_seq).upper()
    logging.debug(str(cds.name))
    for i in range(0, len(cds_string), 3):
        codon = cds_string[i:i + 3]
        if len(codon) == 3:
            aa = DNA_to_aa[codon]
            codon_dict[aa][codon] += 1
            logging.debug('aa: ' + str(aa) + ', Codon: ' + str(codon) + ', Count: ' + str(codon_dict[aa][codon]))

    return codon_dict

def count_amino_acids(amino_acid_dict, cds):
    """
    Take an amino_acid_dict (like in global_vars), and a cds.
    Count the amino acid use in the CDS and add it to the amino_acid_dict
    Return the amino_acid_dict
        """

    cds_string = str(cds.dna_seq).upper()
    for i in range(0, len(cds_string), 3):
        codon = cds_string[i:i + 3]
        if len(codon) == 3:
            amino_acid = DNA_to_aa[codon]
            amino_acid_dict[amino_acid] += 1

    return amino_acid_dict

def check_aa_and_codon_dicts_match_counts(codon_dict, amino_acid_dict):

    def check_aa_equals_total_codons(triplet_dict, amino_acid_no):
        # total triplets, then check it is equal to aa no
        total = 0
        for triplet in triplet_dict:
            total += triplet_dict[triplet]

        if total == amino_acid_no:
            return True
        else:
            return False

    for aa in codon_dict:
        amino_acid_no = amino_acid_dict[aa]
        triplet_dict = codon_dict[aa]
        if check_aa_equals_total_codons(triplet_dict, amino_acid_no) == False:
            print("------  WARNING, AMINO ACID COUNT DOES NOT EQUAL TRIPLET COUNT -----")
            print('AA count = ', amino_acid_no)
            print('triplet_dict = ', triplet_dict)

def convert_to_relative_frequency(codon_dict):
    relative_dict = deepcopy(empty_codon_table)

    for aa in relative_dict:
        total = 0
        for codon in relative_dict[aa]:
            total += codon_dict[aa][codon]
        for codon in relative_dict[aa]:
            if total != 0:
                relative_dict[aa][codon] = round(codon_dict[aa][codon]/total, 3)
            else:
                relative_dict[aa][codon] = float('nan')

    return relative_dict

def save_pickle(file, name, dir):
    file_name = dir + name + '.p'
    pickle.dump(file, open(file_name, "wb"))
    logging.info('Saved ' + str(name) + ' as pickle at ' + str(file_name))

def load_pickle(self, filename, dir):
    file = dir + filename
    loaded_pickle = pickle.load(open(file, "rb"))
    logging.info('Loaded pickle from ' + str(file))
    return loaded_pickle


class CDS():

    def __init__(self, feature, genome, preseq_distance=40):

        # Sequence information
        self.feature = feature
        self.seq_info = feature.extract(genome)
        self.qualifiers = feature.qualifiers
        self.name = self.qualifiers.get("product")
        self.location = feature.location
        self.dna_seq = self.seq_info.seq
        self.protein_seq = self.qualifiers.get("translation")

        # Position
        self.gene_start = self.location.start.position
        self.gene_end = self.location.end.position
        self.gene_strand = feature.strand
        self.preseq_start = self.gene_start - preseq_distance
        self.preseq_end = self.gene_start
        if self.gene_strand == 1:
            self.preseq = genome.seq[self.preseq_start:self.preseq_end]
        elif self.gene_strand == -1:
            self.preseq = genome.seq[self.preseq_start:self.preseq_end].reverse_complement()

        # Codon usage
        self.codon_usage_total = deepcopy(empty_codon_table)
        self.codon_usage_relative = deepcopy(empty_codon_table)
        self.amino_acid_usage = deepcopy(amino_acids_dict)

    def check_seq_only_contains_dna(self):
        passed = True
        list_of_bases = ['A', 'C', 'G', 'T']

        for letter in self.dna_seq:
            if letter.upper() not in list_of_bases:
                passed = False
                logging.warning('CDS ' + str(self.name) + ' contains non DNA letter ' + str(letter) + ' - removed CDS from analysis')

        return passed

class Genome():

    def __init__(self):

        self.genbank = ''
        self.min_gene_size = 0
        self.cds_list = []

        self.total_codon_usage_dict = deepcopy(empty_codon_table)
        self.relative_codon_usage_dict = deepcopy(empty_codon_table)
        self.amino_acid_dict = deepcopy(amino_acids_dict)

    def index_cds(self):
        """ Index_cds goes through all features in genbank, if cds and bigger than min size adds to self.cds_list"""
        logging.info('')
        logging.info("Index CDS's, minimum gene size is " + str(self.min_gene_size))
        count = 0

        for i,feature in enumerate(self.genbank.features):
            if feature.type=='CDS':
                new_cds = CDS(feature, self.genbank)
                if len(str(new_cds.dna_seq)) >= self.min_gene_size:
                    if new_cds.check_seq_only_contains_dna() == True:
                        self.cds_list.append(new_cds)
                        count+=1
                logging.debug(str(new_cds.name))
        logging.info('Indexed ' + str(count) + " CDS's, total now " + str(len(self.cds_list)))

    def analyse_codon_usage(self):
        logging.info('')
        logging.info('Analyse Genome Wide Codon Usage')

        for cds in self.cds_list:
            self.total_codon_usage_dict = count_codons(self.total_codon_usage_dict, cds)
            self.amino_acid_dict = count_amino_acids(self.amino_acid_dict, cds)
            logging.debug('CDS name: ' + str(cds.name))
            logging.debug('Current codon table= ' + str(self.total_codon_usage_dict))
            logging.debug('Current amino acid table= ' + str(self.amino_acid_dict))

        check_aa_and_codon_dicts_match_counts(self.total_codon_usage_dict, self.amino_acid_dict)

        self.relative_codon_usage_dict = convert_to_relative_frequency(self.total_codon_usage_dict)

        logging.info('Total codon table= ' + str(self.total_codon_usage_dict))
        logging.info('Total amino acid table= ' + str(self.amino_acid_dict))
        logging.info('Relative codon usage table= ' + str(self.relative_codon_usage_dict))

    def load_all_genbanks_in_folder(self, folder):
        logging.info('Load all genbank files in dir ' + str(folder))
        logging.info('Files in dir: ' + str(os.listdir(folder)))
        for filename in os.listdir(folder):
            if filename.endswith(".gb"):
                self.load_genbank(folder, filename)

    def save_relative_codon_table_as_json(self, file_name, dir=''):
        file = file_name + '.json'

        with open(dir+file, 'w') as fp:
            json.dump(self.relative_codon_usage_dict, fp)

        logging.info('Saved relative codon usage dict as a json at ' + str(dir) + str(file))

    def load_genbank(self, folder, genbank):
        logging.info('Load gb: ' + str(genbank))
        self.genbank = SeqIO.read(open(folder+genbank,"r"), "genbank")
        self.index_cds()













if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    name = 'Escherichia_coli_K12'
    dir = '/Users/Will/Documents/codonOpt_data/GenomeTmp/Escherichia_coli_K12/Escherichia_coli_str_K12_substr_MG1655/'
    codon_table_dir = '/Users/Will/Documents/codonOpt_data/Codon_tables/'
    pickle_dir = "/Users/Will/Documents/codonOpt_data/PickledGenomes/"

    genome = Genome()
    genome.load_all_genbanks_in_folder(dir)
    genome.analyse_codon_usage()
    genome.save_relative_codon_table_as_json(name, dir=codon_table_dir)
    save_pickle(genome, name, pickle_dir)

