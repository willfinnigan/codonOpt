from codonOpt.global_vars import DNA_to_aa, empty_codon_table, amino_acids_dict, triplets_only_table
from codonOpt.AnalysisRedesignTools import MFE, GCTools
from Bio import SeqIO, Entrez
import logging
from copy import deepcopy
import json
import os.path
import pickle
from scipy.stats import chisquare
Entrez.email = 'wjafinnigan@gmail.com'



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

    logging.debug('codon_dict: ' + str(codon_dict))

    for aa in relative_dict:
        total = 0
        for codon in relative_dict[aa]:
            logging.debug('AA: ' + str(aa))
            logging.debug('Codon: ' + str(codon))
            total += codon_dict[aa][codon]

        for codon in relative_dict[aa]:
            if total != 0:
                relative_dict[aa][codon] = round(codon_dict[aa][codon]/total, 3)
            else:
                relative_dict[aa][codon] = float('nan')

    return relative_dict

def save_pickle(file, name, dir):
    file_name = dir + name + '.pickle'
    file_handle = open(file_name, "wb")
    pickle.dump(file, file_handle)
    file_handle.close()
    logging.info('Saved ' + str(name) + ' as pickle at ' + str(file_name))

def load_pickle(filename, dir):
    file = dir + filename  + '.pickle'
    file_handle = open(file, "rb")
    loaded_pickle = pickle.load(file_handle)
    file_handle.close()
    logging.info('Loaded pickle from ' + str(file))
    return loaded_pickle

def count_codons(codon_dict, seq):
    """
        Take an codon_dict (like in the format of empty_codon_table), and a cds.
        Count the triplet use in the seq and add it to the codon_dict
        Return the codon_dict
        """

    seq = str(seq.upper())

    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if len(codon) == 3:
            aa = DNA_to_aa[codon]
            codon_dict[aa][codon] += 1
            logging.debug('aa: ' + str(aa) + ', Codon: ' + str(codon) + ', Count: ' + str(codon_dict[aa][codon]))

    return codon_dict

def count_amino_acids(amino_acid_dict, seq):
    """
    Take an amino_acid_dict (like in global_vars), and a seq.
    Count the amino acid use in the CDS and add it to the amino_acid_dict
    Return the amino_acid_dict
        """

    seq = str(seq.upper())
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if len(codon) == 3:
            amino_acid = DNA_to_aa[codon]
            amino_acid_dict[amino_acid] += 1

    return amino_acid_dict

def count_codon_context(context_dict, seq):
    seq = str(seq.upper())
    context_codon = seq[0:3]

    for i in range(3, len(seq), 3):
        codon = seq[i:i + 3]

        if len(codon) == 3:
            aa = DNA_to_aa[codon]
            context_dict[context_codon][aa][codon] += 1
            context_codon = seq[i:i + 3]

    return context_dict



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
        self.codon_usage_relative = ''
        self.amino_acid_usage = deepcopy(amino_acids_dict)

        # Initiation Region Secondary structure
        self.ss_regions = []
        self.ss_mfes = []

        # Secondary structure windows
        self.ss_window_size_and_move = []
        self.ss_windows = []

        # GC windows
        self.gc_window_sizes= []
        self.gc_windows = []

        # Codon usage at start of the gene
        self.codon_usage_start_gene_sizes = []
        self.codon_usage_start_gene_tables = []
        self.amino_acid_usage_start_gene_tables = []

        # Min max
        self.min_max_start_sizes = []
        self.min_max_start_calcs = []

    def check_seq_only_contains_dna(self):
        passed = True
        list_of_bases = ['A', 'C', 'G', 'T']

        for letter in self.dna_seq:
            if letter.upper() not in list_of_bases:
                passed = False
                logging.warning('CDS ' + str(self.name) + ' contains non DNA letter ' + str(letter) + ' - removed CDS from analysis')

        return passed

    def calc_initiation_secondary_structures(self, ss_initiation_regions):
        logging.debug('----- Calculate Initiation SS for CDS----- ')

        self.ss_regions = ss_initiation_regions
        mfe=0
        for start,end in self.ss_regions:

            logging.debug('Start: ' + str(start) + " End: " + str(end))

            if start < 0:

                length_preseq = len(self.preseq)
                logging.debug("length preseq: " + str(length_preseq))

                if 0-start > length_preseq:
                    logging.warning('Initiation region start is bigger than then saved presequence length')
                    seq = ''
                    mfe=1

                else:
                    start_preseq = length_preseq+start
                    logging.debug("start preseq for mfe: " + str(start_preseq))

                    seq = self.preseq[start_preseq:]
                    logging.debug("preseq for mfe: " + str(seq))

            else:
                seq=''

            seq += self.dna_seq[0:end]
            logging.debug("all seq: " + str(seq))

            if mfe == 0:
                mfe = MFE.calc_mfe(str(seq))
                self.ss_mfes.append(mfe)

                logging.debug("MFE: " + str(mfe))
                logging.debug("")

    def calc_secondary_structure_windows(self, ss_window_and_move):
        logging.debug('----- Calculate SS Windows for CDS----- ')

        self.ss_window_size_and_move = ss_window_and_move
        self.ss_windows = []

        for window_size, window_move in self.ss_window_size_and_move:

            logging.debug(self.dna_seq)
            logging.debug('Window_size: ' + str(window_size) + " Window_move: " + str(window_move))
            mfe_windows = MFE.calc_mfe_windowed(self.dna_seq, window_size=window_size, window_move=window_move)
            self.ss_windows.append(mfe_windows)

            logging.debug("Number of windows measured: " + str(len(mfe_windows[0])))
            logging.debug("MFEs: " + str(mfe_windows[1]))
            logging.debug("")

    def calc_gc_windows(self, gc_window_sizes):
        logging.debug('----- Calculate GC windows for CDS----- ')

        self.gc_window_sizes = gc_window_sizes
        self.gc_windows = []

        for window_size in self.gc_window_sizes:
            logging.debug(self.dna_seq)
            logging.debug('Window_size: ' + str(window_size))

            gc_windows = GCTools.map_GC(self.dna_seq, window_size)
            self.gc_windows.append(gc_windows)

            logging.debug("Number of windows measured: " + str(len(gc_windows[0])))
            logging.debug("MFEs: " + str(gc_windows[1]))
            logging.debug("")

    def calc_codon_usage_at_start(self, start_sizes):
        logging.debug('----- Count codons and aa at start of gene----- ')

        self.codon_usage_start_gene_sizes = start_sizes
        self.codon_usage_start_gene_tables = []
        self.amino_acid_usage_start_gene_tables = []

        for size in self.codon_usage_start_gene_sizes:
            logging.debug('Size of start area: ' + str(size))
            codon_table = deepcopy(empty_codon_table)
            amino_acid_table = deepcopy(amino_acids_dict)

            seq = self.dna_seq[0:size]

            codon_table = count_codons(codon_table, seq)
            amino_acid_table = count_amino_acids(amino_acid_table, seq)
            logging.debug(str(codon_table))
            logging.debug(str(amino_acid_table))

            self.codon_usage_start_gene_tables.append(codon_table)
            self.amino_acid_usage_start_gene_tables.append(amino_acid_table)

class Genome():

    def __init__(self, name):

        self.name = name
        self.genbank = ''
        self.min_gene_size = 0
        self.cds_list = []

        self.total_codon_usage_dict = deepcopy(empty_codon_table)
        self.relative_codon_usage_dict = ''
        self.amino_acid_dict = deepcopy(amino_acids_dict)

        self.initiation_ss_regions = [(-10,10), (-20,20), (-30,30), (0,10), (0,20), (0,30), (0,40)]
        self.ss_window_size_and_move = [(10,5), (20, 5), (30, 5), (40, 5), (50, 5), (60, 5)]

        self.gc_window_sizes = [10,20,50,100]

        self.codons_at_start_of_genes_sizes = [15,30,60]
        self.codon_tables_at_start = []
        self.relative_codon_tables_at_start = []
        self.amino_acid_tabels_at_start = []

        self.codon_context_tables = deepcopy(triplets_only_table)
        for triplet in self.codon_context_tables:
            self.codon_context_tables[triplet] = deepcopy(empty_codon_table)


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

            # Count codons in cds
            cds.codon_usage_total = count_codons(cds.codon_usage_total, cds.dna_seq)
            cds.amino_acid_usage = count_amino_acids(cds.amino_acid_usage, cds.dna_seq)

            # Calc relative codon usage in cds
            check_aa_and_codon_dicts_match_counts(cds.codon_usage_total, cds.amino_acid_usage)
            cds.codon_usage_relative = convert_to_relative_frequency(cds.codon_usage_total)

            # Add codon and aa usage to global tables
            self.total_codon_usage_dict_in_seq = count_codons(self.total_codon_usage_dict, cds.dna_seq)
            self.amino_acid_dict = count_amino_acids(self.amino_acid_dict, cds.dna_seq)

            logging.debug('CDS name: ' + str(cds.name))
            logging.debug('Current codon table= ' + str(self.total_codon_usage_dict))
            logging.debug('Current amino acid table= ' + str(self.amino_acid_dict))

        # Calc global relative codon usage
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

    def analyse_initiation_secondary_structure(self):
        logging.info('')
        logging.info('Analyse Initiation Secondary Structure')

        for cds in self.cds_list:
            cds.calc_initiation_secondary_structures(self.initiation_ss_regions)

    def analyse_secondary_structures_in_windows(self):
        logging.info('')
        logging.info('Analyse Secondary Structure Windows')

        for cds in self.cds_list:
            cds.calc_secondary_structure_windows(self.ss_window_size_and_move)

    def calc_gc_windows(self):
        logging.info('')
        logging.info('Analyse GC Windows')

        for cds in self.cds_list:
            cds.calc_gc_windows(self.gc_window_sizes)

    def calc_codon_usage_at_start_of_genes(self):
        # Not including the start codon

        logging.info('')
        logging.info('Codon usage and AA usage at start of genes')

        self.codon_tables_at_start = []
        self.relative_codon_tables_at_start = []
        self.amino_acid_tabels_at_start = []

        for cds in self.cds_list:
            cds.calc_codon_usage_at_start(self.codons_at_start_of_genes_sizes)

        for size in self.codons_at_start_of_genes_sizes:
            new_codon_table = deepcopy(empty_codon_table)
            new_amino_acid_table = deepcopy(amino_acids_dict)

            for cds in self.cds_list:
                seq = cds.dna_seq[3:size]
                new_codon_table = count_codons(new_codon_table, seq)
                new_amino_acid_table = count_amino_acids(new_amino_acid_table, seq)

            check_aa_and_codon_dicts_match_counts(new_codon_table, new_amino_acid_table)

            new_relative_codon_table = convert_to_relative_frequency(new_codon_table)

            self.codon_tables_at_start.append(new_codon_table)
            self.relative_codon_tables_at_start.append(new_amino_acid_table)
            self.amino_acid_tabels_at_start.append(new_relative_codon_table)

            logging.debug(str(new_codon_table))
            logging.debug(str(new_amino_acid_table))
            logging.debug(str(new_relative_codon_table))

    def analyse_codon_context(self):

        for cds in self.cds_list:
            self.codon_context_tables = count_codon_context(self.codon_context_tables, cds.dna_seq)

        for triplet in self.codon_context_tables:
            logging.debug("Context Tripler: " + str(triplet))
            logging.debug(str(self.codon_context_tables[triplet]))



if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    name = 'Thermus_thermophilus_pTT27'
    dir = '/Users/Will/Dropbox/Python_Projects/codonOpt/codonOpt/test/example_data/GenomeTmp/Thermus thermophilus/Thermus_thermophilus_HB8/'
    codon_table_dir = '/Users/Will/Dropbox/Python_Projects/codonOpt/codonOpt/test/example_data/Codon_tables/'
    pickle_dir = "/Users/Will/Dropbox/Python_Projects/codonOpt/codonOpt/test/example_data/PickledGenomes/"

    genome = Genome(name)
    genome.load_all_genbanks_in_folder(dir)
    genome.analyse_codon_usage()
    genome.analyse_initiation_secondary_structure()
    genome.save_relative_codon_table_as_json(name, dir=codon_table_dir)
    save_pickle(genome, name, pickle_dir)
    genome = load_pickle(name, pickle_dir)
    print(genome.name)



