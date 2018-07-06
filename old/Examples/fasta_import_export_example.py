import codonOpt
import logging

# logging.basicConfig(level=logging.INFO)

# set up objects for codon optimisation
codon_table = codonOpt.CodonTable(codon_tables_dir='Data/', json_file='Escherichia_coli_K12.json', low_cuttoff=0.1)
seq_gen = codonOpt.Sequence_Generator(codon_table)

# load a fasta of protein sequences as a dictionary - dict_protein_seqs
dict_protein_seqs = codonOpt.import_protein_fasta('example_proteins.fasta', input_dir='Data/')

# iterate through dictionary running seq_gen optimisation functions
# save resulting dna in a dictionary - dict_dna_seqs
dict_dna_seqs = {}
for id in dict_protein_seqs:
    protein_seq = dict_protein_seqs[id]
    dna_seq = seq_gen.optimise(protein_seq)

    dict_dna_seqs[id] = dna_seq

# export dna sequences saved in dictionary, into a new fasta file (filename is timestamped)
codonOpt.export_dna_dict_fasta(dict_dna_seqs, output_dir='Data/')


