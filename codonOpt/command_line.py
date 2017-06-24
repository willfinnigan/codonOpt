from codonOpt.SeqMake import *

working_dir = input('What is the working dir? (Leave blank for current dir)')
codon_table_json = input('Enter codon table json file')
low_cuttoff = input('Enter low cuttoff for codon usage (0-1, defult=0.1)')
protein_seqs = input('Enter protein sequence fasta')
file_format = 'fasta'


input_dir = working_dir
output_dir = working_dir
codon_tables_dir = working_dir

dna_output = "optimised_dna_" + protein_seqs

codon_table = CodonTable(codon_tables_dir=codon_tables_dir, json_file=codon_table_json, low_cuttoff=low_cuttoff)
sequences = import_protein_seqs(input_dir, protein_seqs, file_format)
codon_optimise_seq_list(sequences, codon_table.lea_codon_dict)

for sequence in sequences:
        print()
        print(sequence.seq_name)
        print(sequence.dna_seq)

output_dna_sequence_list_as_fasta(output_dir, dna_output, sequences)


