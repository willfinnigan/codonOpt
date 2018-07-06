import codonOpt
import logging

logging.basicConfig(level=logging.INFO)

codon_table = codonOpt.CodonTable(codon_tables_dir='Data/', json_file='Escherichia_coli_K12.json', low_cuttoff=0.1)
seq_gen = codonOpt.Sequence_Generator(codon_table)

dna_seq = codonOpt.DNA('ATGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGCGTGCAGTGGTGTTTGAAGGAAACAAAGAGCGCGTTGCCGTCAAAGAAGTGAATGCTCCGCAGGGTCTTCAGCACCCGAGGTTGGATGCACTTGTGCGTGTACACTTAGCTGGGATTTGTGGCTCCGACCTCCATCTGTACCACGGTAAAATCCCTGTGCTCCCAGGCTCAGTCTTGGGCCATGAATTCGTTGGTCAGGTTGAAGCCGTCGGTGAAGGTATCCAGGACTTACAACCGGGTGATTGGGTGGTCGGGCCATTTCACATCGCTTGCGGCACTTGTCCATATTGTCGTCGTCATCAGTATAACTTGTGTGAACGTGGCGGTGTGTACGGTTATGGCCCGATGTTCGGTAATCTGCAGGGAGCCCAAGCAGAGATTTTACGCGTGCCGTTCAGCAACGTCAACCTGCGTAAACTGCCTCCGAATCTGTCCCCGGAACGCGCTATCTTTGCGGGTGATATCCTGAGCACCGCATACGGAGGACTCATCCAGGGTCAGCTTCGCCCGGGGGATTCCGTAGCAGTTATTGGAGCGGGCCCTGTAGGGTTAATGGCCATCGAAGTGGCACAGGTCCTTGGCGCATCTAAAATCCTGGCCATTGACCGTATCCCGGAGCGTCTGGAACGCGCGGCGTCCCTCGGCGCAATTCCAATTAATGCCGAACAGGAAAATCCGGTCCGCCGTGTGCGTTCGGAGACTAACGATGAGGGCCCGGATTTAGTTCTGGAAGCGGTAGGCGGTGCGGCAACCCTTTCATTAGCGCTGGAAATGGTTCGCCCTGGCGGTCGCGTTTCTGCGGTGGGGGTCGATAACGCGCCATCATTCCCGTTCCCCCTTGCGTCTGGCCTGGTAAAAGATTTAACGTTTCGTATCGGCCTCGCAAACGTGCATCTCTACATCGATGCCGTCTTAGCCCTGTTAGCCTCTGGCCGTCTCCAACCCGAGCGTATTGTCAGCCACTATTTACCTCTGGAAGAAGCGCCGCGTGGCTATGAATTATTCGACCGTAAAGAAGCTTTAAAGGTGCTGCTGGTGGTGCGCGGGGGTGGTAGTGGAGATTACAAAGACGATGATGATAAGTGAAGGAGGAGGAGGAGGAGGAGGAGGTAA')

# Only change rare codons in the sequence
dna_seq = seq_gen.only_remove_rare_codons(dna_seq)

print(dna_seq)