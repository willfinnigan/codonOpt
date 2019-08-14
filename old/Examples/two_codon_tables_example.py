import codonOpt
import logging

logging.basicConfig(level=logging.INFO)


# make e coli codon table
codon_table = codonOpt.CodonTable(codon_tables_dir='Data/', json_file='Escherichia_coli_K12.json', low_cuttoff=0.095)

# remove rare thermus codons
codon_table.remove_rare_codons_for_another_organism(codon_tables_dir='Data/', json_file='Tth codon table.json', low_cuttoff=0.095)

# set up sequence generator
seq_gen = codonOpt.Sequence_Generator(codon_table)

# codon optimise
protein_seq = "MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**"

dna_seq = seq_gen.optimise(protein_seq, log=True)

print(dna_seq)





""" 
--- Optional ---
This example results in a very high gc rich sequence being produced.  
Possibly the codon low cuttoffs need to be lowered?
One option might be to try and reduce the gc content using the 'remove_high_GC_windows' function.
Also runs of 4 G's or C's in a row could be avoided using the 'motif_removal' function
"""
dna_seq = seq_gen.remove_high_GC_windows(dna_seq, 100, 70)
dna_seq = seq_gen.remove_high_GC_windows(dna_seq, 100, 70)
dna_seq = seq_gen.motif_removal(dna_seq, ["GGGG", "CCCC"])
print(dna_seq)

