import codonOpt
import logging
import time

logging.basicConfig(level=logging.DEBUG)

codon_table = codonOpt.CodonTable(json_file='Data/Escherichia_coli_K12.json', low_cuttoff=0.1)
seq_gen = codonOpt.Sequence_Generator(codon_table, check_dna_to_protein=True, logging=True)

protein_seq = codonOpt.Protein("MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**")
dna_seq = seq_gen.optimise(protein_seq)

dna_seq = seq_gen.motif_removal(dna_seq, ['ATGG'])

print(dna_seq)
print(dna_seq)