import codonOpt
import logging

logging.basicConfig(level=logging.DEBUG)

""" --- Make a sequence --- """
codon_table = codonOpt.CodonTable(json_file='Data/Escherichia_coli_K12.json', low_cuttoff=0.1)
seq_gen = codonOpt.Sequence_Generator(codon_table, check_dna_to_protein=True, logging=True)

protein_seq = codonOpt.Protein("MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**")
dna_seq = seq_gen.optimise(protein_seq)

print(dna_seq)



""" --- Evaluation ---"""
eval_motifs = codonOpt.evaluate_motifs(str(dna_seq), ['AAGG', 'TTCC'])
print(eval_motifs)

# windows are lists of (window_size, window_move, threshold)
eval_gc = codonOpt.evaluate_gc_content(str(dna_seq), windows=((10,1,70), (20,1,70), (50,1,70), (100,1,70)))
print(eval_gc)

# windows are lists of (window_size, window_move, threshold)
eval_mfe = codonOpt.evaluate_mfe(str(dna_seq), window_sizes=((20,5,-5), (50,10,-20), (100,50,-40)))
print(eval_mfe)