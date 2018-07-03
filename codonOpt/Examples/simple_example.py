import codonOpt
import logging
from codonOpt.AnalysisRedesignTools import GCTools

logging.basicConfig(level=logging.DEBUG)

codon_table = codonOpt.CodonTable(codon_tables_dir='Data/', json_file='Escherichia_coli_K12.json', low_cuttoff=0.1)
seq_gen = codonOpt.Sequence_Generator(codon_table)

protein_seq = "MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**"
dna_seq = seq_gen.optimise(protein_seq, log=True)

print(dna_seq)


gc = GCTools.show_bp_windows_over_threshold(dna_seq, 100, 64)
no_over = GCTools.percentage_bp_in_windows_over_gc_threshold(dna_seq, 100, 64)
print(no_over)

"""
---Optional---
Further functions allow the dna_seq to be fine tuned or further optimised
"""

"""
# Remove things like restriction sites or other motifs
motifs_to_remove = ["ATTtt", 'GGaAt', "TaAAt"]
dna_seq = seq_gen.motif_removal(dna_seq, motifs_to_remove)

# Lower GC in windows to below a certain level
dna_seq = seq_gen.remove_high_GC_windows(dna_seq, 100, 67)

# Minimise the RNA folding energy (MFE) using viennaRNA, in windows, to below a certain level
dna_seq = seq_gen.minimise_mfe_windows(dna_seq, -10, 30, 5)

# Minimise the RNA folding energy (MFE) using viennaRNA, of the five prime of the RNA.
dna_seq = seq_gen.minimise_mfe_five_prime(dna_seq, -5)

print(dna_seq)
"""