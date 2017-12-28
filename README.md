# codonOpt
Example:

1.  Set up a codon table, and a sequence generator to use that codon table

    (`)codon_table = CodonTable(codon_tables_dir=codon_tables_dir, json_file='Escherichia_coli_K12.json', low_cuttoff=0.1)(`)
    (`)seq_gen = Sequence_Generator(codon_table)(`)
    
2.  Optimise your protein sequence    
    
    (`)protein_seq = "MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSD*"(`)
    (`)dna_seq = seq_gen.optimise(protein_seq)(`)
    
3.  Remove any restriction sites or other sequence motifs you don't want

    (`)motifs_to_remove = ["ATAA", 'GGAGA', "TACAG"](`)
    (`)dna_seq = seq_gen.motif_removal(dna_seq, motifs_to_remove)(`)
    
4.  Run other types of redesign routines to minimise folding energy, lower gc content ect..

    (`)dna_seq = seq_gen.remove_high_GC_windows(dna_seq, 100, 67)(`)
    (`)dna_seq = seq_gen.minimise_mfe_windows(dna_seq, -10, 30)(`)
    (`)dna_seq = seq_gen.minimise_mfe_five_prime(dna_seq, -5)(`)
