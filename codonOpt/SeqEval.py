from codonOpt.ViennaRNA_wrapper import RNA_fold_output


""" Motifs"""
def find_motifs(dna_string, motif_string):

    motifs_binary = '0' * len(dna_string)
    count = 0

    for i in range(len(dna_string)):
        window = dna_string[i:i+len(motif_string)]

        if window == motif_string:
            motifs_binary = motifs_binary[:i] + ('1'*len(motif_string)) + motifs_binary[(i+len(motif_string)):]
            count += 1


    return_dict = {'motifs_binary' : motifs_binary,
                   'count' : count}

    return return_dict

def evaluate_motifs(dna_string, list_motif_strings):

    motif_evaluation_dict = {}
    count = 0

    for motif in list_motif_strings:
        motif_evaluation_dict[motif] = find_motifs(dna_string, motif)
        count += 1

    motif_evaluation_dict['total_count'] = count

    return motif_evaluation_dict

""" GC content """
def calc_GC_content(sequence):
    # This function will calculate and return the GC content of a sequence

    seq_length = len(sequence)
    count_of_GC = 0

    for i in range(seq_length):
        nucleotide = sequence[i]
        if nucleotide == "G" or nucleotide == "C":
            count_of_GC += 1

    GC_content = (count_of_GC / seq_length)*100
    GC_content = round(GC_content, 2)

    return GC_content

def map_GC(dna_string, window_size, window_move=1):
    # This function calculates GC content in windows, moving along the sequence 1bp at a time.
    bp_list = []
    gc_content_list = []

    for i in range(0, len(dna_string)-window_size, window_move):
        window = dna_string[i:i+window_size]
        gc_window = calc_GC_content(window)
        gc_content_list.append(gc_window)
        bp_list.append(i)

    return (bp_list, gc_content_list)

def gc_windows_above_threshold(list_gc_windows, threshold):

    return sum(gc < threshold for gc in list_gc_windows)

def evaluate_gc_content(dna_string, windows=((10,1,70), (20,1,70), (50,1,70), (100,1,70))):

    return_dict = {}
    return_dict['total_gc'] = calc_GC_content(dna_string)

    count = 0
    for window in windows:
        gc_windows =  map_GC(dna_string, window[0], window_move=window[1])[1]
        num_above_threshold = gc_windows_above_threshold(gc_windows, window[2])
        count += num_above_threshold

        return_dict['gc_' + str(window[0]) + '_windows'] = gc_windows
        return_dict[str(window[0]) + 'bp_windows_above_' + str(window[2])] = num_above_threshold

    return_dict['total_gc_windows_above_threshold'] = count

    return return_dict

""" Secondary Structure """
def calc_mfe(dna_string):
    # This function will calculate the minimum folding energy (MFE) of a given sequence - uses ViennaRNA to do this
    energy = RNA_fold_output(dna_string).energy
    return energy

def calc_mfe_windowed(dna_string, window_size=30, window_move=5):
    # This function will generate a list of MFE's.
    # It will move along a given DNA sequence calculating MFE's a sequence then length of window_size, moving by window_move each time.

    bp_list = []
    mfe_list = []

    for i in range(0, len(dna_string)-window_size, window_move):
        window = dna_string[i:i+window_size]
        energy_window = calc_mfe(window)
        mfe_list.append(energy_window)
        bp_list.append(i)

    return (bp_list, mfe_list)

def mfe_windows_above_threshold(list_of_window_energies, energy_threshold):

    return sum(energy < energy_threshold for energy in list_of_window_energies)

def evaluate_mfe(dna_string, window_sizes=((10,5,-10), (20,10,-10), (100,50,-10))):

    return_dict={}
    return_dict['total_mfe'] = calc_mfe(dna_string)

    count = 0
    for sizes in window_sizes:

        mfe_windows = calc_mfe_windowed(dna_string, window_size=sizes[0], window_move=sizes[1])[1]
        num_above_threshold = mfe_windows_above_threshold(mfe_windows, sizes[2])
        count += num_above_threshold

        return_dict['mfe_' + str(sizes[0]) + '_windows'] = mfe_windows
        return_dict[str(sizes[0]) + 'bp_windows_below_' + str(sizes[2])] = num_above_threshold

    return_dict['total_mfe_windows_above_threshold'] = count

    return return_dict


""" Codon Usage """
# RSCU
# Least squares?
# Some sort of inverse sigmoid function?



""" How repetitive is the sequence """
# Repeats, complementary sequence repeats ect
# Look at IDT



""" Domain boundaries"""
# Need to use something like pfam to find domain boundaries
# Align this with codon usage, secondary structure, sd sites
# pc liklihood of being a pause site?





# TODO Add codon usage metrics
# TODO Add repeat sequence metrics
# TODO Neural network to say whether this gene comes from X organism or not!!!


""" Other things to look at 
Consensus splice sites
Cryptic splice sites
SD sequences
TATA boxes
Termination signals
Artificial recombination sites

RNA instability motifs
Ribosomal entry sites
Repetitive sequences

Premature poly(A) sites
"""