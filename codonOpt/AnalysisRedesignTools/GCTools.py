import random
from codonOpt.AnalysisRedesignTools.GeneralFunctions import translate




def calc_GC_content(sequence):
    """ This function will calculate and return the GC content of a sequence

    """

    seq_length = len(sequence)
    count_of_GC = 0

    for i in range(seq_length):
        nucleotide = sequence[i]
        if nucleotide == "G" or nucleotide == "C":
            count_of_GC += 1

    GC_content = (count_of_GC / seq_length)*100
    GC_content = round(GC_content, 2)

    return GC_content

def map_GC(dna_seq, window_size):
    """  This function calculates GC content in windows, moving along the sequence 1bp at a time.
    It returns a tuple of two lists.  1. bp's in middle of windows, 2.  window GC content
    """

    seq_length = len(dna_seq)
    if window_size % 2 != 0:
        window_size += 1

    half_window = int(window_size / 2)

    y_axis = []
    x_axis = []

    for i in range(half_window, seq_length-half_window):
        window = dna_seq[i-half_window:i+half_window]
        gc_window = calc_GC_content(window)
        y_axis.append(gc_window)
        x_axis.append(i)

    return (x_axis, y_axis)

def map_GC_windows_return_start_and_end(dna_seq, window_size=100):
    map_of_GC = map_GC(dna_seq, window_size)
    gc_high_region_list = []

    for i in range(len(map_of_GC)):
        half_window = window_size / 2
        middle = map_of_GC[i][0]
        start = middle - half_window
        end = middle + half_window

        gc_high_region = (start,end)
        print("GC high region = ", start, "to", end)
        gc_high_region_list.append(gc_high_region)

    return gc_high_region_list

def reduce_GC(dna_seq, target_GC, seq_gen):
    gc_content = calc_GC_content(dna_seq)
    translated = translate(dna_seq)

    while gc_content > target_GC:
        aa_no = random.randint(0, len(translated) - 1)
        aa = translated[aa_no]
        triplet_start = 3 * aa_no

        new_triplet = seq_gen.optimise(aa)
        new_dna_seq = dna_seq[0:triplet_start] + new_triplet + dna_seq[triplet_start + 3:len(dna_seq)]
        new_gc_content = calc_GC_content(new_dna_seq)

        if new_gc_content < gc_content:
            dna_seq = new_dna_seq
            gc_content = new_gc_content

    return dna_seq

