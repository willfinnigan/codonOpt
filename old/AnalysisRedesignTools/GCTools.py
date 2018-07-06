import logging
import random

from old.GeneralFunctions import translate


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

    logging.debug('bps: ' + str(map_of_GC[0]))
    logging.debug('gcs: ' + str(map_of_GC[1]))

    gc_high_region_list = []

    for i in range(len(map_of_GC[0])):
        half_window = window_size / 2
        middle = map_of_GC[0][i]
        start = middle - half_window
        end = middle + half_window

        gc_high_region_list.append((int(start),int(end)))


    return gc_high_region_list, map_of_GC[1]

def reduce_GC(dna_seq, target_GC, seq_gen):
    gc_content = calc_GC_content(dna_seq)
    translated = translate(dna_seq)

    count = 0
    initial_target = target_GC

    while gc_content > target_GC:

        count += 1
        aa_no = random.randint(0, len(translated) - 1)
        aa = translated[aa_no]
        triplet_start = 3 * aa_no

        new_triplet = seq_gen.optimise(aa)
        new_dna_seq = dna_seq[0:triplet_start] + new_triplet + dna_seq[triplet_start + 3:len(dna_seq)]
        new_gc_content = calc_GC_content(new_dna_seq)
        logging.debug('New GC content = ' + str(new_gc_content))

        if new_gc_content < gc_content:
            dna_seq = new_dna_seq
            gc_content = new_gc_content

        if count > 1000:
            target_GC += 1
            logging.info('After 1000 iterations, forced to raise target gc to ' + str(target_GC))
            count = 0



    return dna_seq

def number_gc_windows_higher_than_threshold(dna_seq, window_size, threshold):

    bp_list, gc_list = map_GC_windows_return_start_and_end(dna_seq, window_size)
    return sum(gc > threshold for gc in gc_list)

def show_bp_windows_over_threshold(dna_seq, window_size, threshold):

    bp_list, gc_list = map_GC_windows_return_start_and_end(dna_seq, window_size)

    logging.debug('bp_list = ' + str(bp_list))
    logging.debug('gc_list = ' + str(gc_list))

    high_bp = []
    high_gc = []

    for i in range(len(gc_list)):
        if gc_list[i] > threshold:
            logging.debug("Greater than threshold: gc of " + str(gc_list[i]) + ' at ' + str(bp_list[i]))
            high_bp.append(bp_list[i])
            high_gc.append(gc_list[i])


    gc_var = []
    for i in range(len(dna_seq)):
       gc_var.append(0)

    logging.debug('gc_var zeros = ' + str(gc_var))
    logging.debug('High_bp list = ' + str(high_bp))

    for start, end in high_bp:
        for i in range(start, end):
            gc_var[i] = 1

    logging.debug('gc_var = ' + str(gc_var))

    return gc_var

def number_bp_windows_over_threshold(dna_seq, window_size, threshold):
    gc_str = show_bp_windows_over_threshold(dna_seq, window_size, threshold)
    no_over = sum(gc == 1 for gc in gc_str)

    logging.debug('Bps in window over threshold = ' + str(no_over))

    return no_over

def percentage_bp_in_windows_over_gc_threshold(dna_seq, window_size, threshold):
    no_over = number_bp_windows_over_threshold(dna_seq, window_size, threshold)
    total_bp = len(dna_seq)

    pc_over = round(((no_over / total_bp) * 100), 2)
    logging.debug('Percentage bp over gc threshold in windows = ' + str(pc_over))

    return pc_over


