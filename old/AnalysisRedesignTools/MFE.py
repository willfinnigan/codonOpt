import logging
import random

import matplotlib.pyplot as plt
from codonOpt.AnalysisRedesignTools.GeneralFunctions import return_window_in_frame

from codonOpt.ViennaRNA_wrapper import RNA_fold_output


def random_length(a_list):
    length = list(range(len(a_list)))
    random.shuffle(length)
    return list(length)

def calc_mfe(dna_seq):
    # This function will calculate the minimum folding energy (MFE) of a given sequence - uses ViennaRNA to do this
    energy = RNA_fold_output(dna_seq).energy
    return energy

def calc_mfe_windowed(dna_seq, window_size=30, window_move=5):
    # This function will generate a list of MFE's.
    # It will move along a given DNA sequence calculating MFE's a sequence then length of window_size, moving by window_move each time.

    bp_list = []
    mfe_list = []

    for i in range(0, len(dna_seq)-window_size, window_move):
        window = dna_seq[i:i+window_size]
        energy_window = calc_mfe(window)
        mfe_list.append(energy_window)
        bp_list.append(i)

    return (bp_list, mfe_list)

def minimise_mfe(dna_seq, start, end, seq_gen, energy_limit, iterations=100):

    start_in_frame, end_in_frame = return_window_in_frame(dna_seq, start, end)
    energy_window = calc_mfe(dna_seq[start:end])

    count = 0
    while count <= iterations:
        if energy_window <= energy_limit:
            new_section = seq_gen.optimise_from_dna(dna_seq[start_in_frame:end_in_frame])
            dna_seq = dna_seq[0:start_in_frame] + new_section + dna_seq[end_in_frame:]
            window = dna_seq[start:end]
            energy_window = calc_mfe(window)

            count += 1

        elif energy_window > energy_limit:
            if count == 0:
                logging.info('MFE = ' + str(energy_window) + ', already lower than limit of ' + str(energy_limit))
            else:
                logging.info("In " + str(count) + " iterations, energy is now " + str(energy_window))
            count = iterations + 1

        if count == iterations:
            logging.warning("Warning, could not reduce energy below limit of " + str(energy_limit) +
                            " before max iterations " + str(iterations) + " reached")
            count += 1

    return dna_seq

def minimise_mfe_in_windows(dna_seq, energy_limit, window_size, window_move, seq_gen, iterations=100):
        logging.info('')
        logging.info("---  Scan for high secondary structure in windows of " + str(window_size) +
                     " bp, for energy more than " + str(energy_limit) + ", moving by " +
                     str(window_move) + " bp each time ---")

        list_starts, list_energies = calc_mfe_windowed(dna_seq, window_size=window_size, window_move=window_move)

        for i in random_length(list_energies):

            if list_energies[i] <= energy_limit:
                start = list_starts[i]
                end = list_starts[i] + window_size
                energy = list_energies[i]

                logging.info("Window at " + str(start) + " to " + str(end) + " energy is " +
                             str(energy) + "..reduce to below " + str(energy_limit))

                dna_seq = minimise_mfe(dna_seq, start, end, seq_gen, energy_limit, iterations=iterations)

        return dna_seq

def plot_mfe_windows(dna_seq, window_size=30, window_move=5):

    x_axis, y_axis = calc_mfe_windowed(dna_seq, window_size=window_size, window_move=window_move)
    fig = plt.figure()
    plt.title('Page One')
    plt.plot(x_axis, y_axis)
    plt.close()

    return fig

def count_windows_above_threshold(list_of_window_energies, energy_threshold):

    number_of_above_threshold_windows = sum(energy < energy_threshold for energy in list_of_window_energies)

    return number_of_above_threshold_windows













if __name__ == "__main__":
    seq = 'ATGAGAGGCGAGCGACGAGCGACGAGCGATCAGCTAGCTACTAC'
    mfe = calc_mfe(seq)
    print(mfe)

    fig1 = plot_mfe_windows(seq)
    fig2 = plot_mfe_windows(seq)
    save_pdf([fig1, fig2])
