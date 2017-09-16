import logging
from codonOpt.AnalysisRedesignTools.GeneralFunctions import return_window_in_frame



def remove(dna_seq, start, motif, sequence_generator, max_iteration=500):

    motif = motif.upper()

    end = start + len(motif)
    window = dna_seq[start:end]

    start_in_frame, end_in_frame = return_window_in_frame(dna_seq, start, end)

    count = 0
    while (window == motif) and (count < max_iteration):
        count += 1

        dna_seq = sequence_generator.redesign_section(dna_seq, start_in_frame, end_in_frame)
        window = dna_seq[start:end]

        if count % 50 == 0:
            logging.warning('Reached iteration ' + str(count) + " while removing motif " + str(motif))

    logging.info('Position ' + str(start) + ' removal iterations: ' + str(count))
    if count == max_iteration:
        logging.warning(
            "Could not remove motif " + str(motif) + ", reached" + str(max_iteration) + " iterations")

    return dna_seq

def remove_motif(dna_seq, list_of_motif_starts, motif_seq, sequence_generator):
    """ Takes a list of motif starts sites and removes them all"""

    motif_seq = motif_seq.upper()

    logging.info("Removing motif " + str(motif_seq) + " from " + str(len(list_of_motif_starts)) + " locations")

    for motif_start in list_of_motif_starts:
        dna_seq = remove(dna_seq, motif_start, motif_seq, sequence_generator)

    logging.info('')
    return dna_seq

def find_motif_starts(dna_seq, motif):

    motif = motif.upper()

    logging.info('')
    logging.info('---Find and remove motif ' + str(motif) + '---')
    logging.info("Checking for presence of motif... " + str(motif))

    list_of_motif_starts = []

    for i in range(len(dna_seq)):
        window = dna_seq[i:i + len(motif)]

        if window == motif:
            logging.debug("Motif: " + str(motif) + " found at position " + str(i))
            list_of_motif_starts.append(i)

    return list_of_motif_starts