import logging
from Bio import SeqIO

# ------------------------------------
# Functions to import/export sequences
# ------------------------------------
def import_protein_seqs(location, file, file_type):

    logging.info('')
    logging.info('---Import protein sequences from ' + str(file_type) + '---')

    dict_of_sequences = {}

    for record in SeqIO.parse(location + file, file_type):
        name = record.id
        if name in dict_of_sequences:
            name += '(1)'
        logging.info('Imported ' + str(name))

        seq = record.seq

        dict_of_sequences[name] = seq

    logging.info('')

    return dict_of_sequences

def output_dna_sequences(location_output, file_name_output, dict_of_sequences):
        file = open(location_output + file_name_output, 'w')

        for name in dict_of_sequences:
            file.write('\n' + '>' + name)
            file.write('\n' + dict_of_sequences[name] + '\n')