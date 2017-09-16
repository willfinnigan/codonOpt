from codonOpt.global_vars import DNA_to_aa, empty_codon_table, amino_acids_dict
from Bio import SeqIO, Entrez
import logging
from copy import deepcopy
import json
import os.path
Entrez.email = 'wjafinnigan@gmail.com'

def get_genome(id, dir, name):
    name = remove_spaces_and_dots(name)
    file_name = dir + name + '.gb'

    handle = Entrez.efetch(db='nucleotide',id=id,rettype='gbwithparts' )
    print('File name = ', str(file_name))
    local_file = open(file_name, 'w')
    local_file.write(handle.read())
    handle.close()
    local_file.close()

def look_up_representative_assembly(name, only_search_representative_genomes):
    if only_search_representative_genomes == True:
        search_term = name + ' [ORGN] AND representative [PROP]'
    else:
        search_term = name + ' [ORGN]'

    logging.info('Search term = ' + str(search_term))

    handle = Entrez.esearch(db='assembly', term=search_term)
    result = Entrez.read(handle)
    handle.close()

    logging.info(str(result))
    logging.info(str(result['Count']))
    logging.info(str(result["IdList"]))

    return result

def get_ids_from_assembly_lookup_result(result):

    id_dict = {}
    name = ''


    for id in result["IdList"]:
        hand_two=Entrez.elink(dbfrom='assembly', db='nucleotide', id=id)
        result_two = Entrez.read(hand_two)
        logging.info('Assembly no: ' + str(id))
        logging.info(str(result_two[0]))
        ids_two = (result_two[0]["LinkSetDb"][1]["Link"])
        logging.info(str(result_two[0]["LinkSetDb"][1]["LinkName"]))
        logging.info(str(ids_two))


        name_id_dict = {}
        tax_id = ''
        for id_two in ids_two:

            handle_t = Entrez.esummary(db='nucleotide', id=id_two['Id'])
            result_t = Entrez.read(handle_t)
            handle_t.close()
            logging.info(str(result_t))

            tax_id = result_t[0]['TaxId']

            name = result_t[0]['Title']
            name = remove_spaces_and_dots(name)
            name_id_dict[name] = id_two['Id']

        logging.info(str(name_id_dict))
        logging.info(str())

        name = look_up_tax_name(tax_id)

        id_dict[name] = name_id_dict

    return id_dict

def get_genomes(dir, id_list, file_names_list):
    if not os.path.exists(dir):
        os.makedirs(dir)

    for i in range(len(id_list)):
        id = id_list[i]
        name=file_names_list[i]

        get_genome(id, dir, name)

def get_genomes_for_set_of_strains(id_dict, base_dir):

    for name in id_dict:
        print()
        print()
        print('----', name, '----')
        dir = base_dir + name + '/'
        list_of_names = []
        list_of_ids = []

        for gb_name in id_dict[name]:
            print(gb_name, ':', 'ID number = ', id_dict[name][gb_name])
            list_of_names.append(gb_name)
            list_of_ids.append(id_dict[name][gb_name])

        get_genomes(dir, list_of_ids, list_of_names)

def look_up_tax_name(tax_id):
    handle = Entrez.efetch(db='taxonomy', id=tax_id)
    result = Entrez.read(handle)

    name = result[0]['ScientificName']
    logging.info(str(name))

    name = remove_spaces_and_dots(name)

    return name

def get_genomes_for_org(dir, name, only_search_representative_genomes=True):
    dir = dir + name + '/'

    result = look_up_representative_assembly(name_of_org, only_search_representative_genomes)

    id_dict = get_ids_from_assembly_lookup_result(result)
    # id_dict has the format {'Org name' : {'Gb name' : id}}

    get_genomes_for_set_of_strains(id_dict, dir)

def remove_spaces_and_dots(string):

    new_string = ''
    for letter in string:

        if letter == ' ':
            letter = '_'

        if letter == '.':
            letter = ''

        if letter == ':':
            letter = ''

        if letter == '-':
            letter = ''

        new_string += letter

    return new_string





if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)

    dir = '/Users/Will/Documents/codonOpt_data/GenomeTmp/'
    representative = True

    name_of_org = 'Thermus thermophilus'

    get_genomes_for_org(dir, name_of_org, only_search_representative_genomes=representative)


