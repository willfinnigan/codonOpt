import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__)) # This is your Project Root

DNA_to_aa = {
    "GCT" :"A", "GCG" :"A", "GCC" :"A", "GCA" :"A",
    "AGA" :"R", "CGA" :"R", "CGT" :"R", "AGG" :"R", "CGC" :"R", "CGG" :"R",
    "AAT" :"N", "AAC" :"N",
    "GAT" :"D", "GAC" :"D",
    "TGT" :"C", "TGC" :"C",
    "TAA" :"*", "TAG" :"*", "TGA" :"*",
    "CAA" :"Q", "CAG" :"Q",
    "GAA" :"E", "GAG" :"E",
    "GGT" :"G", "GGA" :"G", "GGC" :"G", "GGG" :"G",
    "CAT" :"H", "CAC" :"H",
    "ATA" :"I", "ATT" :"I", "ATC" :"I",
    "TTA" :"L", "TTG" :"L", "CTT" :"L", "CTG" :"L", "CTC" :"L", "CTA" :"L",
    "AAA" :"K", "AAG" :"K",
    "ATG" :"M",
    "TTT" :"F", "TTC" :"F",
    "CCA" :"P", "CCT" :"P", "CCG" :"P", "CCC" :"P",
    "TCA" :"S", "AGT" :"S", "TCT" :"S", "TCG" :"S", "AGC" :"S", "TCC" :"S",
    "ACA" :"T", "ACT" :"T", "ACG" :"T", "ACC" :"T",
    "TGG" :"W",
    "TAT" :"Y", "TAC" :"Y",
    "GTA" :"V", "GTT" :"V", "GTC" :"V", "GTG" :"V"}


amino_acids_list = ['V', 'G', 'F', 'E', 'N', 'P', 'Q', 'M', 'K', 'T', 'S', 'W', 'A', 'R', 'D', 'L', 'Y', 'H', 'I', 'C', '*']

empty_codon_table = {'E' : {'GAA': 0, 'GAG': 0},
                    'V' : {'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0},
                    'R' : {'CGC': 0, 'AGG': 0, 'CGT': 0, 'CGG': 0, 'CGA': 0, 'AGA': 0},
                    'M' : {'ATG': 0},
                    'T' : {'ACG': 0, 'ACT': 0, 'ACA': 0, 'ACC': 0},
                    'G' : {'GGC': 0, 'GGT': 0, 'GGA': 0, 'GGG': 0},
                    'C' : {'TGT': 0, 'TGC': 0},
                    'A' : {'GCC': 0, 'GCG': 0, 'GCA': 0, 'GCT': 0},
                    'Y' : {'TAT': 0, 'TAC': 0},
                    'P' : {'CCC': 0, 'CCG': 0, 'CCA': 0, 'CCT': 0},
                    'N' : {'AAT': 0, 'AAC': 0},
                    'W' : {'TGG': 0},
                    'H' : {'CAT': 0, 'CAC': 0},
                    'D' : {'GAC': 0, 'GAT': 0},
                    'Q' : {'CAG': 0, 'CAA': 0},
                    'F' : {'TTC': 0, 'TTT': 0},
                    '*' : {'TAA': 0, 'TGA': 0, 'TAG': 0},
                    'L' : {'CTA': 0, 'CTC': 0, 'CTG': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0},
                    'S' : {'TCT': 0, 'TCC': 0, 'AGC': 0, 'TCA': 0, 'AGT': 0, 'TCG': 0},
                    'I' : {'ATA': 0, 'ATC': 0, 'ATT': 0},
                    'K' : {'AAG': 0, 'AAA': 0}}

amino_acids_dict = {"A" : 0,
            "R" : 0,
            "N" : 0,
            "D" : 0,
            "C" : 0,
            "*" : 0,
            "Q" : 0,
            "E" : 0,
            "G"	: 0,
            "H" : 0,
            "I" : 0,
            "L" : 0,
            "K" : 0,
            "M" : 0,
            "F" : 0,
            "P" : 0,
            "S"	: 0,
            "T" : 0,
            "W" : 0,
            "Y" : 0,
            "V" : 0}

triplets_only_table = total = {"GCT" : 0, "GCG" : 0, "GCC" : 0, "GCA": 0,
            "AGA" : 0, "CGA" : 0, "CGT" : 0, "AGG" : 0, "CGC" : 0, "CGG" : 0,
            "AAT" : 0, "AAC" : 0,
            "GAT" : 0, "GAC" : 0,
            "TGT" : 0, "TGC" : 0,
            "TAA" : 0, "TAG" : 0, "TGA" : 0,
            "CAA" : 0, "CAG" : 0,
            "GAA" : 0, "GAG" : 0,
            "GGT" : 0, "GGA" : 0, "GGC" : 0, "GGG" : 0,
            "CAT" : 0, "CAC" : 0,
            "ATA" : 0, "ATT" : 0, "ATC" : 0,
            "TTA" : 0, "CTA" : 0, "TTG" : 0, "CTT" : 0, "CTG" : 0, "CTC" : 0,
            "AAA" : 0, "AAG" : 0,
            "ATG" : 0,
            "TTT" : 0, "TTC" : 0,
            "CCA" : 0, "CCT" : 0, "CCG" : 0, "CCC" : 0,
            "TCA" : 0, "AGT" : 0, "TCT" : 0, "TCG" : 0, "AGC" : 0, "TCC" : 0,
            "ACA" : 0, "ACT" : 0, "ACG" : 0, "ACC" : 0,
            "TGG" : 0,
            "TAT" : 0, "TAC" : 0,
            "GTA" : 0, "GTT" : 0, "GTC" : 0, "GTG" : 0}




