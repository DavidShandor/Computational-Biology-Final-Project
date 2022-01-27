import os
import re
from Bio import SeqIO, pairwise2

gencode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

bac_gencode = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '_', 'TGA': '_',
    'TTG': 'L', 'TCG': 'S', 'TAG': '_', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}


def parser(files):
    """

    :param files: the genbank files we wish to analyse
    :return: the features extracted from the given genbank file
    """
    with open(files, "r") as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        record_gb = next(gen)  # content of 1st record
    features = record_gb.features
    return features


def get_dna_sequence(dna):
    """

    :param dna: the DNA sequence we wish to analyse (here its covid)
    :return: an STR representing the sequence
    """
    with open(dna, "r") as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        record_gb = next(gen)  # content of 1st record
    print("total features: " + str(len(record_gb.features)))
    count = 0
    for feature in record_gb.features:
        if feature.type == 'CDS':
            count += 1
            # print(feature)
    print("number of CDS features: " + str(count))
    return str(record_gb.seq)


def features_compare(feat1, feat2):
    """

    :param feat1: first pair of the feature sequence we want to compare (july cov)
    :param feat2: second pair of the feature sequence we want to compare (january cov)
    :return: the similar features of these 2 sequences
    """
    feat1_list = []
    feat2_list = []
    similar = []
    dif = []
    feat1_count = count_genes(feat1, feat1_list)
    feat2_count = count_genes(feat2, feat2_list)
    print(feat1_count, feat2_count)

    for feature in feat1_list:
        if feature in feat2_list:
            similar.append(feature)
        else:
            dif.append(feature)
    print("similar genes:")
    print(similar)
    print("different genes:")
    print(dif)


def count_genes(features, feat_list):
    print("These are the features:")
    print("================================== \n")
    num_genes = 0
    for f in features:
        if "gene" in f.qualifiers.keys():
            name = f.qualifiers["gene"][0]
            print(f.qualifiers["gene"][0])
            feat_list.append(name)
    #for i in range(len(features)):
        #f = features[i]
        #if f.type == 'gene' or f.type == 'CDS':
            # feat_list = f.qualifiers['locus_tag']
            #print(f.qualifiers)
            # feat_list.append(f.qualifiers)
            #num_genes += 1
    return feat_list


def count_codons(seq):
    """

    :param seq: this is the DNA sequence we want to analyse
    :return: a Dictionary with the amounts of every CODON and how many times it appeared in the DNA sequence
    """
    codon = {}
    protein = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    codon = dict.fromkeys(protein)
    for triplet in protein:
        count = protein.count(triplet)
        codon[triplet] = count
    print(codon)
    return codon


def count_mutation_by_type():
    codons = ('A', 'C', 'T', 'G')
    gen_dict = {}
    for key, value in gencode.items():
        synon = 0
        for codon in codons:
            for position in range(3):
                if key[position] == codon:
                    continue
                else:
                    temp = list(key)
                    temp[position] = codon
                    temp = "".join(temp)
                    if value == gencode[temp]:
                        synon += 1
        gen_dict[key] = synon
    return gen_dict


def main():
    # PART III
    ## read file and import details
    january = 'CoronaJanuary2022.gb'
    july = 'CoronaJuly2020.gb'
    # making sure that the path is valid
    assert (os.path.exists(january))
    assert (os.path.exists(july))

    # 3.1

    corona_dna_sequence = get_dna_sequence(july)
    count_codons(corona_dna_sequence)
    print(count_mutation_by_type())
    print("END OF PART 1")
    print("================================== \n")

    # 3.2
    feat_jan = parser(january)
    feat_july = parser(july)
    features_compare(feat_jan, feat_july)


if __name__ == '__main__':
    main()
