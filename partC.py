import os
import re
from Bio import SeqIO, pairwise2


def parser(files):
    """

    :param files: the genbank files we wish to analyse
    :return: the features extracted from the given genbank file
    """
    with open(files, "r") as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        record_gb = next(gen)  # content of 1st record
    features = record_gb.features
    print(len(features))
    print(features)
    return features


def get_dna_sequence(dna):
    """

    :param dna: the DNA sequence we wish to analyse (here its covid)
    :return: an STR representing the sequence
    """
    with open(dna, "r") as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        record_gb = next(gen)  # content of 1st record
    print(str(record_gb.seq))
    return str(record_gb.seq)


def features_compare(feat1, feat2):
    """

    :param feat1: first pair of the feature sequence we want to compare (july cov)
    :param feat2: second pair of the feature sequence we want to compare (january cov)
    :return: the similar features of these 2 sequences
    """
    for i, feature in enumerate(feat1[1:3]):
        print('i: ', i)
        print('feature.type: ', feature.type)
        print('feature.location: ', feature.location)
        print('feature.strand: ', feature.strand)
        print('feature.qualifiers: ', feature.qualifiers)
        print('feature.id: ', feature.id)
        print('=====================')


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

    # 3.2
    feat_jan = parser(january)
    feat_july = parser(july)
    features_compare(feat_jan, feat_july)


if __name__ == '__main__':
    main()
