import os
import re
import sys
from typing import NamedTuple

import Bio
import matplotlib.pyplot as plt
from Bio import SeqIO


class GeneInfo(NamedTuple):
    gc_percentage: int
    gene_name: str
    start: int
    end: int
    strand: int


def highest_gc_in_list(genes_info: list[GeneInfo]) -> GeneInfo:
    """
    finds the gene with the highest GC percentage
    :param genes_info: a list of genes with info about them
    :return: the highest GC percentage gene info in the list as GeneInfo
    """
    highest_gc_percentage = 0
    highest_gc_percentage_gene_index = 0

    for index, gene in enumerate(genes_info):
        if gene.gc_percentage > highest_gc_percentage:
            highest_gc_percentage = gene.gc_percentage
            highest_gc_percentage_gene_index = index

    return genes_info[highest_gc_percentage_gene_index]


def lowest_gc_in_list(genes_info: list[GeneInfo]) -> GeneInfo:
    """
    finds the gene with the lowest GC percentage
    :param genes_info: a list of genes with info about them
    :return: the lowest GC percentage gene info in the list as GeneInfo
    """
    lowest_gc_percentage = 100
    lowest_gc_percentage_gene_index = 0

    for index, gene in enumerate(genes_info):
        if gene.gc_percentage > lowest_gc_percentage:
            lowest_gc_percentage = gene.gc_percentage
            lowest_gc_percentage_gene_index = index

    return genes_info[lowest_gc_percentage_gene_index]


def five_highest_gc(five_highest_gc_genes: list,
                    feature: Bio.SeqFeature.SeqFeature,
                    current_gene_seq_gc_percentage: float) -> list:
    """
    gets a list of the highest gc percentage genes finds the lowest one
    and replaces it with the new value
    :param five_highest_gc_genes: a list with five genes (their info)
    :param feature: all the info about the gene
    :param current_gene_seq_gc_percentage: the gene gc percentage
    :return: a list with 5 genes
    """
    lowest_value = lowest_gc_in_list(five_highest_gc_genes)
    five_highest_gc_genes.remove(lowest_value)
    five_highest_gc_genes.append(GeneInfo(gc_percentage=current_gene_seq_gc_percentage,
                                          gene_name=feature.qualifiers.get('gene'),
                                          start=feature.location.start.position,
                                          end=feature.location.end.position,
                                          strand=feature.strand))
    return five_highest_gc_genes


def five_lowest_gc(five_lowest_gc_genes: list,
                   feature,
                   current_gene_seq_gc_percentage: float) -> list:
    highest_value = highest_gc_in_list(five_lowest_gc_genes)
    five_lowest_gc_genes.remove(highest_value)
    five_lowest_gc_genes.append(GeneInfo(gc_percentage=current_gene_seq_gc_percentage,
                                         gene_name=feature.qualifiers.get('gene'),
                                         start=feature.location.start.position,
                                         end=feature.location.end.position,
                                         strand=feature.strand))
    return five_lowest_gc_genes


class GatherInfoAboutAGenome:
    def __init__(self,
                 gene_bank_file: str):

        self.features = None
        self.record_gb = None

        assert (os.path.exists(gene_bank_file)), 'gene bank file does not exist'

        with open(gene_bank_file, "r") as input_handle:
            gen = SeqIO.parse(input_handle, "genbank")
            self.record_gb = next(gen)  # content of 1st record

        pattern = '[^ACGTacgt]'  # pattern of non-nucleotide
        test_string = str(self.record_gb.seq)
        result = re.match(pattern, test_string)
        if result is not None:
            print(f'There are letters that are not nucleotides: {result}')

        self.features = self.record_gb.features

    def get_all_genes_type_and_amount(self) -> dict:
        """
        Count genes by type (CDS/rRNA/misc_RNA etc...)
        :return: a dictionary with the key gene type and the value its amount
        ** 73612 line in gb file (misc_feature number 46) there is a gene that has CDS and misc_features
        """

        gene_type_and_amount = {}
        for feature in self.features:
            if feature.type == 'source':
                pass
            elif gene_type_and_amount.get(feature.type) is None:
                gene_type_and_amount[feature.type] = 1
            else:
                gene_type_and_amount[feature.type] = gene_type_and_amount[feature.type] + 1

        return gene_type_and_amount

    def characterization_of_gene_lengths(self) -> tuple[dict, dict, list, list, list]:
        """
        Counts the length and amount of protein/non protein genes
        Find the maximum, minimum and average length of protein/non protein genes
        Preserves the information about the lengths of the genes (all genes/protein/non protein)
        :return: tuple with 2 dictionaries one with the info about the protein/ non protein genes (their amount, length, etc..)
                3 lists with the information about the lengths of the protein/non protein/all genes
        """

        protein_lengths = []
        non_protein_lengths = []
        all_genes_lengths = []

        protein_info = {
            'amount': 0,
            'length': 0,
            'min_length': sys.maxsize,
            'max_length': -sys.maxsize - 1
        }

        non_protein_info = {
            'amount': 0,
            'length': 0,
            'min_length': sys.maxsize,
            'max_length': -sys.maxsize - 1
        }

        for feature in self.features:
            if not(feature.type == 'gene' or feature.type == 'source'):
                gen_length = feature.location.end.position - feature.location.start.position
                all_genes_lengths.append(gen_length)

                if feature.type == 'CDS':
                    protein_info['amount'] = protein_info['amount'] + 1
                    protein_info['length'] = protein_info['length'] + gen_length
                    if gen_length < protein_info['min_length']:
                        protein_info['min_length'] = gen_length
                    if gen_length > protein_info['max_length']:
                        protein_info['max_length'] = gen_length

                    protein_lengths.append(gen_length)
                else:
                    non_protein_info['amount'] = non_protein_info['amount'] + 1
                    non_protein_info['length'] = non_protein_info['length'] + gen_length
                    if gen_length < non_protein_info['min_length']:
                        non_protein_info['min_length'] = gen_length
                    if gen_length > non_protein_info['max_length']:
                        non_protein_info['max_length'] = gen_length

                    non_protein_lengths.append(gen_length)

        protein_info['average_length'] = protein_info['length'] / protein_info['amount']
        non_protein_info['average_length'] = non_protein_info['length'] / non_protein_info['amount']

        return protein_info, non_protein_info, protein_lengths, non_protein_lengths, all_genes_lengths

    def build_histograms(self,
                         histogram_title: str,
                         histogram_value: list,
                         x_label: str = 'Length',
                         y_label: str = 'Number of genes',
                         x_low_lim: int = 0,
                         x_high_lim: int = 4000,
                         y_low_lim: int = 0,
                         y_high_lim: int = 200,
                         bins_num: int = 300,
                         show_grid: bool = True,
                         graph_color: str = 'blue'):
        """
        Builds histogram using Matplotlib
        :param histogram_title: the title
        :param histogram_value: the values for the histogram as list
        :param x_label: x axis label as str
        :param y_label: y axis label as str
        :param x_low_lim: x axis minimum value as int
        :param x_high_lim: x axis max value as int
        :param y_low_lim: y axis minimum value as int
        :param y_high_lim: y axis max value as int
        :param bins_num: number of bins in the graph as int
        :param show_grid: show grid as boolean
        :param graph_color: graph color as str
        """

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(histogram_title)
        plt.xlim(x_low_lim, x_high_lim)
        plt.ylim(y_low_lim, y_high_lim)
        plt.grid(show_grid)
        plt.hist(x=histogram_value, bins=bins_num, facecolor=graph_color)
        plt.show()

    def calculate_gc_percentage_in_genes(self,
                                         number_of_proteins: int) -> tuple[float, float, list, list, list]:
        """
        1. Calculates the GC percentage in the whole gene
        2. Calculates the average GC in the proteins
        3. Retains all GC percentage values of all proteins
        4. Finds the 5 highest GC percentage genes
        5. Finds the 5 lowest GC percentage genes
        :param number_of_proteins: the number of the protein genes
        :return: a tuple with -
                1. GC percentage in the whole gene as float
                2 average GC in the proteins as float
                3. list of all GC percentage values of all proteins
                4. 5 highest GC percentage genes as list
                5. 5 lowest GC percentage genes as list
        """

        # TODO: check if all letters are valid (G/C/T/A)
        proteins_gc_percentage = []
        five_highest_gc_genes = []
        five_lowest_gc_genes = []

        proteins_total_gc_counter = 0

        counter_g = self.record_gb.seq.upper().count('G')
        counter_c = self.record_gb.seq.upper().count('C')
        gc_percentage_full_gene = (counter_g + counter_c) / len(self.record_gb.seq)

        for feature in self.features:
            if feature.type == 'CDS':
                # TODO: check how to read the seq (start or start -1)
                current_gene_seq = self.record_gb.seq[feature.location.start.position:
                                                      feature.location.end.position].upper()
                counter_g = current_gene_seq.count('G')
                counter_c = current_gene_seq.count('C')
                proteins_total_gc_counter = proteins_total_gc_counter + counter_g + counter_c

                single_protein_gc_percentage = (counter_g + counter_c) / len(current_gene_seq)
                proteins_gc_percentage.append(single_protein_gc_percentage)

        gc_average_percentage_proteins = proteins_total_gc_counter / number_of_proteins

        for feature in self.features:
            if feature.type == 'gene':
                # TODO: check how to read the seq (start or start -1)
                current_gene_seq = self.record_gb.seq[feature.location.start.position:
                                                      feature.location.end.position].upper()
                counter_g = current_gene_seq.count('G')
                counter_c = current_gene_seq.count('C')
                proteins_total_gc_counter = proteins_total_gc_counter + counter_g + counter_c

                current_gene_seq_gc_percentage = (counter_g + counter_c) / len(current_gene_seq)

                if len(five_highest_gc_genes) < 5:
                    five_highest_gc_genes.append(GeneInfo(gc_percentage=current_gene_seq_gc_percentage,
                                                          gene_name=feature.qualifiers.get('gene'),
                                                          start=feature.location.start.position,
                                                          end=feature.location.end.position,
                                                          strand=feature.strand))
                if len(five_lowest_gc_genes) < 5:
                    five_lowest_gc_genes.append(GeneInfo(gc_percentage=current_gene_seq_gc_percentage,
                                                         gene_name=feature.qualifiers.get('gene'),
                                                         start=feature.location.start.position,
                                                         end=feature.location.end.position,
                                                         strand=feature.strand))

                if current_gene_seq_gc_percentage > lowest_gc_in_list(five_highest_gc_genes).gc_percentage:
                    five_highest_gc_genes = five_highest_gc(five_highest_gc_genes, feature, current_gene_seq_gc_percentage)
                if current_gene_seq_gc_percentage < highest_gc_in_list(five_lowest_gc_genes).gc_percentage:
                    five_lowest_gc_genes = five_lowest_gc(five_lowest_gc_genes, feature, current_gene_seq_gc_percentage)


        return gc_percentage_full_gene, gc_average_percentage_proteins, proteins_gc_percentage, five_highest_gc_genes, five_lowest_gc_genes
