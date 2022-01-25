import operator
import os
import re
import sys
from typing import NamedTuple

import Bio
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from data_generator import GeneticDataGenerator


class GeneInfo(NamedTuple):
    gc_percentage: float
    gene_name: str
    start: int
    end: int
    strand: int


class DataFileErrorGenes(NamedTuple):
    gene_name: str
    start: int
    end: int
    strand: int
    error: str


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

    try:
        gene_name = feature.qualifiers.get('gene')[0]
    except TypeError:
        gene_name = ''

    five_highest_gc_genes.append(GeneInfo(gc_percentage=current_gene_seq_gc_percentage,
                                          gene_name=gene_name,
                                          start=feature.location.start.position,
                                          end=feature.location.end.position,
                                          strand=feature.strand))
    return five_highest_gc_genes


def five_lowest_gc(five_lowest_gc_genes: list,
                   feature: Bio.SeqFeature.SeqFeature,
                   current_gene_seq_gc_percentage: float) -> list:
    """
    gets a list of the lowest gc percentage genes finds the lowest one
    and replaces it with the new value
    :param five_lowest_gc_genes: a list with five genes (their info)
    :param feature: all the info about the gene
    :param current_gene_seq_gc_percentage: the gene gc percentage
    :return: a list with 5 genes
    """

    highest_value = highest_gc_in_list(five_lowest_gc_genes)
    five_lowest_gc_genes.remove(highest_value)

    try:
        gene_name = feature.qualifiers.get('gene')[0]
    except TypeError:
        gene_name = ''

    five_lowest_gc_genes.append(GeneInfo(gc_percentage=current_gene_seq_gc_percentage,
                                         gene_name=gene_name,
                                         start=feature.location.start.position,
                                         end=feature.location.end.position,
                                         strand=feature.strand))
    return five_lowest_gc_genes


def convert_to_protein(gene_seq: Bio.Seq.Seq,
                       strand: int,
                       table: int) -> tuple[bool, str]:
    """
    convert gene sequence to protein
    :param gene_seq: the gene sequence as Bio.Seq
    :param strand: the strand 1/-1
    :param table: the table to use as int
    :return: a tuple[ error occurred false/true, the translation/error]
    """
    if strand == -1:
        gene_seq = gene_seq.reverse_complement()

    try:
        return False, str(gene_seq.translate(table=table,
                                             cds=True))
    except TranslationError as error:
        return True, str(error)


def extract_data_from_list_of_objects_to_list_of_lists(list_with_data: list) -> list[list]:
    """
    makes a list with lists from a list with objects (function made for
    using pandas dataframe creator)
    :param list_with_data: the list to convert
    :return: list of lists
    """
    all_data = []

    for single_data_object in list_with_data:
        single_data_list = []
        for data in single_data_object:
            single_data_list.append(data)
        all_data.append(single_data_list)

    return all_data


def get_all_genes_type_and_amount(obj: GeneticDataGenerator) -> dict:
    """
    Count genes by type (CDS/rRNA/misc_RNA etc...)
    :return: a dictionary with the key gene type and the value its amount
    ** 73612 line in gb file (misc_feature number 46) there is a gene that has CDS and misc_features
    """
    s = obj.gb_df['type'].value_counts()

    return dict(s)


def calc_df_stats(_s: pd.Series) -> dict:
    return{
        'total length': _s.sum(),
        'max length': _s.max(),
        'min length': _s.min(),
        'average length': _s.mean(),
        'median length': _s.median(),
    }


def characterization_of_gene_lengths(obj: GeneticDataGenerator) -> \
        tuple[pd.DataFrame, dict, dict, list, list, list]:
    """
    Counts the length and amount of protein/non-protein genes
    Find the maximum, minimum and average length of protein/non-protein genes
    Preserves the information about the lengths of the genes (all genes/protein/non-protein)
    :return: tuple with 2 dictionaries one with the info about the protein/ non-protein genes (their amount, length, etc.)
            3 lists with the information about the lengths of the protein/non-protein/all genes
    """
    # calc the length of the gene sequence
    obj.gb_df['sequence length'] = obj.gb_df['end'] - obj.gb_df['start']

    # divide genes to proteins and non-proteins
    df_protein = obj.gb_df[obj.gb_df['type'] == 'CDS']
    df_non_protein = obj.gb_df[~obj.gb_df['type'].isin(['CDS', 'gene', 'source'])]

    protein_stats = calc_df_stats(df_protein['sequence length'])
    non_protein_stats = calc_df_stats(df_non_protein['sequence length'])

    # with open(self.answers_file, 'a') as answersFile:
    #     answersFile.write('\n2.Statistics for proteins/non proteins\n')
    #     answersFile.write('Proteins -\n')
    #     for key, value in protein_info.items():
    #         answersFile.write(f'\t"{key}" : "{value}"\n')
    #     answersFile.write('Non proteins -\n')
    #     for key, value in protein_info.items():
    #         answersFile.write(f'\t"{key}" : "{value}"\n')
    #     answersFile.write('** histograms for all genes/protein/non protein will be shown on screen\n')

    return df_protein, protein_stats, non_protein_stats, \
        list(df_protein['sequence length']), \
        list(df_non_protein['sequence length']),\
        list(obj.gb_df['sequence length'])


def build_histograms(histogram_title: str,
                     histogram_value: list,
                     x_label: str,
                     y_label: str,
                     bins_num: int = 200,
                     show_grid: bool = 'True',
                     graph_color: str = 'blue',
                     x_low_lim: int = 0,
                     x_high_lim: int = 4000,
                     y_low_lim: int = 0,
                     y_high_lim: int = 200):
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


# count occurrence of substrings' list in a sequence
def count_occ_in_seq(seq, occ):
    c = sum(map(lambda x: x in occ, seq))
    return c/len(seq)*100


def calculate_gc_percentage_in_genes(obj: GeneticDataGenerator,
                                     df_prot: pd.DataFrame) -> tuple[list: pd.DataFrame, list: float]:
    gc = ['G', 'C']
    full_gc_percent = count_occ_in_seq(obj.sequence.upper(), gc)

    df_prot['gene gc%'] = count_occ_in_seq(df_prot['sequence'].upper(), gc)
    avg_prot_gc_percent = df_prot['gene gc%'].mean()

    df_prot = df_prot.sort_values(by=['gene gc%'])
    largest_5 = df_prot.nlargest(5, 'gene gc%')
    smallest_5 = df_prot.nsmallest(5, 'gene gc%')

    return (df_prot, largest_5, smallest_5), (full_gc_percent, avg_prot_gc_percent)

    # """
    # 1. Calculates the GC percentage in the whole gene
    # 2. Calculates the average GC in the proteins
    # 3. Retains all GC percentage values of all proteins
    # 4. Finds the 5 highest GC percentage genes
    # 5. Finds the 5 lowest GC percentage genes
    # :param number_of_proteins: the number of the protein genes
    # :return: a tuple with -
    #         1. GC percentage in the whole gene as float
    #         2 average GC in the proteins as float
    #         3. list of all GC percentage values of all proteins
    #         4. 5 highest GC percentage genes as list
    #         5. 5 lowest GC percentage genes as list
    # """


    # # TODO add all_genes_data to description and add GeneInfo class the coverted gene (coverts to protein/ ...)

#     with open(self.answers_file, 'a') as answersFile:
#         answersFile.write('\n3.Calculation of GC percentage in genes\n')
#         answersFile.write(f'\t Full gene GC percentage: {gc_percentage_full_gene}\n')
#         answersFile.write(f'\t Average GC percentage in proteins: {gc_average_percentage_proteins}\n')
#         answersFile.write('Five highest GC percentage genes-\n')
#         for gene in five_highest_gc_genes:
#             answersFile.write(f'Gene name: {gene[1]}\n')
#             answersFile.write(f'\t GC percentage: {gene[0]}\n')
#             answersFile.write(f'\t Start: {gene[2]}\n')
#             answersFile.write(f'\t End: {gene[3]}\n')
#             answersFile.write(f'\t Strand: {gene[4]}\n')
#         answersFile.write('Five lowest GC percentage genes-\n')
#         for gene in five_lowest_gc_genes:
#             answersFile.write(f'Gene name: {gene[1]}\n')
#             answersFile.write(f'\t GC percentage: {gene[0]}\n')
#             answersFile.write(f'\t Start: {gene[2]}\n')
#             answersFile.write(f'\t End: {gene[3]}\n')
#             answersFile.write(f'\t Strand: {gene[4]}\n')
#         answersFile.write('** histogram for the GC percentage in proteins will be shown on screen\n')

#
# def consistent_checks_data_file(self):
#     """
#     checks the data consistent in the data file -
#     1. checks if all the genes that are translated to proteins length is a double of three times
#     2. checks if the translation of the gene to proteins in the data file is correct
#     3. Finds more errors using Bio translate function
#     :return: a list with the genes that raised errors
#     """
#
#     data_file_errors = []
#
#     for index, feature in enumerate(self.features):
#         cds_index = index + 1
#         if feature.type == 'gene' and self.features[cds_index].type == 'CDS':
#             if feature.location.start.position == self.features[cds_index].location.start.position and feature.location.end.position == self.features[cds_index].location.end.position:
#                 error_found, gene_to_protein = convert_to_protein(gene_seq=self.record_gb.seq[feature.location.start.position:
#                                                                                               feature.location.end.position],
#                                                                   strand=self.features[cds_index].location.strand,
#                                                                   table=self.features[cds_index].qualifiers[
#                                                                       'transl_table'][0])
#                 if error_found:
#                     error = gene_to_protein
#                     data_file_errors.append(DataFileErrorGenes(gene_name=feature.qualifiers.get('gene')[0],
#                                                                start=feature.location.start.position,
#                                                                end=feature.location.end.position,
#                                                                strand=feature.strand,
#                                                                error=error))
#                 else:
#                     if gene_to_protein != self.features[cds_index].qualifiers['translation'][0]:
#                         data_file_errors.append(DataFileErrorGenes(gene_name=feature.qualifiers.get('gene')[0],
#                                                                    start=feature.location.start.position,
#                                                                    end=feature.location.end.position,
#                                                                    strand=feature.strand,
#                                                                    error='Translations do not match'))
#             else:
#                 pass
#
#     with open(self.answers_file, 'a') as answersFile:
#         answersFile.write('\n4.Consistent check in the data file\n')
#         answersFile.write('Data file errors -\n')
#         for error in data_file_errors:
#             answersFile.write(f'Gene name: {error[0]}\n')
#             answersFile.write(f'\t Start: {error[1]}\n')
#             answersFile.write(f'\t End: {error[2]}\n')
#             answersFile.write(f'\t Strand: {error[3]}\n')
#             answersFile.write(f'\t Error: {error[4]}\n')
#         answersFile.write('** gene_exceptions.csv file with the errors was created\n')
#
#     return data_file_errors
#
#
# def create_csv_file(self,
#                     columns: list[str],
#                     data: list[DataFileErrorGenes],
#                     csv_file_name: str = None):
#     """
#     create a csv file
#     :param columns: columns names as list
#     :param data: the data to enter as list
#     :param csv_file_name: the name of the file to create
#     """
#     all_data = extract_data_from_list_of_objects_to_list_of_lists(data)
#
#     df = pd.DataFrame(data=all_data, columns=columns)
#     df.to_csv(csv_file_name, index=False)


