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

class GatherInfoAboutAGenome:
    def __init__(self,
                 gene_bank_file: str,
                 answers_file: str):

        self.features = None
        self.record_gb = None
        self.answers_file = answers_file

        assert (os.path.exists(gene_bank_file)), 'gene bank file does not exist'

        with open(answers_file, 'w') as answersFile:
            answersFile.write('Part A answers\n')

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

        with open(self.answers_file, 'a') as answersFile:
            answersFile.write('\n1.Number of elements of every type in the file\n')
            for gene, amount in gene_type_and_amount.items():
                answersFile.write(f'\t"{gene}" : "{amount}"\n')

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
            'total_length': 0,
            'min_length': sys.maxsize,
            'max_length': -sys.maxsize - 1
        }

        non_protein_info = {
            'amount': 0,
            'total_length': 0,
            'min_length': sys.maxsize,
            'max_length': -sys.maxsize - 1
        }

        for feature in self.features:
            if not (feature.type == 'gene' or feature.type == 'source'):
                gen_length = feature.location.end.position - feature.location.start.position
                all_genes_lengths.append(gen_length)

                if feature.type == 'CDS':
                    protein_info['amount'] = protein_info['amount'] + 1
                    protein_info['total_length'] = protein_info['total_length'] + gen_length
                    if gen_length < protein_info['min_length']:
                        protein_info['min_length'] = gen_length
                    if gen_length > protein_info['max_length']:
                        protein_info['max_length'] = gen_length

                    protein_lengths.append(gen_length)
                else:
                    non_protein_info['amount'] = non_protein_info['amount'] + 1
                    non_protein_info['total_length'] = non_protein_info['total_length'] + gen_length
                    if gen_length < non_protein_info['min_length']:
                        non_protein_info['min_length'] = gen_length
                    if gen_length > non_protein_info['max_length']:
                        non_protein_info['max_length'] = gen_length

                    non_protein_lengths.append(gen_length)

        protein_info['average_length'] = protein_info['total_length'] / protein_info['amount']
        non_protein_info['average_length'] = non_protein_info['total_length'] / non_protein_info['amount']

        with open(self.answers_file, 'a') as answersFile:
            answersFile.write('\n2.Statistics for proteins/non proteins\n')
            answersFile.write('Proteins -\n')
            for key, value in protein_info.items():
                answersFile.write(f'\t"{key}" : "{value}"\n')
            answersFile.write('Non proteins -\n')
            for key, value in protein_info.items():
                answersFile.write(f'\t"{key}" : "{value}"\n')
            answersFile.write('** histograms for all genes/protein/non protein will be shown on screen\n')

        return protein_info, non_protein_info, protein_lengths, non_protein_lengths, all_genes_lengths

    def build_histograms(self,
                         histogram_title: str,
                         histogram_value: list,
                         x_label: str,
                         y_label: str,
                         x_low_lim: int,
                         x_high_lim: int,
                         y_low_lim: int,
                         y_high_lim: int,
                         bins_num: int,
                         show_grid: bool,
                         graph_color: str):
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

        # TODO add all_genes_data to description and add GeneInfo class the coverted gene (coverts to protein/ ...)
        proteins_gc_percentage = []
        five_highest_gc_genes = []
        five_lowest_gc_genes = []
        all_genes_data = []

        proteins_total_gc_counter = 0

        counter_g = self.record_gb.seq.upper().count('G')
        counter_c = self.record_gb.seq.upper().count('C')
        gc_percentage_full_gene = (counter_g + counter_c) / len(self.record_gb.seq)

        for feature in self.features:
            if feature.type == 'CDS':
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
                current_gene_seq = self.record_gb.seq[feature.location.start.position:
                                                      feature.location.end.position].upper()
                counter_g = current_gene_seq.count('G')
                counter_c = current_gene_seq.count('C')
                proteins_total_gc_counter = proteins_total_gc_counter + counter_g + counter_c

                current_gene_seq_gc_percentage = (counter_g + counter_c) / len(current_gene_seq)

                try:
                    gene_name = feature.qualifiers.get('gene')[0]
                except TypeError:
                    gene_name = ''

                if len(five_highest_gc_genes) < 5:
                    five_highest_gc_genes.append(GeneInfo(gc_percentage=current_gene_seq_gc_percentage,
                                                          gene_name=gene_name,
                                                          start=feature.location.start.position,
                                                          end=feature.location.end.position,
                                                          strand=feature.strand))
                if len(five_lowest_gc_genes) < 5:
                    five_lowest_gc_genes.append(GeneInfo(gc_percentage=current_gene_seq_gc_percentage,
                                                         gene_name=gene_name,
                                                         start=feature.location.start.position,
                                                         end=feature.location.end.position,
                                                         strand=feature.strand))

                if current_gene_seq_gc_percentage > lowest_gc_in_list(five_highest_gc_genes).gc_percentage:
                    five_highest_gc_genes = five_highest_gc(five_highest_gc_genes, feature,
                                                            current_gene_seq_gc_percentage)
                if current_gene_seq_gc_percentage < highest_gc_in_list(five_lowest_gc_genes).gc_percentage:
                    five_lowest_gc_genes = five_lowest_gc(five_lowest_gc_genes, feature, current_gene_seq_gc_percentage)

                all_genes_data.append(GeneInfo(gc_percentage=current_gene_seq_gc_percentage,
                                               gene_name=gene_name,
                                               start=feature.location.start.position,
                                               end=feature.location.end.position,
                                               strand=feature.strand))

        five_highest_gc_genes = sorted(five_highest_gc_genes, key=operator.attrgetter('gc_percentage'))
        five_lowest_gc_genes = sorted(five_lowest_gc_genes, key=operator.attrgetter('gc_percentage'))

        with open(self.answers_file, 'a') as answersFile:
            answersFile.write('\n3.Calculation of GC percentage in genes\n')
            answersFile.write(f'\t Full gene GC percentage: {gc_percentage_full_gene}\n')
            answersFile.write(f'\t Average GC percentage in proteins: {gc_average_percentage_proteins}\n')
            answersFile.write('Five highest GC percentage genes-\n')
            for gene in five_highest_gc_genes:
                answersFile.write(f'Gene name: {gene[1]}\n')
                answersFile.write(f'\t GC percentage: {gene[0]}\n')
                answersFile.write(f'\t Start: {gene[2]}\n')
                answersFile.write(f'\t End: {gene[3]}\n')
                answersFile.write(f'\t Strand: {gene[4]}\n')
            answersFile.write('Five lowest GC percentage genes-\n')
            for gene in five_lowest_gc_genes:
                answersFile.write(f'Gene name: {gene[1]}\n')
                answersFile.write(f'\t GC percentage: {gene[0]}\n')
                answersFile.write(f'\t Start: {gene[2]}\n')
                answersFile.write(f'\t End: {gene[3]}\n')
                answersFile.write(f'\t Strand: {gene[4]}\n')
            answersFile.write('** histogram for the GC percentage in proteins will be shown on screen\n')

        return gc_percentage_full_gene, gc_average_percentage_proteins, proteins_gc_percentage, five_highest_gc_genes, five_lowest_gc_genes, all_genes_data

    def consistent_checks_data_file(self):
        """
        checks the data consistent in the data file -
        1. checks if all the genes that are translated to proteins length is a double of three times
        2. checks if the translation of the gene to proteins in the data file is correct
        3. Finds more errors using Bio translate function
        :return: a list with the genes that raised errors
        """

        data_file_errors = []

        for index, feature in enumerate(self.features):
            cds_index = index + 1
            if feature.type == 'gene' and self.features[cds_index].type == 'CDS':
                if feature.location.start.position == self.features[cds_index].location.start.position and feature.location.end.position == self.features[cds_index].location.end.position:
                    error_found, gene_to_protein = convert_to_protein(gene_seq=self.record_gb.seq[feature.location.start.position:
                                                                                                  feature.location.end.position],
                                                                      strand=self.features[cds_index].location.strand,
                                                                      table=self.features[cds_index].qualifiers[
                                                                          'transl_table'][0])
                    if error_found:
                        error = gene_to_protein
                        data_file_errors.append(DataFileErrorGenes(gene_name=feature.qualifiers.get('gene')[0],
                                                                   start=feature.location.start.position,
                                                                   end=feature.location.end.position,
                                                                   strand=feature.strand,
                                                                   error=error))
                    else:
                        if gene_to_protein != self.features[cds_index].qualifiers['translation'][0]:
                            data_file_errors.append(DataFileErrorGenes(gene_name=feature.qualifiers.get('gene')[0],
                                                                       start=feature.location.start.position,
                                                                       end=feature.location.end.position,
                                                                       strand=feature.strand,
                                                                       error='Translations do not match'))
                else:
                    pass

        with open(self.answers_file, 'a') as answersFile:
            answersFile.write('\n4.Consistent check in the data file\n')
            answersFile.write('Data file errors -\n')
            for error in data_file_errors:
                answersFile.write(f'Gene name: {error[0]}\n')
                answersFile.write(f'\t Start: {error[1]}\n')
                answersFile.write(f'\t End: {error[2]}\n')
                answersFile.write(f'\t Strand: {error[3]}\n')
                answersFile.write(f'\t Error: {error[4]}\n')
            answersFile.write('** gene_exceptions.csv file with the errors was created\n')

        return data_file_errors

    def create_csv_file(self,
                        columns: list[str],
                        data: list[DataFileErrorGenes],
                        csv_file_name: str = None):
        """
        create a csv file
        :param columns: columns names as list
        :param data: the data to enter as list
        :param csv_file_name: the name of the file to create
        """
        all_data = extract_data_from_list_of_objects_to_list_of_lists(data)

        df = pd.DataFrame(data=all_data, columns=columns)
        df.to_csv(csv_file_name, index=False)


