import Bio
import regex as re
import numpy as np
import pandas as pd
from Bio.Seq import translate, reverse_complement
from typing import NamedTuple
import matplotlib.pyplot as plt
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from genetic_data_file import Hidrophobic, bac_gencode, Nucleotides
from matplotlib.patches import ConnectionPatch
from data_generator import GeneticDataGenerator
from Bio.Data.CodonTable import TranslationError
# from Bio.SeqUtils import


def get_all_genes_type_and_amount(df: pd.DataFrame,
                                  col: str = 'type') -> dict:
    """
    Count genes by type (CDS/rRNA/misc_RNA etc...)
    :return: a dictionary with the key gene type and the value its amount
    ** 73612 line in gb file (misc_feature number 46) there is a gene that has CDS and misc_features
    """
    return dict(df[col].value_counts())


def calc_df_stats(_s: pd.Series) -> dict:
    return{
        'total_length': _s.sum(),
        'max_length': _s.max(),
        'min_length': _s.min(),
        'average_length': _s.mean(),
        'median_length': _s.median(),
    }


def calc_list_stats(_s: list) -> dict:
    return{
        'total_length': sum(_s),
        'max_length': max(_s),
        'min_length': min(_s),
        'average_length': np.average(_s),
        'median_length': np.median(_s),
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


def build_histograms(hist_title: str,
                     hist_value: list,
                     x_label: str,
                     y_label: str,
                     bins_num: int | list,
                     show_grid: bool = 'True',
                     graph_color: str = 'blue'):
    # x_low_lim: int = 0,
    # x_high_lim: int = 4000,
    # y_low_lim: int = 0,
    # y_high_lim: int = 200):
    """
    Builds histogram using Matplotlib
    :param hist_title: the title
    :param hist_value: the values for the histogram as list
    :param x_label: x axis label as str
    :param y_label: y axis label as str
    :param bins_num: number of bins in the graph as int
    :param show_grid: show grid as boolean
    :param graph_color: graph color as str
    """
    # :param x_low_lim: x axis minimum value as int
    # :param x_high_lim: x axis max value as int
    # :param y_low_lim: y axis minimum value as int
    # :param y_high_lim: y axis max value as int

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(hist_title)
    # plt.xlim(x_low_lim, x_high_lim)
    # plt.ylim(y_low_lim, y_high_lim)
    plt.grid(show_grid)
    plt.hist(x=hist_value, bins=bins_num, facecolor=graph_color)
    plt.tight_layout()
    plt.show()


# count occurrence of substrings' list in a sequence
def count_occ_in_seq(seq, occ):
    gc_sum = sum(map(lambda x: x in occ, seq))
    return gc_sum/len(seq)*100


def calculate_gc_percentage_in_genes(obj: GeneticDataGenerator,
                                     df_prot: pd.DataFrame) -> tuple[list:pd.DataFrame, list: float]:
    gc = ['G', 'C']
    full_gc_percent = count_occ_in_seq(obj.sequence.upper(), gc)
    df_prot = df_prot.assign(gc=df_prot.apply(lambda x: count_occ_in_seq(x.sequence, gc), axis=1))
    df_prot = df_prot.rename({'gc': 'gene gc%'}, axis=1)
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
# TODO : --------------------------  START HERE with Eyal --------------------------
class DataFileErrorGenes(NamedTuple):
    gene_name: str
    start: int
    end: int
    strand: int
    error: str


def convert_to_protein(gene_seq: str,
                       table: int,
                       strand: int) -> tuple[bool, str]:
    """
    convert gene sequence to protein
    :param gene_seq: the gene sequence as Bio.Seq
    :param strand: the strand 1/-1
    :param table: the table to use as int
    :return: a tuple[ error occurred false/true, the translation/error]
    """
    # print(table, strand, ind)
    if strand == -1:
        gene_seq = reverse_complement(gene_seq)

    try:
        return False, str(translate(gene_seq, table=table,
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


def consistent_checks_data_file(df: pd.DataFrame) -> list:
    """
    checks the data consistent in the data file -
    1. checks if all the genes that are translated to proteins length is a double of three times
    2. checks if the translation of the gene to proteins in the data file is correct
    3. Finds more errors using Bio translate function
    :return: a list with the genes that raised errors
    """

    data_file_errors = []

    # check genes numbers, e.g. if there are 4536 genes (spoiler: NOT) hahaha
    # genes_number = get_all_genes_type_and_amount(df)
    #
    # if genes_number['gene'] != (sum(genes_number.values()) - genes_number['gene'] - 1):
    #     gen_id = df.locus_tag.mode()
    #     # df_gene = df[df['locus_tag'] == gen_id[0]]
    #     # print(df_gene.head)
    #     print(gen_id)

    for ind in range(len(df.index)):
        if df.loc[ind, 'type'] == 'gene' and df.loc[ind+1, 'type'] == 'CDS':
            if df.loc[ind, 'start'] == df.loc[ind+1, 'start'] and df.loc[ind, 'end'] == df.loc[ind+1, 'end']:
                # print(df.loc[ind, 'transl_table'])
                error_found, gene_to_protein = convert_to_protein(gene_seq=df.loc[ind, 'sequence'],
                                                                  strand=df.loc[ind, 'strand'],
                                                                  table=df.loc[ind+1, 'transl_table'][0])
                if error_found:
                    error = gene_to_protein
                    data_file_errors.append(DataFileErrorGenes(gene_name=df.loc[ind, 'locus_tag'],
                                                               start=df.loc[ind, 'start'],
                                                               end=df.loc[ind, 'end'],
                                                               strand=df.loc[ind, 'strand'],
                                                               error=error))
                else:
                    if gene_to_protein != df.loc[ind+1, 'translation'][0]:
                        data_file_errors.append(DataFileErrorGenes(gene_name=df.loc[ind, 'locus_tag'],
                                                                   start=df.loc[ind, 'start'],
                                                                   end=df.loc[ind, 'end'],
                                                                   strand=df.loc[ind, 'strand'],
                                                                   error='Translations do not match'))
            else:
                pass

    print(data_file_errors)

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


# def uni_dataframe_details(df: pd.DataFrame,
#                           col: str):
#     """
#     Gathered information about the DF base on the request column
#     @param df:bla
#     @param col:bla
#     @return:asd
#     """
#
#     null_count = df[col].isna().sum()
#     name_df = df[df[col].notna()]
#     name_series = name_df[col]
#     name_unique = name_series.nunique()
#     name_count = len(list(name_series))
#     name_series.drop_duplicates(inplace=True)
#     total_gb = name_count + null_count
#     duplicate_names = name_count - name_unique
#
#     print('Total genes in UniProtKB file: ', total_gb)
#     print('Number of named genes: ', name_count)
#     print('Number of genes without name: ', null_count)
#     print('Number of unique genes name: ', name_unique)
#     print('Number of Duplicates genes by name: ', duplicate_names)
#
#     return name_df, list(name_series), name_unique, name_count, total_gb, duplicate_names


# def features_analysis(features: list):
#     """Analyze the features of the gene, mainly by genes name"""
#
#     features_prot_names = []
#     feature_prot = []
#     gb_null_counter = 0
#
#     for feature in features:
#         if feature.type == 'CDS':
#             feature_prot.append(feature)
#             if feature.qualifiers.get('gene') is not None:
#                 features_prot_names.append(feature.qualifiers.get('gene')[0])
#             else:
#                 gb_null_counter += 1
#
#     total_gb_prot = gb_null_counter + len(features_prot_names)
#
#     print('Total proteins genes in GeneBank file: ', total_gb_prot)
#     print('Number of named proteins in GB: ', len(features_prot_names))
#     print('Number of unnamed proteins in GB: ', gb_null_counter)
#
#     return feature_prot, features_prot_names, gb_null_counter, total_gb_prot


def compare_files_data(first_df: pd.DataFrame,
                       second_df: pd.DataFrame) -> tuple[list: pd.DataFrame]:
    # same protein in both
    same_prot = first_df[first_df['locus_tag'] == second_df['locus']]
    # protein exist only in genebank file
    gb_only = first_df[first_df['locus_tag'] != second_df['locus']]
    # protein exist only in UniPortKB file
    uni_only = second_df[second_df['locus'] != first_df['lucos_tag']]

    same_len = len(same_prot.index)
    gb_only_len = len(gb_only.index)
    uni_only_len = len(uni_only.index)
    same_dup = len(first_df.index) - same_len - gb_only_len
    print('Same proteins in both files: ', same_len)
    print('Proteins in GeneBank only: ', gb_only_len)
    print('Proteins in UniProt only: ', uni_only_len)
    print('Number of Duplicates in Same proteins list: ', same_dup)

    return same_prot, gb_only, uni_only


# def get_same_diff_from_gb(features: list,
#                           same_prot: list): pass
    # gb_same = []
    # gb_diff = []
    # gb_count = 0
    # for feature in features:
    #     if feature.type == 'CDS':
    #         if feature.qualifiers.get('gene') is not None:
    #             gb_count += 1
    #             if feature.qualifiers.get('gene')[0] in same_prot:
    #                 gb_same.append(feature)
    #             else:
    #                 gb_diff.append(feature)
    #
    # return gb_same, gb_diff, gb_count


# def create_diff_dataframe(df: pd.DataFrame,
#                           uni_only: list,
#                           col: str):
#     # Data frame of the differences, e.g. gene only in UniPort.
#     df_diff = df[df[col].isin(uni_only)]
#     return df_diff.drop_duplicates(subset=col)


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct, v=val)

    return my_autopct


def plot_gene_pie(_data, _labels, _titles, _exp, _angle, _width):
    fig = plt.figure(figsize=(15, 10))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    fig.subplots_adjust(wspace=0)

    # large pie chart parameters
    ax1.pie(_data[0], autopct=make_autopct(_data[0]), startangle=_angle,
            labels=_labels[0], explode=_exp, shadow=True)

    # small pie chart parameters
    ax2.pie(_data[1], autopct=make_autopct(_data[1]), startangle=_angle, textprops={'size': 'smaller'},
            labels=_labels[1], radius=0.5, shadow=True, explode=_exp)

    ax1.set_title(_titles[0])
    ax2.set_title(_titles[1])

    # use ConnectionPatch to draw lines between the two plots
    # get the wedge data
    theta1, theta2 = ax1.patches[0].theta1, ax1.patches[0].theta2
    center, r = ax1.patches[0].center, ax1.patches[0].r

    # draw top connecting line
    x = r * np.cos(np.pi / 180 * theta2) + center[0]
    y = np.sin(np.pi / 180 * theta2) + center[1]
    con = ConnectionPatch(xyA=(- _width / 2, .5), xyB=(x, y),
                          coordsA="data", coordsB="data", axesA=ax2, axesB=ax1)
    con.set_color([0, 0, 0])
    con.set_linewidth(2)
    ax2.add_artist(con)

    # draw bottom connecting line
    x = r * np.cos(np.pi / 180 * theta1) + center[0]
    y = np.sin(np.pi / 180 * theta1) + center[1]
    con = ConnectionPatch(xyA=(- _width / 2, -.5), xyB=(x, y), coordsA="data",
                          coordsB="data", axesA=ax2, axesB=ax1)
    con.set_color([0, 0, 0])
    ax2.add_artist(con)
    con.set_linewidth(2)

    plt.show()


# TODO: GET THIS FUNCTION WORKS!
# def hell_func_to_plot_graph():
#     uni_file_count = len(prot_names)
#     total_uni_file = uni_file_count + uni_null_count
#     uni_file_dup = uni_file_count - count_unique
#
#     uni_chart = [[uni_file_count, uni_null_count], [count_unique, uni_file_dup]]
#     labels = [['Named', 'Unnamed'], ['Unique', 'Duplicates']]
#     titles = ['Total Genes in UniPort', 'Named Genes']
#     explode = [0.1, 0]
#     angle = -180 * (uni_file_count / total_uni_file)
#     width = .2
#
#     plot_gene_pie(uni_chart, labels, titles, explode, angle, width)
#
#     # show GeneBank file details
#
#     total_named_prot_gb = len(features_prot)
#     total_unprot_gb = len(features) - total_gb_prot
#     total_gb = total_gb_prot + total_unprot_gb
#     gb_chart = [[total_gb_prot, total_unprot_gb], [total_named_prot_gb, gb_null_counter]]
#     labels = [['Protiens', 'Other Types'], ['Named', 'Unnamed']]
#     titles = ['Total Genes in GeneBank', 'Protiens Genes']
#     explode = [0.1, 0]
#     angle = -180 * (total_gb_prot / total_gb)
#     width = .2
#
#     plot_gene_pie(gb_chart, labels, titles, explode, angle, width)
#
#     compare_chart = [len(same_prot), len(gb_only), len(uni_only)]
#     _labels = ['Same Protiens', 'GB Unique Protiens', 'UniProtKB Unique Protiens']
#     plt.bar(_labels, compare_chart, color=['r', 'g', 'b'], width=0.5)
#     plt.title('Files Comparison')
#     plt.tight_layout()
#     plt.show()


trans_len = []
hidro_prec = []


def count_hidro(seq):
    hidro_prec.append(count_occ_in_seq(seq, Hidrophobic))


def find_all_pos(pos, seq):
    r1 = re.findall(r'\d*[0-9]\.\.\d*[0-9]', pos)  # find all coordinators of the seq
    res = []
    for r in r1:
        first, second = re.split(r'\.\.', r)  # split the coord and get the seq
        # first +=0 if first == 0 else first -= 1
        trans_seq = seq[int(first):int(second) + 1]
        count_hidro(trans_seq)
        res.append(trans_seq)
        trans_len.append(len(trans_seq))

    return res


def create_transmembrane_df(df: pd.DataFrame):
    trans_df = df[df['Transmembrane'].notna()]
    trans_df['Transmembrane Sequences'] = trans_df.apply(lambda x: find_all_pos(x.Transmembrane, x.Sequence), axis=1)
    return trans_df


def calc_gf_in_uni(df: pd.DataFrame,
                   new_col: str,
                   col: str,
                   sublist: list,
                   f):
    df[new_col] = df.apply(lambda x: f(x[col], sublist), axis=1)

    b_gc = list(df['GC%'])

    return df, max(b_gc), min(b_gc), np.median(b_gc), np.average(b_gc)


def count_mutation_by_type(_dict: dict = bac_gencode,
                           nuc: list = Nucleotides) -> dict:

    gen_dict = {}
    for key, value in _dict.items():
        synon = 0
        for codon in nuc:
            for position in range(3):
                if key[position] == codon:
                    continue
                else:
                    temp = list(key)
                    temp[position] = codon
                    temp = "".join(temp)
                    if value == _dict[temp]:
                        synon += 1
        gen_dict[key] = synon
    return gen_dict


def calc_selection(seq1: str,
                   seq2: str) -> tuple[list: int, str]:
    limit = 0.95
    seq1 = CodonSeq(seq1)
    seq2 = CodonSeq(seq2)
    dn, ds = cal_dn_ds(seq1, seq2)
    dn_ds_ratio = float(dn / ds)
    select = 'positive' if dn_ds_ratio >= 1 else 'neutral' if limit < dn_ds_ratio < 1 else 'negative'
    print("dN:%0.3f " % dn)
    print("dS:%0.3f " % ds)
    print("dN/dS:%0.3f " % dn_ds_ratio)
    print(f'The selection is: {select}')

    return dn, ds, dn_ds_ratio, select


calc_selection('atggtgctcagcgacgcagaatggcagttggtgctgaacatctgggcgaaggtggaagct',
               'atggggctcagcgacggggaatggcagttggtgctgaatgcctgggggaaggtggaggct')
