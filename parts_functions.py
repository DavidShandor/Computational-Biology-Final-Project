from typing import NamedTuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import regex as re
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import translate, reverse_complement
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from matplotlib.patches import ConnectionPatch
from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment
from automated_answer_file import AutomatedAnswerFile
from data_generator import GeneticDataGenerator
from genetic_data_file import Hidrophobic, bac_gencode, Nucleotides

trans_len = []
hidro_prec = []

# TODO: Document all and write answers into files.
#  NOTE: for each object we create there is another file (but maybe you can concat them?)


def get_all_genes_type_and_amount(df: pd.DataFrame,
                                  col: str = 'type',
                                  automated_answer_file: AutomatedAnswerFile = None) -> dict:
    """
    Count genes by type (CDS/rRNA/misc_RNA etc...)
    :return: a dictionary with the key gene type and the value its amount
    ** 73612 line in gb file (misc_feature number 46) there is a gene that has CDS and misc_features
    """
    automated_answer_file.write_answer_from_dict(answer_dict=dict(df[col].value_counts()),
                                                 section='1',
                                                 answer_description='All genes types and their count-',
                                                 question_number='1',
                                                 unwanted_fields_in_dict=['source'])
    return dict(df[col].value_counts())


def calc_df_stats(_s: pd.Series,
                  _description: str) -> dict:
    """
    @param _description:
    @param _s:
    @return:
    """

    return{
        'description': _description,
        'total': round(sum(_s), 2),
        'max': round(max(_s), 2),
        'min': round(min(_s), 2),
        'average': round(np.average(_s), 2),
        'median': round(np.median(_s), 2)
    }


def calc_list_stats(_s: list,
                    _description: str,
                    section: str = None,
                    automated_answer_file: AutomatedAnswerFile = None) -> dict:
    """
    @param automated_answer_file:
    @param section:
    @param _description:
    @param _s:
    @return:
    """
    stats = {
        'description': _description,
        'total': round(sum(_s), 2),
        'max': round(max(_s), 2),
        'min': round(min(_s), 2),
        'average': round(np.average(_s), 2),
        'median': round(np.median(_s), 2)
    }

    automated_answer_file.write_answer_from_dict(answer_dict=stats,
                                                 section=section,
                                                 data_description=_description,
                                                 unwanted_fields_in_dict=['description'])

    return stats


def characterization_of_gene_lengths(_df: pd.DataFrame,
                                     automated_answer_file: AutomatedAnswerFile = None) -> tuple[pd.DataFrame, list: dict, list: list]:
    """
    Counts the length and amount of protein/non-protein genes
    Find the maximum, minimum and average length of protein/non-protein genes
    Preserves the information about the lengths of the genes (all genes/protein/non-protein)
    :return: tuple with 2 dictionaries one with the info about the protein/ non-protein genes (their amount, length, etc.)
            3 lists with the information about the lengths of the protein/non-protein/all genes
    """
    # calc the length of the gene sequence
    _df['sequence length'] = _df['end'] - _df['start']

    # divide genes to proteins and non-proteins
    df_protein = _df[_df['type'] == 'CDS']
    df_non_protein = _df[~_df['type'].isin(['CDS', 'gene', 'source'])]

    protein_stats = calc_df_stats(df_protein['sequence length'], 'Proteins info')
    non_protein_stats = calc_df_stats(df_non_protein['sequence length'], 'non-Proteins info')

    automated_answer_file.write_answer_from_dict(answer_dict=protein_stats,
                                                 section='2',
                                                 answer_description='Statistics for protein and non proteins-',
                                                 data_description='Proteins info:',
                                                 unwanted_fields_in_dict=['description', 'total_length'])

    automated_answer_file.write_answer_from_dict(answer_dict=non_protein_stats,
                                                 data_description='Non proteins info:',
                                                 unwanted_fields_in_dict=['description', 'total_length'],
                                                 more_info='** histograms for all genes/protein/non protein will be shown on screen')

    return df_protein, (protein_stats, non_protein_stats),\
        (df_protein['sequence length'], df_non_protein['sequence length'], df_protein['sequence length'].append(df_non_protein['sequence length']))


def build_histograms(hist_title: str,
                     hist_value: list,
                     x_label: str,
                     y_label: str,
                     bins_num: int | list | range = None,
                     show_grid: bool = 'True',
                     graph_color: str = 'blue',
                     rows: int = None,
                     columns: int = None,
                     cell: int = None,
                     y_low_lim: int = 0,
                     y_high_lim: int = 0,
                     x_low_lim: int = 0,
                     x_high_lim: int = 0):

    """
    Builds histogram using Matplotlib
    @param y_high_lim:
    @param y_low_lim:
    @param hist_title: the title
    @param hist_value: the values for the histogram as list
    @param x_label: x axis label as str
    @param y_label: y axis label as str
    @param bins_num: number of bins in the graph as int
    @param show_grid: show grid as boolean
    @param graph_color: graph color as str
    @param rows:
    @param columns:
    @param cell:
    @param x_low_lim: x axis minimum value as int
    @param x_high_lim: x axis max value as int
    @param y_low_lim: y axis minimum value as int
    @param y_high_lim: y axis max value as int
    """
    if rows and columns:
        plt.subplot(rows, columns, cell)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(hist_title)

    if x_high_lim:
        plt.xlim(x_low_lim, x_high_lim)

    if y_high_lim:
        plt.ylim(y_low_lim, y_high_lim)

    plt.grid(show_grid)
    plt.hist(x=hist_value, bins=bins_num, facecolor=graph_color)

    if not(rows and columns):
        plt.tight_layout()
        plt.show()


# count occurrence of substrings' list in a sequence
def count_occ_in_seq(seq, occ):
    gc_sum = sum(map(lambda x: x in occ, seq))
    return gc_sum/len(seq)*100


def calculate_gc_percentage_in_genes(obj: GeneticDataGenerator,
                                     df_prot: pd.DataFrame,
                                     automated_answer_file: AutomatedAnswerFile = None) -> tuple[list:pd.DataFrame, list: float]:
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
    """

    @param automated_answer_file:
    @param obj:
    @param df_prot:
    @return:
    """
    gc = ['G', 'C']
    full_gc_percent = count_occ_in_seq(obj.sequence.upper(), gc)
    df_prot = df_prot.assign(gc=df_prot.apply(lambda x: count_occ_in_seq(x.sequence, gc), axis=1))
    df_prot = df_prot.rename({'gc': 'gene gc%'}, axis=1)
    avg_prot_gc_percent = df_prot['gene gc%'].mean()

    df_prot = df_prot.sort_values(by=['gene gc%'])
    largest_5 = df_prot.nlargest(5, 'gene gc%')
    smallest_5 = df_prot.nsmallest(5, 'gene gc%')

    answer_for_part_three_a = f'GC percentage in the whole genome: {full_gc_percent}'
    answer_for_part_three_b = f'GC average percentage in the proteins: {avg_prot_gc_percent}'

    automated_answer_file.write_answer_from_string(answer=answer_for_part_three_a,
                                                   section='3',
                                                   answer_description='Calculation of GC percentage in genes-')
    automated_answer_file.write_answer_from_string(answer=answer_for_part_three_b)
    automated_answer_file.write_answer_from_dataframe(answer_dataframe=largest_5,
                                                      data_description='Five highest GC percentage genes:',
                                                      specific_wanted_columns_in_dataframe=['gene', 'locus_tag', 'start', 'end', 'strand', 'gene gc%'])
    automated_answer_file.write_answer_from_dataframe(answer_dataframe=smallest_5,
                                                      data_description='Five lowest GC percentage genes:',
                                                      specific_wanted_columns_in_dataframe=['gene', 'locus_tag', 'start', 'end', 'strand', 'gene gc%'],
                                                      more_info='** histogram for the GC percentage in proteins will be shown on screen')
    return (df_prot, largest_5, smallest_5), (full_gc_percent, avg_prot_gc_percent)


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


def consistent_checks_data_file(df: pd.DataFrame,
                                automated_answer_file: AutomatedAnswerFile = None) -> list:
    """
    checks the data consistent in the data file -
    1. checks if all the genes that are translated to proteins length is a double of three times
    2. checks if the translation of the gene to proteins in the data file is correct
    3. Finds more errors using Bio translate function
    :return: a list with the genes that raised errors
    """

    data_file_errors = []

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

    if data_file_errors:
        for index, data in enumerate(data_file_errors):
            data_to_dict = {
                'gene_name': data[0],
                'start': data[1],
                'end': data[2],
                'strand': data[3],
                'error': data[4]
            }
            if index == 0:
                automated_answer_file.write_answer_from_dict(answer_dict=data_to_dict,
                                                             section='4',
                                                             answer_description='Errors found in data file (gene bank file)-',
                                                             data_description=f'Error number {index+1}:')
            else:
                automated_answer_file.write_answer_from_dict(answer_dict=data_to_dict,
                                                             data_description=f'Error number {index+1}:')

        create_csv_file(columns=['gene_name', 'start', 'end', 'strand', 'error'],
                        data=data_file_errors,
                        csv_file_name='gene_exceptions.csv')

    return data_file_errors


def create_csv_file(columns: list[str],
                    data: list[DataFileErrorGenes],
                    csv_file_name: str):
    """
    create a csv file
    :param columns: columns names as list
    :param data: the data to enter as list
    :param csv_file_name: the name of the file to create
    """
    all_data = extract_data_from_list_of_objects_to_list_of_lists(data)

    df = pd.DataFrame(data=all_data, columns=columns)
    df.to_csv(csv_file_name, index=False)


def fix_gb_uni_diff(gb_df: pd.DataFrame,
                    uni_df: pd.DataFrame,
                    col1: str, col2: str,
                    from_c: str, to_c: str) -> tuple[list: pd.DataFrame]:
    """
    Matches the different sequence identifiers and creates uniform identification.
    (GB Locus format: BSU_###, Uni: BSU###, need to remove "_")
    @param gb_df:
    @param uni_df:
    @param col1:
    @param col2:
    @param from_c:
    @param to_c:
    @return:
    """

    def adjust_string(value) -> str:
        """

        @param value: the string to fix
        @return: fixed string after remove irrelevant chars
        """
        if type(value) != str:
            return str(value).replace(from_c, to_c)
        return value.replace(from_c, to_c)

    gb_df[col1] = gb_df[col1].apply(lambda x: adjust_string(x))
    uni_df[col2] = uni_df[col2].apply(lambda x: adjust_string(x))

    return gb_df, uni_df


def compare_files_data(first_df: pd.DataFrame,
                       second_df: pd.DataFrame,
                       col1: str,
                       col2: str,
                       automated_answer_file: AutomatedAnswerFile = None) -> tuple[list: pd.DataFrame, list: int]:

    first_s = set(first_df[col1])
    second_s = set(second_df[col2])

    # same protein in both
    same_gene_list = first_s.intersection(second_s)
    same_gene_df = first_df[first_df[col1].isin(same_gene_list)]

    # protein exist only in genebank file
    first_only_list = first_s.difference(second_s)
    first_only = first_df[first_df[col1].isin(first_only_list)]

    # protein exist only in UniPortKB file
    second_only_list = second_s.difference(first_s)
    second_only = second_df[second_df[col2].isin(second_only_list)]

    same_len = len(same_gene_list)
    first_only_len = len(first_only_list)
    second_only_len = len(second_only_list)
    first_diff = len(first_df.index) - same_len - first_only_len
    second_diff = len(second_df.index) - same_len - second_only_len

    answer_for_two_a = (f'Same genes in both files: {same_len}\n'
                        f'\tProteins in first file only: {first_only_len}\n'
                        f'\tProteins in second file only: {second_only_len}\n'
                        f'\tFirst data file number of duplicates or unreviewed genes: {first_diff}\n'
                        f'\tSecond data file number of duplicates or unreviewed genes: {second_diff}')

    automated_answer_file.write_answer_from_string(answer=answer_for_two_a,
                                                   section='1',
                                                   answer_description='Protein comparison by gene locus tag-',
                                                   more_info='** pie chart for the difference between two databases files',
                                                   question_number='2')
    # print('Same genes in both files: ', same_len)
    # print('Genes in first file only: ', first_only_len)
    # print('Proteins in second file only: ', second_only_len)
    # print(f'First data file number of duplicates or unreviewed genes: {first_diff}\n'
    #       f'Second data file number of duplicates or unreviewed genes: {second_diff}')

    return same_gene_df, first_only, second_only,\
        (same_len, first_only_len, first_diff, second_only_len, second_diff)


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct, v=val)

    return my_autopct


def plot_pie(_data, _labels, _titles, _exp, _angle, _width):
    fig = plt.figure(figsize=(15, 10))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    fig.subplots_adjust(wspace=0)

    # large pie chart parameters
    ax1.pie(_data[0], autopct=make_autopct(_data[0]), startangle=_angle,
            explode=_exp[0], labels=_labels[0], shadow=True)

    # small pie chart parameters
    ax2.pie(_data[1], autopct=make_autopct(_data[1]), startangle=_angle, textprops={'size': 'smaller'},
            explode=_exp[1], labels=_labels[1], radius=0.5, shadow=True)

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


def compare_graph(len_list: list):

    same_len, first_only_len, first_diff, second_only_len, second_diff = len_list

    total_first = same_len+first_only_len+first_diff
    total_second = same_len+second_only_len+second_diff

    # First file

    labels = [['First File', 'Second File'], ['Same', 'First only']]
    title = [['Files Comparison', 'First in-Depth'], ['Files Comparison', 'Second in-Depth']]
    explode = [[0.1, 0], [0, 0.1]]
    data = [[total_first, total_second], [same_len, first_only_len]]

    if first_diff:
        labels = [['First File', 'Second File'], ['Same', 'First only', 'First Duplicates or Unreviewed']]
        explode = [[0.1, 0], [0, 0.1, 0.1]]
        data = [[total_first, total_second], [same_len, first_only_len, first_diff]]

    angle = -90 * ((same_len + first_only_len + first_diff) / total_first)
    plot_pie(_data=data, _labels=labels, _titles=title[0],
             _exp=explode, _angle=angle, _width=0.2)

    # Second file
    data = [[total_second, total_first], [same_len, second_only_len, second_diff]]
    labels = [['Second File', 'First File'], ['Same', 'Second only', 'Second Duplicates or Unreviewed']]
    explode = [[0.1, 0], [0, 0.1, 0.1]]
    angle = -180 * ((same_len + second_only_len + second_diff) / total_second)
    plot_pie(_data=data, _labels=labels, _titles=title[1],
             _exp=explode, _angle=angle, _width=0.2)

    compare_chart = [same_len, first_only_len, second_only_len]
    _labels = ['Same Proteins', 'First File Unique Proteins', 'Second File Unique Proteins']
    plt.bar(_labels, compare_chart, color=['r', 'g', 'b'], width=0.1)
    plt.title('Files Comparison')
    plt.tight_layout()
    plt.show()


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
    if ds != 0 and dn != 0:
        dn_ds_ratio = float(dn / ds)
    else:
        dn_ds_ratio = 1
    select = 'positive' if dn_ds_ratio > 1 else 'neutral' if limit < dn_ds_ratio <= 1 else 'negative'
    print("dN: %0.3f " % dn)
    print("dS: %0.3f " % ds)
    print("dN/dS: %0.3f " % dn_ds_ratio)
    print(f'The selection is: {select}')

    return dn, ds, dn_ds_ratio, select


def get_seq_by_prot(dna_seq, prot_seq):

    temp = ''
    for ind, ltr in enumerate(prot_seq):
        if ltr == '-':
            temp += '---'
        else:
            temp += dna_seq[ind*3:ind*3+3]

    if temp[:-3] in ['TAA', 'TAG', 'TGA']:
        temp = temp[:len(temp)-3]
    # modulo = len(temp) % 3
    # if modulo != 0:

    return temp


def protein_to_dnds(first_file_gene_seq: pd.DataFrame,
                    second_file_gene_seq: pd.DataFrame,
                    num: int = 5):

    first = first_file_gene_seq.sort_values(by=['gene'])
    second = second_file_gene_seq.sort_values(by=['gene'])

    seq_a = first['sequence'][:5]
    protein_a = first['translation'][:5].str[0]

    seq_b = second['sequence'][:5]
    protein_b = second['translation'][:5].str[0]

    for s_a, prot_a, s_b, prot_b in zip(seq_a, protein_a, seq_b, protein_b):
        alignments = align.localxx(prot_a, prot_b)
        align1, align2, score, begin, end = alignments[0]  # "unzip" the tuple
        print(format_alignment(align1=align1, align2=align2, score=score, begin=begin, end=end))
        seq1 = get_seq_by_prot(s_a, align1)
        seq2 = get_seq_by_prot(s_b, align2)
        print(seq1)
        print(seq2)
        print()
        calc_selection(seq1, seq2)





