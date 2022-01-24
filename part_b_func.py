import os
from Bio import SeqIO
import pandas as pd
import numpy as np
import regex as re
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch

table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

Hidrophobic = ('A', 'F', 'L', 'I', 'M', 'P', 'W')
Hidrophilic = ('G', 'S', 'Y', 'Q', 'K', 'D', 'T', 'C', 'N', 'R', 'H', 'E')


def read_uni_dataframe(filepath: str):
    """
    Read file and generate pandas dataframe
    @param filepath: file path to load
    @return: Data frame contain file data
    """
    assert (os.path.exists(filepath))  # sanity check
    return pd.read_excel(filepath)


def uni_dataframe_details(df: pd.DataFrame,
                          col: str):
    """
    Gathered information about the DF base on the request column
    @param df:bla
    @param col:bla
    @return:asd
    """

    null_count = df[col].isna().sum()
    name_df = df[df[col].notna()]
    name_series = name_df[col]
    name_unique = name_series.nunique()
    name_count = len(list(name_series))
    name_series.drop_duplicates(inplace=True)
    total_gb = name_count + null_count
    duplicate_names = name_count - name_unique

    print('Total genes in UniProtKB file: ', total_gb)
    print('Number of named genes: ', name_count)
    print('Number of genes without name: ', null_count)
    print('Number of unique genes name: ', name_unique)
    print('Number of Duplicates genes by name: ', duplicate_names)

    return name_df, list(name_series), name_unique, name_count, total_gb, duplicate_names


def parse_genebank_file(filepath: str):
    """Parse the GeneBank file and return its features"""
    assert (os.path.exists(filepath))

    with open(filepath, "r") as input_handle:
        gen = SeqIO.parse(input_handle, "genbank")
        record_gb = next(gen)

    return record_gb.features


def features_analysis(features: list):
    """Analyze the features of the gene, mainly by genes name"""

    features_prot_names = []
    feature_prot = []
    gb_null_counter = 0

    for feature in features:
        if feature.type == 'CDS':
            feature_prot.append(feature)
            if feature.qualifiers.get('gene') is not None:
                features_prot_names.append(feature.qualifiers.get('gene')[0])
            else:
                gb_null_counter += 1

    total_gb_prot = gb_null_counter + len(features_prot_names)

    print('Total proteins genes in GeneBank file: ', total_gb_prot)
    print('Number of named proteins in GB: ', len(features_prot_names))
    print('Number of unnamed proteins in GB: ', gb_null_counter)

    return feature_prot, features_prot_names, gb_null_counter, total_gb_prot


def compare_files_data(gb_name: list,
                       uni_name: list):
    # same protein in both
    same_prot = set(gb_name).intersection(set(uni_name))
    # protein exist only in genebank file
    gb_only = set(gb_name).difference(set(uni_name))
    # protein exist only in UniPortKB file
    uni_only = set(uni_name).difference(set(gb_name))

    # uni_only = pd.Series(list(uni_only))
    # uni_same = list(same_prot)

    same_len = len(same_prot)
    gb_only_len = len(gb_only)
    uni_only_len = len(uni_only)
    same_dup = len(gb_name) - same_len - gb_only_len
    print('Same proteins in both files: ', same_len)
    print('Proteins in GeneBank only: ', gb_only_len)
    print('Proteins in UniProt only: ', uni_only_len)
    print('Number of Duplicates in Same proteins list: ', same_dup)

    return same_prot, gb_only, uni_only, same_len, gb_only_len, uni_only_len, same_dup


def get_same_diff_from_gb(features: list,
                          same_prot: list):
    gb_same = []
    gb_diff = []
    gb_count = 0
    for feature in features:
        if feature.type == 'CDS':
            if feature.qualifiers.get('gene') is not None:
                gb_count += 1
                if feature.qualifiers.get('gene')[0] in same_prot:
                    gb_same.append(feature)
                else:
                    gb_diff.append(feature)

    return gb_same, gb_diff, gb_count


def create_diff_dataframe(df: pd.DataFrame,
                          uni_only: list,
                          col: str):
    # Data frame of the differences, e.g. gene only in UniPort.
    df_diff = df[df[col].isin(uni_only)]
    return df_diff.drop_duplicates(subset=col)


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


def count_occ_in_seq(seq, occ):
    c = sum(map(lambda x: x in occ, seq))
    return c / len(seq) * 100


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


def plot_hist(_data,
              _bins,
              xlabel: str,
              ylabel: str,
              title: str):
    plt.hist(_data, bins=_data)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, fontweight="bold")
    plt.tight_layout()
    plt.show()


def calc_gf_in_uni(df: pd.DataFrame,
                   new_col: str,
                   col: str,
                   sublist: list,
                   f):
    df[new_col] = df.apply(lambda x: f(x[col], sublist), axis=1)

    b_gc = list(df['GC%'])

    return df, max(b_gc), min(b_gc), np.median(b_gc), np.average(b_gc)


