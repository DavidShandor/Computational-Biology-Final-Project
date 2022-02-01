import matplotlib.pyplot as plt
import pandas as pd

from data_generator import GeneticDataGenerator
import parts_functions as func
from automated_answer_file import AutomatedAnswerFile
import warnings

warnings.simplefilter('ignore')

if __name__ == '__main__':
    print('Start project. Only histograms and graphs will be shown on the screen.\n'
          'All the data will be written into "automated_answer_file.txt" file.')
    automated_answer_file = AutomatedAnswerFile('automated_answer_file.txt')

    # ------------ File meta data--------------------
    # Index(['type', 'strand', 'start', 'end', 'sequence', 'db_xref', 'gene',
    #        'locus_tag', 'old_locus_tag', 'function', 'experiment', 'note',
    #        'codon_start', 'transl_table', 'product', 'protein_id', 'translation'],
    #       dtype='object')
    # type ={'gene': 4536, 'CDS': 4237, 'misc_RNA': 93,
    #       'misc_feature': 89, 'tRNA': 86, 'rRNA': 30,
    #       'ncRNA': 2, 'source': 1}

    # PART A
    # ----------------------------------

    gb_file = 'BS168.gb'
    cols = ['organism', 'mol_type', 'strain', 'sub_species', 'type_material',
            'pseudo', 'ncRNA_class', 'ribosomal_slippage', 'inference', 'EC_number']

    part_a = GeneticDataGenerator(genebank_file=gb_file,
                                  cols_to_drop=cols)

    all_gene_dict = func.get_all_genes_type_and_amount(df=part_a.gb_df,
                                                       automated_answer_file=automated_answer_file)

    df_proteins, stats, gen_len = func.characterization_of_gene_lengths(_df=part_a.gb_df,
                                                                        automated_answer_file=automated_answer_file)

    hist_title = ['Protein length', 'Non Protein length', 'All genes length']
    hist_val = gen_len
    colors = ['blue', 'green', 'yellow']
    for i, (_title, val, color) in enumerate(zip(hist_title, hist_val, colors), start=1):
        func.build_histograms(hist_title=_title,
                              hist_value=val,
                              x_label='Length',
                              y_label='Number of genes',
                              rows=1,
                              columns=3,
                              cell=i,
                              x_high_lim=17000,
                              y_high_lim=4500,
                              graph_color=color)

    plt.tight_layout()
    plt.show()

    df_prot, gc_percent = func.calculate_gc_percentage_in_genes(part_a, df_proteins,
                                                                automated_answer_file=automated_answer_file)

    bins = range(0, 100)
    func.build_histograms(hist_title='Gene GC percentage',
                          hist_value=list(df_prot[0].loc[:, 'gene gc%']),
                          x_label='GC percentage',
                          y_label='Number of genes',
                          bins_num=bins)

    col_to_5 = ['gene', 'start', 'end', 'strand', 'gene gc%', 'locus_tag', 'translation']

    data_file_errors = func.consistent_checks_data_file(part_a.gb_df,
                                                        automated_answer_file=automated_answer_file)

    # PART B
    # --------------- File meta data ---------------
    # # Index(['Entry', 'Entry name', 'Status', 'Protein names', 'Gene names',
    # #        'Length', 'Transmembrane', 'Sequence', 'name', 'locus'],
    # #       dtype='object')

    unifile = 'uniprot_file.xlsx'
    part_b = GeneticDataGenerator(unifile=unifile)

    part_b.uni_df.rename({'Gene names  (primary )': 'name',
                          'Gene names  (ordered locus )': 'locus'}, axis=1, inplace=True)

    fixed_gb_df, part_b.uni_df = func.fix_gb_uni_diff(df_prot[0], part_b.uni_df,
                                                      col1='locus_tag',
                                                      col2='locus',
                                                      from_c="_", to_c="")

    # len list = [same_len, first_only_len, first_diff, second_only_len, second_diff]
    same_prot_df, f_only_df, s_only_df, len_list = func.compare_files_data(first_df=fixed_gb_df,
                                                                           second_df=part_b.uni_df,
                                                                           col1='locus_tag',
                                                                           col2='locus',
                                                                           automated_answer_file=automated_answer_file)

    # show visualization of the data:
    func.compare_graph(len_list)

    trans_df = func.create_transmembrane_df(part_b.uni_df)

    trans_len = func.trans_len
    hidro_prec = func.hidro_prec
    trans_len_stats = func.calc_list_stats(trans_len,
                                           _description='Transmembrane Sequences Length Stats:',
                                           section='2',
                                           automated_answer_file=automated_answer_file)

    hidro_prec_stat = func.calc_list_stats(hidro_prec,
                                           _description='Transmembrane Sequences Hydrophobic Percentage Stats:',
                                           automated_answer_file=automated_answer_file)

    titles = ['Transmembrane Sequences Length Distribution', 'Transmembrane Sequences Hydro-acids% Distribution']
    val = [trans_len, hidro_prec]
    x_l = ['Transmembrane Sequences Length', 'Hydrophobic(%)']
    bins = [range(10, 40), range(0, 100)]

    for t, v, x, b in zip(titles, val, x_l, bins):
        func.build_histograms(hist_title=t,
                              hist_value=v,
                              x_label=x,
                              y_label='Number of Sequences',
                              bins_num=b)

    A_df = same_prot_df
    B_df = pd.merge(same_prot_df, trans_df, how='inner', left_on='locus_tag', right_on='locus')

    A_df_locus = set(A_df['locus_tag'])
    B_df_locus = set(B_df['locus_tag'])

    A_not_B_locus = A_df_locus.difference(B_df_locus)
    A_not_B_locus = A_df[A_df['locus_tag'].isin(A_not_B_locus)]

    A_gc_percentage = A_df['gene gc%']
    B_gc_percentage = B_df['gene gc%']
    A_not_B_locus_gc_percentage = A_not_B_locus['gene gc%']

    func.calc_list_stats(A_df['gene gc%'],
                         _description='Statistics for group A genes GC percentage',
                         section='3',
                         automated_answer_file=automated_answer_file)
    func.calc_list_stats(B_gc_percentage,
                         _description='Statistics for group B genes GC percentage',
                         automated_answer_file=automated_answer_file)
    func.calc_list_stats(A_not_B_locus_gc_percentage,
                         _description='Statistics for group A not B genes GC percentage',
                         automated_answer_file=automated_answer_file)

    histograms_titles = ['Group A', 'Group B', 'Group A not B']
    histograms_data = [A_gc_percentage, B_gc_percentage, A_not_B_locus_gc_percentage]
    colors = ['red', 'blue', 'green']
    bins = range(100)
    for index, (title, data, color) in enumerate(zip(histograms_titles, histograms_data, colors), start=1):
        func.build_histograms(hist_title=title,
                              hist_value=data,
                              x_label='Percentage',
                              y_label='Number of Genes',
                              bins_num=bins,
                              graph_color=color,
                              rows=2,
                              columns=2,
                              cell=index,
                              y_low_lim=0, y_high_lim=500)

    plt.subplot(2, 2, 4)
    plt.xlabel('Percentage ')
    plt.ylabel('Number of Genes')
    plt.title('B and A not B')
    plt.grid(True)
    plt.ylim(0, 500)
    plt.hist(x=[B_gc_percentage, A_not_B_locus_gc_percentage],
             bins=bins, color=['blue', 'yellow'], label=['B', 'A not B'])
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()

    # PART C
    # ------------------------------

    january = 'january2022.gb'
    july = 'july2020.gb'

    # list of cols to drop
    drop_jan = ['organism', 'mol_type', 'isolate', 'host', 'db_xref',
                'country', 'collection_date', 'ribosomal_slippage', 'codon_start']
    drop_jul = ['organism', 'mol_type', 'isolate', 'host', 'db_xref',
                'country', 'collection_date', 'ribosomal_slippage', 'codon_start',
                'inference', 'function', 'gene_synonym']

    january = GeneticDataGenerator(genebank_file=january,
                                   cols_to_drop=drop_jan, verbose=False)
    july = GeneticDataGenerator(genebank_file=july,
                                cols_to_drop=drop_jul, verbose=False)

    covid_synonyms = func.count_mutation_by_type(automated_answer_file=automated_answer_file)

    # len list = [same_len, first_only_len, first_diff, second_only_len, second_diff]
    same_gene_df, first_only_df, second_only_df, len_list = func.compare_files_data(first_df=january.gb_df,
                                                                                    second_df=july.gb_df,
                                                                                    col1='gene',
                                                                                    col2='gene',
                                                                                    automated_answer_file=automated_answer_file)
    # show comparison visualization
    func.compare_graph(len_list)

    # get 5 genes
    first_file_gene_seq = july.gb_df[july.gb_df['gene'].isin(same_gene_df['gene'])].dropna(subset='translation')
    second_file_gene_seq = january.gb_df[january.gb_df['gene'].isin(same_gene_df['gene'])].dropna(subset='translation')

    dnds_df = func.protein_to_dnds(first_file_gene_seq, second_file_gene_seq, num=12,
                                   automated_answer_file=automated_answer_file)
