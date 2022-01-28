import pandas as pd

from data_generator import GeneticDataGenerator
import parts_functions as func
from automated_answer_file import AutomatedAnswerFile
import warnings
warnings.simplefilter('ignore')


if __name__ == '__main__':
    automated_answer_file = AutomatedAnswerFile('automated_answer_file.txt')
    # TODO: doc all the files: main, func, data generator, etc., and write answers into files (check in the file order).
    print('*--------------- Part A ---------------*')

    # ------------ DO NOT DELETE THIS COMMENTS --------------------
    # Index(['type', 'strand', 'start', 'end', 'sequence', 'db_xref', 'gene',
    #        'locus_tag', 'old_locus_tag', 'function', 'experiment', 'note',
    #        'codon_start', 'transl_table', 'product', 'protein_id', 'translation'],
    #       dtype='object')
    # type ={'gene': 4536, 'CDS': 4237, 'misc_RNA': 93,
    #       'misc_feature': 89, 'tRNA': 86, 'rRNA': 30,
    #       'ncRNA': 2, 'source': 1}

    print('Question 1')

    gb_file = 'BS168.gb'
    cols = ['organism', 'mol_type', 'strain', 'sub_species', 'type_material',
            'pseudo', 'ncRNA_class', 'ribosomal_slippage', 'inference', 'EC_number']

    part_a = GeneticDataGenerator(genebank_file=gb_file,
                                  answers_file='Answers part A',
                                  cols_to_drop=cols)

    all_gene_dict = func.get_all_genes_type_and_amount(df=part_a.gb_df,
                                                       automated_answer_file=automated_answer_file)
    print(f'All genes and their count in file: {all_gene_dict}')

    print('\nQuestion 2')
    df_proteins, stats, gen_len = func.characterization_of_gene_lengths(_df=part_a.gb_df,
                                                                        automated_answer_file=automated_answer_file)

    print(f'Statistics of each group:')
    print(f'proteins: {stats[0]}\n'
          f'non proteins: {stats[1]}\n')

    print(f'The histograms are shown on screen')

    hist_title = ['Protein length', 'Non Protein length', 'All genes length']
    hist_val = gen_len
    colors = ['blue', 'green', 'yellow']
    bins = list(range(0, 100, 10))

    # TODO: check all histograms , bins are not good fit and also range. make sure its look good
    for _title, val, color in zip(hist_title, hist_val, colors):
        func.build_histograms(hist_title=_title,
                              hist_value=val,
                              x_label='Length',
                              y_label='Number of genes',
                              bins_num=bins,
                              graph_color=color)

    print('\nQuestion 3')

    df_prot, gc_percent = func.calculate_gc_percentage_in_genes(part_a, df_proteins,
                                                                automated_answer_file=automated_answer_file)

    print(f'GC percentage of the whole genome: {gc_percent[0]}')
    print(f'GC average percentage in the proteins: {gc_percent[1]}')

    print(f'The histogram is shown on screen')
    # TODO: check histogram
    func.build_histograms(hist_title='GC percentage proteins',
                          hist_value=list(df_prot[0].loc[:, 'gene gc%']),
                          x_label='GC percentage',
                          y_label='Number of proteins',
                          bins_num=bins)

    col_to_5 = ['gene', 'start', 'end', 'strand', 'gene gc%', 'locus_tag', 'translation']
    print(f'5 gene with the largest GC%:\n\t {df_prot[1].loc[:, col_to_5]}')
    print(f'5 gene with the smallest GC%:\n\t {df_prot[2].loc[:, col_to_5]}')

    print('\nQuestion 4')
    data_file_errors = func.consistent_checks_data_file(part_a.gb_df,
                                                        automated_answer_file=automated_answer_file)
    print(f'Data file errors: {data_file_errors}')

    # TODO: make this function works.
    # func.create_csv_file(columns=['Gene name', 'Start position', 'End position', 'Strand', 'Error'],
    #                      data=data_file_errors,
    #                      csv_file_name='gene_exceptions.csv')

    print('*--------------- Part B ---------------*')
    # --------------- once again: DO NOT DELTE THIS COMMENTS ---------------
    # # Index(['Entry', 'Entry name', 'Status', 'Protein names', 'Gene names',
    # #        'Length', 'Transmembrane', 'Sequence', 'name', 'locus'],
    # #       dtype='object')
    #
    unifile = 'uniprot_file.xlsx'
    part_b = GeneticDataGenerator(unifile=unifile, answers_file='Answer part B')

    part_b.uni_df.rename({'Gene names  (primary )': 'name',
                          'Gene names  (ordered locus )': 'locus'}, axis=1, inplace=True)

    fixed_gb_df, part_b.uni_df = func.fix_gb_uni_diff(df_prot[0], part_b.uni_df,
                                                      col1='locus_tag',
                                                      col2='locus',
                                                      from_c="_", to_c="")

    print('Question 1\n')
    # len list = [same_len, first_only_len, first_diff, second_only_len, second_diff]
    same_prot_df, f_only_df, s_only_df, len_list = func.compare_files_data(first_df=fixed_gb_df,
                                                                           second_df=part_b.uni_df,
                                                                           col1='locus_tag',
                                                                           col2='locus')

    # show visualization of the data:
    func.compare_graph(len_list)

    print('Question 2\n')
    trans_df = func.create_transmembrane_df(part_b.uni_df)

    trans_len = func.trans_len
    hidro_prec = func.hidro_prec
    # TODO: this paragraph need to be written into the answers file
    """
    The internal environment of the structure is Hydrophobic,
    and the proteins in the Trans-membranous move through this structure
    and hence the expectation that the sequences
    that create these proteins will be Hydrophobic or will have a high percentage of
    hydrophobic amino acids .This assumption is similar to the results we received in the histogram
    """
    trans_len_stats = func.calc_list_stats(trans_len, 'Trans-membranous stats')
    print('Transmembrane Sequences Length Stats: \n', trans_len_stats)
    hidro_prec_stat = func.calc_list_stats(hidro_prec, 'Hidrophobic stats')
    print('Transmembrane Sequences Hidrophobic Stats: \n', hidro_prec_stat)

    titles = ['Transmembrane Sequences Length Distribution', 'Transmembrane Sequences Hidro-acids% Distribution']
    val = [trans_len, hidro_prec]
    x_l = ['Transmembrane Sequences Length', 'Hidrophobic(%)']
    bins = [list(range(10, 40, 5)), list(range(0, 100, 10))]

    for t, v, x, b in zip(titles, val, x_l, bins):
        func.build_histograms(hist_title=t,
                              hist_value=v,
                              x_label=x,
                              y_label='Number of Sequences',
                              bins_num=b)

    print('Question 3:\n')
    # TODO: I start this question but didn't finish yet. It's just to draw histograms. leave it to me if you want.
    A_df = same_prot_df
    B_df = pd.merge(same_prot_df, trans_df, how='inner', left_on='locus_tag', right_on='locus')

    print(len(A_df), len(B_df), len(trans_df))

    print('#################### PART C ############################')
    print('Question 1\n')
    january = 'CoronaJanuary2022.gb'
    july = 'CoronaJuly2020.gb'

    # list of cols to drop
    drop_jan = ['organism', 'mol_type', 'isolate', 'host', 'db_xref',
                'country', 'collection_date', 'ribosomal_slippage', 'codon_start']
    drop_jul = ['organism', 'mol_type', 'isolate', 'host', 'db_xref',
                'country', 'collection_date', 'ribosomal_slippage', 'codon_start',
                'inference', 'function', 'gene_synonym']

    january = GeneticDataGenerator(genebank_file=january, answers_file='january',
                                   cols_to_drop=drop_jan, verbose=False)
    july = GeneticDataGenerator(genebank_file=july, answers_file='july',
                                cols_to_drop=drop_jul, verbose=False)

    print(january.gb_df['gene'].value_counts())
    print(july.gb_df['gene'].value_counts())

    covid_synon = func.count_mutation_by_type()
    print(covid_synon)

    print('Question 2\n')
    # len list = [same_len, first_only_len, first_diff, second_only_len, second_diff]
    same_gene_df, first_only_df, second_only_df, len_list = func.compare_files_data(first_df=january.gb_df,
                                                                                    second_df=july.gb_df,
                                                                                    col1='gene', col2='gene')

    # TODO: do the last part on this question.
