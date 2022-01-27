from data_generator import GeneticDataGenerator
import parts_functions as func

import warnings
warnings.simplefilter('ignore')


if __name__ == '__main__':

    # print('*--------------- Part A ---------------*')
    # Index(['type', 'strand', 'start', 'end', 'sequence', 'db_xref', 'gene',
    #        'locus_tag', 'old_locus_tag', 'function', 'experiment', 'note',
    #        'codon_start', 'transl_table', 'product', 'protein_id', 'translation'],
    #       dtype='object')
    # type ={'gene': 4536, 'CDS': 4237, 'misc_RNA': 93,
    #       'misc_feature': 89, 'tRNA': 86, 'rRNA': 30,
    #       'ncRNA': 2, 'source': 1}

    # print('Question 1')
    #
    gb_file = 'BS168.gb'
    cols = ['organism', 'mol_type', 'strain', 'sub_species', 'type_material',
            'pseudo', 'ncRNA_class', 'ribosomal_slippage', 'inference', 'EC_number']
    part_a = GeneticDataGenerator(genebank_file=gb_file,
                                  answers_file='Answers part A',
                                  cols_to_drop=cols)

    # print(part_a.gb_df.index)
    # func.consistent_checks_data_file(part_a.gb_df)
    #
    # all_gene_dict = func.get_all_genes_type_and_amount(part_a.gb_df)
    # print(f'All genes and their count in file: {all_gene_dict}')
    # #
    # print('\nQuestion 2')
    # df_proteins, protein_info, non_protein_info, protein_lengths, non_protein_lengths,\
    #     all_genes_lengths = func.characterization_of_gene_lengths(part_a)
    # #
    # print(f'Statistics of each group:')
    # print(f'proteins: {protein_info}\n'
    #       f'non proteins: {non_protein_info}\n')
    #
    # print(f'The histograms are shown on screen')
    #
    # hist_title = ['Non Protein length', 'Protein length', 'All genes length']
    # hist_val = [non_protein_lengths, protein_lengths, all_genes_lengths]
    # colors = ['blue', 'green', 'yellow']
    # bins = list(range(0, 100, 10))  # TODO: check this out, bins are not good fit
    #
    # for _title, val, color in zip(hist_title, hist_val, colors):
    #     func.build_histograms(hist_title=_title,
    #                           hist_value=val,
    #                           x_label='Length',
    #                           y_label='Number of genes',
    #                           bins_num=bins,
    #                           graph_color=color)
    #
    # print('\nQuestion 3')
    #
    # df_prot, gc_percent = func.calculate_gc_percentage_in_genes(part_a, df_proteins)
    #
    # print(f'GC percentage of the whole genome: {gc_percent[0]}')
    # print(f'GC average percentage in the proteins: {gc_percent[1]}')
    #
    # print(f'The histogram is shown on screen')
    # func.build_histograms(hist_title='GC percentage proteins',
    #                       hist_value=list(df_prot[0].loc[:, 'gene gc%']),
    #                       x_label='GC percentage',
    #                       y_label='Number of proteins',
    #                       bins_num=bins)
    #
    # col_to_5 = ['gene', 'start', 'end', 'strand', 'locus_tag', 'translation']
    # print(f'5 gene with the largest GC%:\n\t {df_prot[1].loc[:, col_to_5]}')
    # print(f'5 gene with the smallest GC%:\n\t {df_prot[2].loc[:, col_to_5]}')
    #
    # print('\nQuestion 4')
    # # data_file_errors = func.consistent_checks_data_file()
    # # print(f'Data file errors: {data_file_errors}')
    # #
    # # func.create_csv_file(columns=['Gene name', 'Start position', 'End position', 'Strand', 'Error'],
    # #                        data=data_file_errors,
    # #                        csv_file_name='gene_exceptions.csv')
    #
    # print('\n\n#################### PART B ############################\n\n')
    #
    # # Index(['Entry', 'Entry name', 'Status', 'Protein names', 'Gene names',
    # #        'Length', 'Transmembrane', 'Sequence', 'name', 'locus'],
    # #       dtype='object')
    #
    unifile = 'uniprot_file.xlsx'
    part_b = GeneticDataGenerator(unifile=unifile, answers_file='Answer part B')

    part_b.uni_df.rename({'Gene names  (primary )': 'name',
                          'Gene names  (ordered locus )': 'locus'}, axis=1, inplace=True)
    #
    print(part_b.uni_df.columns)
    print('Question 1\n')
    same_prot_df, gb_only_df, uni_only_df, duplicates = func.compare_files_data(first_df=part_a.gb_df,
                                                                                second_df=part_b.uni_df,
                                                                                col1='locus_tag',
                                                                                col2='locus')
    #
    # print('Question 2\n')
    # trans_df = func.create_transmembrane_df(part_b.uni_df)
    #
    # trans_len = parts_functions.trans_len
    # hidro_prec = parts_functions.hidro_prec
    """
    The internal environment of the structure is Hydrophobic, 
    and the proteins in the Trans-membranous move through this structure
    and hence the expectation that the sequences 
    that create these proteins will be Hydrophobic or will have a high percentage of  
    hydrophobic amino acids .This assumption is similar to the results we received in the histogram
    """
    # trans_len_stats = parts_functions.calc_list_stats(trans_len)
    # print('Transmembrane Sequences Length Stats: \n', trans_len_stats)
    # hidro_prec_stat = parts_functions.calc_list_stats(hidro_prec)
    # print('Transmembrane Sequences Hidrophobic Stats: \n', hidro_prec_stat)

    # titles = ['Transmembrane Sequences Length Distribution', 'Transmembrane Sequences Hidro-acids% Distribution']
    # val = [trans_len, hidro_prec]
    # x_l = ['Transmembrane Sequences Length', 'Hidrophobic(%)']
    # bins = [list(range(10, 40, 5)), list(range(0, 100, 10))]
    #
    # for t, v, x, b in zip(titles, val, x_l, bins):
    #     parts_functions.build_histograms(hist_title=t,
    #                                      hist_value=v,
    #                                      x_label=x,
    #                                      y_label='Number of Sequences',
    #                                      bins_num=b)
    #
    # print('Question 3:\n')

    # print('#################### PART C ############################')
    # print('Question 1\n')
    # january = 'CoronaJanuary2022.gb'
    # july = 'CoronaJuly2020.gb'
    #
    # # list of cols to drop
    # drop_jan = ['organism', 'mol_type', 'isolate', 'host', 'db_xref',
    #             'country', 'collection_date', 'ribosomal_slippage', 'codon_start']
    # drop_jul = ['organism', 'mol_type', 'isolate', 'host', 'db_xref',
    #             'country', 'collection_date', 'ribosomal_slippage', 'codon_start',
    #             'inference', 'function', 'gene_synonym']
    #
    # january = GeneticDataGenerator(genebank_file=january, answers_file='january',
    #                                cols_to_drop=drop_jan)
    # july = GeneticDataGenerator(genebank_file=july, answers_file='july',
    #                             cols_to_drop=drop_jul)

    # print(january.gb_df['gene'].value_counts())
    # print(july.gb_df['gene'].value_counts())

    # covid_synon = func.count_mutation_by_type()
    # print(covid_synon)

    # print('Question 2\n')
    # same_gene_df, first_only_df, second_only_df = func.compare_files_data(first_df=january.gb_df,
    #                                                                       second_df=july.gb_df,
    #                                                                       col1='gene', col2='gene')
    #
    # print(same_gene_df.head(), first_only_df.head(), second_only_df.head())


