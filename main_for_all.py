from data_generator import GeneticDataGenerator
import part_a_func as a_func

if __name__ == '__main__':

    print('*--------------- Part A ---------------*')
    print('Question 1')

    gb_file = 'BS168.gb'
    part_a = GeneticDataGenerator(genebank_file=gb_file, answers_file='Answers part A')
    part_a.init_data_()

    print(part_a.gb_df.columns)

    all_gene_dict = a_func.get_all_genes_type_and_amount(part_a)
    print(f'All genes and their count in file: {all_gene_dict}')

    print('\nQuestion 2')
    df_proteins, protein_info, non_protein_info, protein_lengths, non_protein_lengths,\
        all_genes_lengths = a_func.characterization_of_gene_lengths(part_a)

    print(f'Statistics of each group:')
    print(f'proteins: {protein_info}\n'
          f'non proteins: {non_protein_info}\n')

    print(f'The histograms are shown on screen')

    hist_title = ['Non Protein length', 'Protein length', 'All genes length']
    hist_val = [non_protein_lengths, protein_lengths, all_genes_lengths]
    colors = ['blue', 'green', 'yellow']

    for _title, val, color in zip(hist_title, hist_val, colors):
        a_func.build_histograms(histogram_title=_title,
                                histogram_value=val,
                                x_label='Length',
                                y_label='Number of genes',
                                graph_color=color)

    print('\nQuestion 3')

    df_prot, largest_5, smallest_5, full_gc_percent, avg_prot_gc_percent = a_func.calculate_gc_percentage_in_genes(part_a, df_proteins)

    print(f'GC percentage of the whole genome: {full_gc_percent}')
    print(f'GC average percentage in the proteins: {avg_prot_gc_percent}')

    print(f'The histogram is shown on screen')
    a_func.build_histograms(histogram_title='GC percentage proteins',
                            histogram_value=list(df_prot['gene gc%']),
                            x_label='GC percentage',
                            y_label='Number of proteins',
                            )

    print(f'5 gene with the largest GC%:\n\t {largest_5}')
    print(f'5 gene with the smallest GC%:\n\t {smallest_5}')


    #
    #
    # print('#################### PART B ############################')
    # print('#################### PART C ############################')

